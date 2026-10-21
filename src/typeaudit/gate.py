# -*- coding: utf-8 -*-
"""C1 可行性门：换属之前先判「这个属能不能审」。

三根轴（每条都来自实测翻车）：
  一、切块规则：不要要求作者名以「大写字母+点」开头（见 textutil）
  二、类型基质：等价链法只适用于**培养物型**；标本型抽出率实测 0%
      实测 166 条目 / 78 属：culture 117 条 → 抽出率 63%；specimen 49 条 → 抽出率 0%
  三、样本真伪：检索命中标题 ≠ 论文处理该属。必须校验「属名在正文的出现次数」。
      实测 `"Bipolaris"[Title]` 命中一篇《Two New Chalcid Wasps…》——
      `Ophelimus bipolaris` 是**膜翅目小蜂**的种加词，正文里 Bipolaris 出现 0 次。

用法：
    python -m typeaudit.gate --file paper.xml --genus Diaporthe
    python -m typeaudit.gate --pmcid 10301816 --genus Diaporthe --cache-dir .cache
    python -m typeaudit.gate --genus Diaporthe --batch        # 联网抓取（需 biopython）
"""
import argparse
import io
import os
import re
import sys
import time

from .textutil import blocks, segment, strip_tags
from .designations import LABEL_VARIANTS, parse_chains

MIN_GENUS_MENTIONS = 10        # 属名在正文的出现次数下限（低于即判假样本）
VARIANTS = ["A", "B", "C", "D"]


def _has_label(blk):
    return any(re.search(re.escape(l), blk, re.I) for l in LABEL_VARIANTS)


def gate(text, genus, raw=None):
    """四变体各算 [块数, 含标签块数, 抽出 ex-type 数]。

    给了 raw（原始 XML）时按单元切分；否则按文本切块。
    """
    rows = {}
    for v in VARIANTS:
        bl = segment(raw, genus) if (raw and v == "B") else blocks(
            text, genus, v, stop=None if v != "D" else __import__(
                "typeaudit.textutil", fromlist=["STOPWORDS"]).STOPWORDS)
        want = [e for e in bl if _has_label(bl[e])]
        got = [e for e in want if parse_chains(bl[e])]
        rows[v] = [len(bl), len(want), len(got)]
    return rows


def type_basis(t):
    """判定该文的模式指定属于培养物型还是标本型。返回 (类, ex-type 数, holotype 数, specimen 数)。"""
    ex = len(re.findall(r"ex-type culture|ex-holotype|ex-type", t, re.I))
    holo = len(re.findall(r"\bholotype\b", t, re.I))
    spec = len(re.findall(r"\bspecimen\b", t, re.I))
    if ex >= 1:
        return "culture", ex, holo, spec
    if holo >= 1:
        return "specimen", ex, holo, spec
    return "unknown", ex, holo, spec


def genus_mentions(t, genus):
    return len(re.findall(r"\b" + re.escape(genus) + r"\b", t))


def full_gate(text, genus):
    """三根轴合成判定。"""
    v = gate(text, genus)
    basis, ex, holo, spec = type_basis(text)
    gm = genus_mentions(text, genus)
    nn = len(re.findall(r"\bsp\.\s*nov\.", text))
    if gm < MIN_GENUS_MENTIONS:
        verdict = "样本无效（属名正文出现 {} 次 < {}）疑标题同名异属，须剔除".format(gm, MIN_GENUS_MENTIONS)
    elif basis == "specimen":
        verdict = "不可审 · 等价链法不适用（标本型：ex-type culture = 0）"
    elif v["B"][1] == 0:
        verdict = "不可审 · 切块失败（变体 B 也无带标签块）"
    elif v["B"][2] == 0:
        verdict = "部分可审 · 切块成功但抽不出 ex-type"
    else:
        verdict = "可审 · B 块 {} 含标签 {} 抽出 {}".format(*v["B"])
    return {"variants": v, "basis": basis, "basis_stats": [ex, holo, spec],
            "genus_mentions": gm, "n_spnov": nn, "verdict": verdict}


# ---------------- 联网抓取（惰性导入 biopython） ----------------
def _entrez(email=None):
    try:
        from Bio import Entrez
    except ImportError:
        raise SystemExit(
            "联网抓取需要 biopython（pip install biopython）。\n"
            "若已有本地文件，可用 --file 直接分析，无需联网。")
    Entrez.email = email or os.environ.get("NCBI_EMAIL")
    if not Entrez.email:
        raise SystemExit("NCBI 要求提供邮箱：用 --email 或设环境变量 NCBI_EMAIL。")
    return Entrez


def _rd(h):
    d = h.read()
    return d.decode("utf-8", "replace") if isinstance(d, bytes) else d


def _cache(cache_dir, name, make):
    os.makedirs(cache_dir, exist_ok=True)
    p = os.path.join(cache_dir, name)
    if os.path.exists(p):
        return io.open(p, encoding="utf-8", errors="replace").read()
    s = make()
    io.open(p, "w", encoding="utf-8", errors="replace").write(s)
    return s


def fetch_pmc(pmcid, cache_dir=".cache", email=None, pause=1.3):
    """抓 PMC 全文 XML（走缓存）。"""
    def make():
        Entrez = _entrez(email)
        time.sleep(pause)
        return _rd(Entrez.efetch(db="pmc", id=str(pmcid), rettype="xml"))
    return _cache(cache_dir, "pmc_{}.xml".format(pmcid), make)


def batch_gate(genus, cache_dir=".cache", email=None, maxpapers=4, pause=1.3):
    """按四段检索阶梯抓该属的 PMC 全文并过门。返回 [(pmcid, full_gate 结果)]。"""
    LADDER = [
        '"{g}"[Title] AND (holotype[Text Word] OR Typus[Text Word]) AND "sp. nov."[Text Word]',
        '"{g}"[Title] AND (holotype[Text Word] OR Typus[Text Word])',
        '"{g}"[Title] AND ("new species"[Title/Abstract] OR "sp. nov."[Title/Abstract] OR taxonomy[Title/Abstract])',
        '"{g}"[Title/Abstract] AND (holotype[Text Word] OR Typus[Text Word])',
    ]
    Entrez = _entrez(email)
    # ⚠️ 阶梯语义：**逐级累积**，不要"命中即停"。
    #    实测 Colletotrichum：最严的 R1 命中 3 条，而 R2 命中 75 条、R4 命中 107 条。
    #    "命中即停"会把 75 篇里的 3 篇当全部（欠采样 96%），而它看起来完全正常。
    ids, seen = [], set()
    for tpl in LADDER:
        def mk(t=tpl):
            time.sleep(pause)
            return _rd(Entrez.esearch(db="pmc", term=t.format(g=genus), retmax=maxpapers, sort="relevance"))
        for i in re.findall(r"<Id>(\d+)</Id>",
                            _cache(cache_dir, "search_{}_{}.xml".format(genus, abs(hash(tpl)) % 99999), mk)):
            if i not in seen:
                seen.add(i); ids.append(i)
        if len(ids) >= maxpapers:
            break
    ids = ids[:maxpapers]
    out = []
    for pid in ids:
        t = strip_tags(fetch_pmc(pid, cache_dir, email, pause))
        out.append((pid, full_gate(t, genus)))
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description="C1 可行性门：判断某属能否用类型材料审计方法")
    ap.add_argument("--genus", required=True, help="属名，如 Diaporthe")
    ap.add_argument("--file", help="本地全文文件（.xml 或 .txt）；给了就不联网")
    ap.add_argument("--pmcid", help="PMC ID；给了就联网抓这一篇")
    ap.add_argument("--batch", action="store_true", help="联网按检索阶梯抓该属若干篇")
    ap.add_argument("--cache-dir", default=".cache", help="缓存目录（默认 .cache）")
    ap.add_argument("--email", default=None, help="NCBI 要求；也可设环境变量 NCBI_EMAIL")
    a = ap.parse_args(argv)

    if a.file:
        text = strip_tags(io.open(a.file, encoding="utf-8", errors="replace").read())
        fg = full_gate(text, a.genus)
        print("文件：{}".format(a.file))
        print("属名正文出现 {} 次 ；sp.nov. {} 次 ；基质 {}".format(
            fg["genus_mentions"], fg["n_spnov"], fg["basis"]))
        for v in VARIANTS:
            print("   变体 {}：块 {:>5} 含标签 {:>5} 抽出 {:>5}".format(v, *fg["variants"][v]))
        print("判定：{}".format(fg["verdict"]))
        return 0
    if a.pmcid:
        text = strip_tags(fetch_pmc(a.pmcid, a.cache_dir, a.email))
        fg = full_gate(text, a.genus)
        print("PMC{}：{}".format(a.pmcid, fg["verdict"]))
        return 0
    if a.batch:
        for pid, fg in batch_gate(a.genus, a.cache_dir, a.email):
            print("PMC{:<10}{}".format(pid, fg["verdict"]))
        return 0
    ap.error("需要 --file / --pmcid / --batch 之一")


if __name__ == "__main__":
    sys.exit(main())
