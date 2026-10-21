# -*- coding: utf-8 -*-
"""全自动审计管道：一条命令从「属名」跑到「报告」。

    python -m typeaudit.audit --genus Colletotrichum --corpus .cache --out out
    python -m typeaudit.audit --genus Colletotrichum --fetch --email you@x.org --out out

六个阶段：
  1 acquire   取得语料（本地目录，或联网抓 PMC 全文）
  2 gate      C1 可行性门（三根轴）
  3 extract   逐物种抽类型指定段 / 等价号链 / 标本号
  4 rules     J-01…J-21 规则库（离线可跑的那部分）
  5 metrics   按**声明式口径**算指标（口径随报告一起输出，数字自解释）
  6 report    落盘 species TSV / findings TSV / metrics JSON + 控制台摘要

口径是**声明式常量**（见 CONVENTIONS）。改口径必须同时改这里与报告，
这样「同一个数字」在不同时间的含义不会漂移——本项目曾因此漏改过 3 处。
"""
import argparse
import csv
import io
import json
import os
import re
import sys
from collections import Counter, defaultdict

from . import designations as DG
from . import gate as G
from . import rules as R
from .specimen import candidates as specimen_candidates, type_sentence as specimen_type_sentence
from .textutil import STOPWORDS, blocks, segment, species_units, strip_tags, units

# ==================== 口径（声明式；报告里会原样带出） ====================
CONVENTIONS = {
    "block_variant": "B_paren",
    "block_variant_note": "种加词后接【大写字母或左括号】；仅要求大写字母会漏掉 "
                          "`Xylaria phyllocharis (Mont.) Fr.` 这类重组合名",
    "label_set": "full(%d)" % len(DG.LABEL_VARIANTS),
    "label_set_note": "designations.LABEL_VARIANTS 全变体表；用窄集会低估含类型指定块数",
    "chain_counting": "all",
    "chain_counting_note": "统计块内**全部** ex-type 链；只取首链会低估多号链数",
    "species_key": "epithet_lower",
}


# ==================== 1. acquire ====================
def acquire(corpus_dir=None, genus=None, fetch=False, cache_dir=".cache",
            email=None, max_papers=12, pause=1.3):
    """返回 [(来源标识, XML/文本内容)]。fetch=True 时联网抓取（需 biopython）。"""
    items = []
    if corpus_dir and os.path.isdir(corpus_dir):
        for f in sorted(os.listdir(corpus_dir)):
            if not f.lower().endswith((".xml", ".txt")):
                continue
            # ⚠️ 必须过滤：缓存目录里还躺着 esearch 的结果文件
            #    （search_*.xml / count_*.xml / s3_*.xml / s4_*.xml）。
            #    它们不是论文，收进来会凭空产生一批"样本无效"噪音（实测 8 条）。
            if re.match(r"^(search_|count_|s[34]_|gb_|gbT_|gbs_)", f):
                continue
            p = os.path.join(corpus_dir, f)
            raw = io.open(p, encoding="utf-8", errors="replace").read()
            if f.lower().endswith(".xml") and not re.search(r"<article|<pmc-articleset", raw):
                continue                      # 不是文章 XML，跳过
            items.append((f, raw))
    if fetch and genus:
        for pid, _fg in G.batch_gate(genus, cache_dir=cache_dir, email=email,
                                     maxpapers=max_papers, pause=pause):
            p = os.path.join(cache_dir, "pmc_{}.xml".format(pid))
            if os.path.exists(p):
                items.append(("pmc_{}".format(pid), io.open(p, encoding="utf-8", errors="replace").read()))
    return items


# ==================== 2+3. gate & extract ====================
def extract_one(source, raw, genus):
    """对单篇：跑门 + 抽物种级记录。返回 (gate 结果, [记录])。"""
    text = strip_tags(raw)
    fg = G.full_gate(text, genus)
    # ⚠️ 切分路径：**结构化语料优先用 XML 单元**（文档序的 <p>/<caption>/<td>/<list-item>），
    #    压平文本上的正则切块只在纯 txt 语料里用。
    #    实测差异：单元路径下 J-07 的一号两主候选从 12 条降到 7 条，且能抓到真实的一号两主；
    #    文本路径的候选大多来自「块跨表格」——同一个号落进多个物种块。
    bl = segment(raw, genus)
    path = "units" if re.search(r"<p[^>]*>|<title[^>]*>", raw or "") else "blocks"
    recs = []
    for ep in sorted(bl):
        blk = bl[ep]
        if not any(re.search(re.escape(l), blk, re.I) for l in DG.LABEL_VARIANTS):
            continue
        chains_all = DG.parse_chains(blk)                     # 口径：全部链
        first = chains_all[0][1] if chains_all else None
        ho = DG.RE_HOLO.search(blk)
        ts = specimen_type_sentence(blk)
        recs.append({
            "source": source, "species": ep,
            "chain": first or [], "chain_len_multi": sum(1 for _, t in chains_all if len(t) >= 2),
            "n_chains": len(chains_all),
            "all_numbers": sorted({x for _, t in chains_all for x in t}),
            "holotype": (ho.group(1).strip() if ho else ""),
            "vouchers": [c[0] + c[1] for c in specimen_candidates(ts)] if ts else [],
            "basis": fg["basis"],
        })
    fg["segmentation"] = path
    return fg, recs


# ==================== 4. rules ====================
def run_rules(items, per_source, genus, all_text):
    """跑离线可跑的 J 规则。返回 findings 列表 + 覆盖报告。"""
    out, ran = [], []

    # J-07：同侧内部，同一号指给两个物种
    for src, recs in per_source.items():
        m = defaultdict(set)
        for r in recs:
            # ⚠️ 只喂**链内声明的号**（ex-type 等价链），不要喂块内全部号——
            #    块会跨表格，导致同一个号出现在多个物种块里，产生大量假阳性。
            for n in r.get("chain", []):
                m[n].add(r["species"])
        f = R.j07_one_number_two_species(m)
        for x in f:
            x["source"] = src
        out += f
        ran.append("J-07")

    # J-08：同批位数孤例（按「来源 × 机构前缀」分组）
    recs_batch = []
    for src, recs in per_source.items():
        for r in recs:
            for n in r["all_numbers"]:
                code, num = R.split_code(n)
                if code and num:
                    recs_batch.append({"group": "{}|{}".format(src, code), "num": num})
    out += R.j08_batch_regularity(recs_batch, group_key="group")
    ran.append("J-08")

    # J-13：跨来源数字级近失配（同物种、两只来源）
    by_species = defaultdict(dict)
    for src, recs in per_source.items():
        for r in recs:
            if r["all_numbers"]:
                by_species[r["species"]][src] = r["all_numbers"]
    n13 = 0
    for sp, d in by_species.items():
        srcs = sorted(d)
        for i in range(len(srcs)):
            for j in range(i + 1, len(srcs)):
                for x in R.j13_near_miss(d[srcs[i]], d[srcs[j]]):
                    x["subject"] = "{} · {}".format(sp, x["subject"])
                    x["evidence"] = "{} vs {}".format(srcs[i], srcs[j])
                    out.append(x); n13 += 1
    ran.append("J-13")

    # J-01：跨来源对账（四态）
    n01 = 0
    for sp, d in by_species.items():
        srcs = sorted(d)
        if len(srcs) < 2:
            continue
        for i in range(len(srcs)):
            for j in range(i + 1, len(srcs)):
                for x in R.j01_cross_source({sp: d[srcs[i]]}, {sp: d[srcs[j]]}):
                    x["evidence"] = "{} vs {}".format(srcs[i], srcs[j])
                    out.append(x); n01 += 1
    ran.append("J-01")

    # J-14：图标签 vs 正文
    n14 = 0
    for src, raw in items:
        fig_txt = " ".join(re.findall(r"<fig[^>]*>.*?</fig>", raw, re.S))
        f = R.j14_figure_vs_body(strip_tags(fig_txt), strip_tags(raw), genus)
        for x in f:
            x["source"] = src
        out += f; n14 += len(f)
    ran.append("J-14")

    # J-15：枚举范围（「类型材料字段枚举」vs「全部记录标识符枚举」）
    enum_declared = {n for recs in per_source.values() for r in recs for n in r["all_numbers"]}
    # J-15 的"较宽枚举"要**收窄到像保藏号的东西**，否则会退化成"正文里所有 code+数字"
    # （实测原文：不收窄会抽出 8356 个记号，含登录号与基因坐标，比较失去意义）。
    from .specimen import ACC_PREFIX, NOT_HERB
    enum_all = set()
    for _, raw in items:
        for code, num in re.findall(r"\b([A-Z][A-Za-z]{1,11})\s*[:：]?\s*(\d{3,9})\b", strip_tags(raw)):
            C = code.upper()
            if C in ACC_PREFIX or C in NOT_HERB or len(C) < 2:
                continue
            enum_all.add(C + num)
    scope = R.j15_enumeration_scope(enum_declared, enum_all,
                                    "类型材料字段里的号", "全文出现的号")
    ran.append("J-15")
    ran.append("J-06")     # J-06（原文模式指定段）由第 3 阶段 extract 实现

    # J-20：本地 .gb 文件里的 Type status（若有）
    gb_hits = []
    for d in (".cache", "data", "."):
        if os.path.isdir(d):
            for f in os.listdir(d):
                if f.endswith(".gb"):
                    gb_hits += R.j20_note_type_status(
                        io.open(os.path.join(d, f), encoding="utf-8", errors="replace").read())
    if gb_hits:
        ran.append("J-20")
    out += gb_hits

    return out, {"ran": sorted(set(ran)), "scope": scope}


# ==================== 5. metrics（按口径） ====================
def metrics(items, gates, per_source, genus):
    all_recs = [r for recs in per_source.values() for r in recs]
    chains_multi = sum(r["chain_len_multi"] for r in all_recs)
    return {
        "genus": genus,
        "conventions": CONVENTIONS,
        "n_sources": len(items),
        "n_blocks": sum(g[1]["n_blocks"] for g in gates),
        "n_blocks_with_label": sum(g[1]["n_blocks_with_label"] for g in gates),
        "n_species_records": len(all_recs),
        "n_records_with_chain": sum(1 for r in all_recs if r["chain"]),
        "n_multi_number_chains": chains_multi,          # 口径：全部链
        "n_distinct_numbers": len({n for r in all_recs for n in r["all_numbers"]}),
        "basis": dict(Counter(g[1]["basis"] for g in gates)),
        "verdicts": dict(Counter(g[1]["verdict"].split("（")[0].split(" ·")[0] for g in gates)),
    }


# ==================== 6. report ====================
def _w(path, header, rows):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with io.open(path, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(header); w.writerows(rows)


def run(genus, corpus_dir=None, fetch=False, out_dir="out", cache_dir=".cache",
        email=None, max_papers=12):
    items = acquire(corpus_dir, genus, fetch, cache_dir, email, max_papers)
    if not items:
        print("没有语料：给出 --corpus <目录> 或加 --fetch 联网抓取。")
        return 1
    print("=" * 96); print("【1-3 采集 / 门 / 抽取】属 = {}，语料 {} 篇".format(genus, len(items))); print("=" * 96)
    gates, per_source = [], {}
    for src, raw in items:
        fg, recs = extract_one(src, raw, genus)
        nb = len(blocks(strip_tags(raw), genus, "B"))
        nlab = sum(1 for e, b in blocks(strip_tags(raw), genus, "B").items()
                   if any(re.search(re.escape(l), b, re.I) for l in DG.LABEL_VARIANTS))
        gates.append((src, dict(fg, n_blocks=nb, n_blocks_with_label=nlab)))
        per_source[src] = recs
        print("   {:<18}{:<9}块{:>4} 含标签{:>4} 记录{:>4}  {}".format(
            src[:16], fg["basis"], nb, nlab, len(recs), fg["verdict"][:34]))

    all_text = " ".join(strip_tags(r) for _, r in items)
    print()
    print("=" * 96); print("【4 J 规则库】"); print("=" * 96)
    found, cov = run_rules(items, per_source, genus, all_text)
    print("   已跑规则：{}".format(" ".join(cov["ran"])))
    print("   未跑规则（需外部源）：{}".format(
        " ".join(r["id"] for r in R.all_rules() if r["id"] not in cov["ran"])))
    print("   ⚙️ J-15 枚举范围：类型材料字段 {} 个号 vs 全文 {} 个号 ⇒ **{} 个号从未被检查**".format(
        cov["scope"]["n_a"], cov["scope"]["n_b"], cov["scope"]["n_only_in_b"]))
    kept = [f for f in found if f["severity"] in ("candidate", "evidence")]
    print("   候选合计 {} 条（ok/untestable 共 {} 条不入表）".format(len(kept), len(found) - len(kept)))
    for k, v in Counter(f["rule"] for f in kept).most_common():
        print("      {}：{} 条候选".format(k, v))

    m = metrics(items, gates, per_source, genus)
    print()
    print("=" * 96); print("【5 口径化指标（数字自解释）】"); print("=" * 96)
    print("   口径：切块 {} ；标签集 {} ；链计数 {}".format(
        CONVENTIONS["block_variant"], CONVENTIONS["label_set"], CONVENTIONS["chain_counting"]))
    for k in ["n_sources", "n_blocks", "n_blocks_with_label", "n_species_records",
              "n_records_with_chain", "n_multi_number_chains", "n_distinct_numbers"]:
        print("   {:<26}{}".format(k, m[k]))
    print("   基质分布 {}".format(m["basis"]))
    print("   判定分布 {}".format(m["verdicts"]))

    os.makedirs(out_dir, exist_ok=True)
    name = genus.lower()
    _w(os.path.join(out_dir, "{}_species.tsv".format(name)),
       ["source", "species", "basis", "chain", "n_chains", "multi", "all_numbers", "holotype", "vouchers"],
       [[r["source"], r["species"], r["basis"], " = ".join(r["chain"]), r["n_chains"],
         r["chain_len_multi"], "|".join(r["all_numbers"]), r["holotype"], "|".join(r["vouchers"])]
        for recs in per_source.values() for r in recs])
    _w(os.path.join(out_dir, "{}_findings.tsv".format(name)),
       ["rule", "source", "subject", "detail", "severity", "evidence"],
       [[f["rule"], f.get("source", ""), f["subject"], f["detail"], f["severity"], f["evidence"]]
        for f in found if f["severity"] in ("candidate", "evidence")])   # ok / untestable 不入表
    with io.open(os.path.join(out_dir, "{}_metrics.json".format(name)), "w", encoding="utf-8") as f:
        json.dump({"metrics": m, "rule_coverage": cov, "findings_n": len(found)}, f, ensure_ascii=False, indent=1)
    print()
    print("   报告落盘：{}/{}_species.tsv  {}_findings.tsv  {}_metrics.json".format(out_dir, name, name, name))
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description="全自动类型材料审计管道")
    ap.add_argument("--genus", required=True)
    ap.add_argument("--corpus", help="本地语料目录（.xml/.txt）")
    ap.add_argument("--fetch", action="store_true", help="联网抓 PMC（需 biopython + --email）")
    ap.add_argument("--out", default="out")
    ap.add_argument("--cache-dir", default=".cache")
    ap.add_argument("--email", default=os.environ.get("NCBI_EMAIL"))
    ap.add_argument("--max-papers", type=int, default=12)
    a = ap.parse_args(argv)
    return run(a.genus, a.corpus, a.fetch, a.out, a.cache_dir, a.email, a.max_papers)


if __name__ == "__main__":
    sys.exit(main())
