# -*- coding: utf-8 -*-
r"""反向审计：从【保藏目录】反推 GenBank —— 目录说它是模式菌株，序列库里有声明吗？

    python tools\reverse_audit.py --genus Colletotrichum [--limit N] [--out out_reverse]

【为什么要做反向】
  现有流水线是单向的：GenBank 记录 → 抽号 → 核目录。它只能发现「记录声明了、目录里没有」。
  反过来问是**另一种错误类型**，而且不依赖 GenBank 先声明什么：
      「目录里标为模式菌株的号」→ GenBank 里有没有任何记录声明它？
  没有 ⇒ 该种的**模式材料序列数据缺失**（从未测序，或测了但没写号）。

【为什么不能只看命中数 —— ATM 陷阱，本机实测】
  NCBI 的自动词映射（ATM）会把「散文式」查询映射到 Taxonomy。而 CBS 类型号
  恰恰登记在 Taxonomy 里，于是：
      "CBS 112996"  → 20,561 条，抓回来一核，**一个都不含该号**（全是 C. acutatum 的记录）
      "CBS 133.44"  → 19,130 条，全是 C. godetiae
      "CBS 125378"  → 14 条，抓回来 3/3 确实含 CBS:125378     ← 这个是真的
  ⇒ **命中数不是证据。** 必须 ① 用字段限定绕开 ATM，② 把记录抓回来实证。

【实测选定的检索式】`"<CODE> <号>"[Text Word]`

  为什么不是 `[Strain]`（本机实测抓到 2 个假阴性）：
      NCBI 的 `[Strain]` 索引**不覆盖 /culture_collection**。只写在 cc、没写进
      /strain 的号会被漏掉。实测：
          C. fructicola CBS 130416 → [Strain]=0，[Text Word]=1（真记录 JQ899290.1）
          C. acericola  CBS 149616 → [Strain]=0，[Text Word]=1（真记录 NR_191204.1）
      一开始用 [Strain] 做审计，这 2 条被误判成「GenBank 无声明」。

  为什么不是 `[All Fields]`：会被 ATM 劫持（见上）。

  `[Text Word]` 两个条件都满足：
      125378 → 14 条，3/3 真含号 ✅      （[Strain] 只有 5 条）
      112996 → 11 条，3/3 真含号 ✅      （[All Fields] 是 20,561 条，全假）
      133.44 → 10 条，3/3 真含号 ✅      （[All Fields] 是 19,130 条，全假）
      999999 →  0 条（负对照）✅

【诚实边界】
  · 这是**字段限定**检索，可能漏掉「号只写在别处」的记录。
  · 结论说的是「**GenBank 里没有记录声明该号**」，不等于该菌株没有序列数据。
  · 「目录标为模式菌株」以 CBS 目录的 status 文本为准（可能本身有误）。
"""
import argparse
import hashlib
import io
import json
import os
import re
import sys
import time

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(REPO, "src"))

TYPE_WORDS = (r"(?i)\b(holotype|isotype|neotype|epitype|paratype|ex[- ]?type|"
              r"type of|ex holotype|authentic strain)\b")


def binomial(status, genus):
    m = re.search(re.escape(genus) + r"\s+([a-z][a-z\-]+)", status or "", re.I)
    if not m:
        return None
    # 保留原文的属名大小写
    g = re.search(re.escape(genus), status, re.I).group(0)
    return "{} {}".format(g.capitalize(), m.group(1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genus", default="Colletotrichum")
    ap.add_argument("--out", default="out_reverse")
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--email", default=os.environ.get("NCBI_EMAIL"))
    ap.add_argument("--check", action="store_true",
                    help="对已跑过的结果做一次汇总（不联网）")
    a = ap.parse_args()

    od = a.out if os.path.isabs(a.out) else os.path.join(REPO, a.out)
    os.makedirs(od, exist_ok=True)
    out_tsv = os.path.join(od, "反向审计.tsv")

    if a.check:
        if not os.path.isfile(out_tsv):
            print("还没有结果：{}".format(out_tsv))
            return 2
        rows = [l.rstrip("\n").split("\t") for l in io.open(out_tsv, encoding="utf-8-sig")
                if l.strip()][1:]
        import collections
        c = collections.Counter(r[5] for r in rows)
        print("反向审计汇总：{}".format(dict(c)))
        for r in rows:
            if r[5] != "已在 GenBank 声明":
                print("   {} ｜ {} ｜ {}".format(r[0][:32], r[1], r[2][:40]))
        return 0

    # ---- 1) 从目录取「带类型地位」的声明 ----
    from typeaudit import sources as S
    srcs = S.load(os.path.join(REPO, "config", "sources.json"))
    cbs = srcs.get("cbs")
    if cbs is None:
        print("❌ 没有 CBS 目录源")
        return 2
    items = []
    for r in getattr(cbs, "rows", []):
        st = (r.get("status") or "").strip()
        if not re.search(re.escape(a.genus), st, re.I):
            continue
        if not re.search(TYPE_WORDS, st):
            continue
        nm = binomial(st, a.genus)
        if not nm:
            continue
        items.append({"species": nm, "cbs_number": (r.get("cbs_number") or "").strip(),
                      "status": st})
    # 去重：同 (种, 号) 只留一次
    seen, uniq = set(), []
    for it in items:
        k = (it["species"], it["cbs_number"])
        if k in seen:
            continue
        seen.add(k)
        uniq.append(it)
    if a.limit:
        uniq = uniq[:a.limit]
    print("目录里 {} 的带类型地位声明：{} 条（去重后）".format(a.genus, len(uniq)))

    # ---- 2) 逐号回查 GenBank ----
    from Bio import Entrez
    if not (a.email or os.environ.get("NCBI_EMAIL")):
        print("❌ 需要 NCBI 邮箱：--email 或环境变量 NCBI_EMAIL")
        return 2
    Entrez.email = a.email or os.environ.get("NCBI_EMAIL")
    Entrez.tool = "typeaudit-reverse"

    cache_dir = os.path.join(od, "queries")
    os.makedirs(cache_dir, exist_ok=True)

    def cached(key, fn):
        p = os.path.join(cache_dir, re.sub(r"\W+", "_", key)[:80] + ".json")
        if os.path.exists(p):
            try:
                return json.load(io.open(p, encoding="utf-8"))
            except Exception:
                pass
        v = fn()
        try:
            json.dump(v, io.open(p, "w", encoding="utf-8"), ensure_ascii=False)
        except Exception:
            pass
        return v

    def esearch(term, retmax=5):
        for i in range(3):
            try:
                time.sleep(1.35)
                h = Entrez.esearch(db="nucleotide", term=term, retmax=retmax)
                x = h.read()
                x = x.decode("utf-8", "replace") if isinstance(x, bytes) else x
                c = re.search(r"<Count>(\d+)</Count>", x)
                return {"count": int(c.group(1)) if c else 0,
                        "ids": re.findall(r"<Id>(\d+)</Id>", x)}
            except Exception as e:
                if i == 2:
                    return {"count": None, "ids": [], "error": type(e).__name__}
                time.sleep(3)

    def efetch_has(ids, num, n=3):
        """抓回前 n 条，检查号是否真的出现在记录文本里。"""
        if not ids:
            return {"checked": 0, "hit": 0}
        for i in range(3):
            try:
                time.sleep(1.35)
                h = Entrez.efetch(db="nucleotide", id=",".join(ids[:n]),
                                  rettype="gb", retmode="text",
                                  seq_start=1, seq_stop=5)
                d = h.read()
                d = d.decode("utf-8", "replace") if isinstance(d, bytes) else d
                recs = [r for r in re.split(r"(?=^LOCUS)", d, flags=re.M)
                        if r.strip().startswith("LOCUS")]
                return {"checked": len(recs),
                        "hit": sum(1 for r in recs if num in r)}
            except Exception:
                time.sleep(3)
        return {"checked": 0, "hit": 0}

    out_rows = []
    for i, it in enumerate(uniq, 1):
        num = it["cbs_number"]
        bare = re.sub(r"^CBS\s*", "", num).strip()
        q = '"CBS {}"'.format(bare)
        print("   [{:>3}/{}] {:<34} {}".format(i, len(uniq), it["species"][:34], num))
        # ⚠️ 用 [Text Word]：它覆盖 /culture_collection、/strain、/note、Title，
        #    且【不触发】ATM。原来用 [Strain] 会漏掉只写在 cc 里的号（实测 2 个假阴性）。
        sw = cached("textw_" + bare, lambda: esearch(q + "[Text Word]"))
        cnt = sw.get("count") or 0
        ids = sw.get("ids") or []
        # ⚠️ 核验结果的缓存键【必须含 ID 列表】。只按号做键的话：
        #    第一轮用 [Strain] 时命中 0 条、核验结果是 {checked:0,hit:0}，会被缓存下来；
        #    换了检索式后明明命中了，却仍读到那份陈旧结果 ⇒ 真记录被判成「未实证」。
        #    本机实测踩到（CBS 130416 / CBS 149616 两条假阴性就来自这里）。
        _idkey = hashlib.sha1(",".join(ids).encode()).hexdigest()[:12] if ids else "empty"
        v = cached("verify_{}_{}".format(bare, _idkey), lambda: efetch_has(ids, bare))
        if cnt == 0:
            # ---- CBS 号查不到声明时，先别急着判「无声明」----
            # 一管菌常同时存于多家机构，作者只填手边那个号。实测 6 条「无声明」里
            # 3 条是这种假阳性（MAFF 240517 / CPC 28612 / CPC 30245）。
            # 分两级证据：
            #   强证据：CBS 记录自己登记了别名（other_coll），且那个别名在 GenBank 被声明
            #   弱证据：CBS 没登记别名，但该种的类型材料记录在 GenBank 声明了【别的号】
            _aliases = []
            try:
                _mp = srcs.get("cbsmap")
                if _mp is not None and hasattr(_mp, "numbers_of"):
                    _aliases = (_mp.numbers_of(num) or {}).get("aliases") or []
            except Exception:
                _aliases = []
            _hit = []
            for _a in _aliases:
                _disp = _a["display"]
                _sw = cached("alt_" + re.sub(r"\W+", "_", _disp),
                             lambda _d=_disp: esearch('"{}"[Text Word]'.format(_d)))
                _n = _sw.get("count") or 0
                if _n > 0:
                    _hit.append((_disp, _n))
            if _hit:
                concl = "经他机构号声明"
                note = "CBS 号未声明；该菌以 {} 在 GenBank 声明（命中 {:d} 条）".format(
                    " / ".join(x[0] for x in _hit), _hit[0][1])
            else:
                # 弱证据：该种在 GenBank 的类型材料记录声明了哪些号
                _ccs = []
                try:
                    _gb = srcs.get("genbank")
                    if _gb is not None and hasattr(_gb, "type_material_cc_of"):
                        _ccs = _gb.type_material_cc_of(it["species"])
                except Exception:
                    _ccs = []
                _other = sorted({c for c in _ccs
                                 if c and not re.match(r"^\s*CBS", c, re.I)})
                if _other:
                    concl = "该种类型材料以他号声明"
                    note = ("CBS 号在 GenBank 无声明，且 CBS 记录未登记别名；"
                            "但该种的类型材料记录声明了 {}（CBS 号缺失，可能是同菌异号）"
                            .format(" / ".join(_other[:4])))
                else:
                    concl = "GenBank 无声明"
                    note = ("[Text Word] 字段检索命中 0；已查别名 {:d} 个均未命中；"
                            "该种类型材料记录未声明任何保藏号").format(len(_aliases))
        elif v.get("hit"):
            concl, note = "已在 GenBank 声明", "抓回 {:d} 条中 {:d} 条真含该号".format(
                v.get("checked", 0), v.get("hit", 0))
        else:
            concl, note = "命中但未实证", "命中 {:d} 条，抓回的都【不含】该号（可能 ATM 干扰）".format(cnt)
        out_rows.append([it["species"], num, it["status"][:80], str(cnt),
                         "{}/{}".format(v.get("hit", 0), v.get("checked", 0)),
                         concl, note])

    with io.open(out_tsv, "w", encoding="utf-8-sig", newline="") as f:
        f.write("种名\t目录号\t目录说的类型地位\t检索命中\t实证核验\t结论\t说明\n")
        for r in out_rows:
            f.write("\t".join(r) + "\n")

    import collections
    c = collections.Counter(r[5] for r in out_rows)
    print()
    print("=" * 78)
    print("反向审计完成：{} 条".format(len(out_rows)))
    for k, v in c.most_common():
        print("   {:<18} {}".format(k, v))
    print("   结果表：{}".format(out_tsv))
    print("=" * 78)
    print()
    print("缺声明的（目录说它是模式菌株，但 GenBank 里没记录声明该号）：")
    for r in out_rows:
        if r[5] == "GenBank 无声明":
            print("   {:<36} {}".format(r[0][:36], r[1]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
