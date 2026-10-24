#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""从 MBList 导出（mycobank_authors.tsv）生成两张查找表。

为什么要它：mb_class.tsv 有 70 MB、mb_current.tsv 有 24 MB，远超版本库的合理体积，
所以不入库；但守卫（tools/collect_findings.py 的 reclassify_aliases）依赖它们，
必须能从原始导出重建。

用法：
    python tools/build_mb_maps.py --src 交付_三属对照_20261022/08_MBList增强/mycobank_authors.tsv

产出：
    data/mb_current.tsv    taxon -> current name
    data/mb_class.tsv      taxon -> classification（MycoBank 的分类串，逗号分隔）
"""
import argparse, io, os


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", required=True)
    ap.add_argument("--outdir", default="data")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    p_cur = os.path.join(a.outdir, "mb_current.tsv")
    p_cls = os.path.join(a.outdir, "mb_class.tsv")
    seen = set()
    n = 0
    with io.open(a.src, encoding="utf-8", errors="replace") as f, \
         io.open(p_cur, "w", encoding="utf-8", newline="\n") as fo, \
         io.open(p_cls, "w", encoding="utf-8", newline="\n") as fc:
        h = f.readline().rstrip("\n").split("\t")
        ix = {k: h.index(k) for k in ("Taxon name", "Current name", "Classification")}
        fo.write("taxon\tcurrent\n")
        fc.write("taxon\tclassification\n")
        for line in f:
            q = line.rstrip("\n").split("\t")
            if len(q) < len(h):
                continue
            nm = q[ix["Taxon name"]].strip()
            if not nm or nm.lower() in seen:
                continue
            seen.add(nm.lower())
            n += 1
            fo.write("{}\t{}\n".format(nm, q[ix["Current name"]].strip()))
            fc.write("{}\t{}\n".format(nm, q[ix["Classification"]].strip()))
    print("已生成 {:,} 条 → {} / {}".format(n, p_cur, p_cls))


if __name__ == "__main__":
    main()
