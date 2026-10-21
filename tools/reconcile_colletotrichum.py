#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""reconcile_colletotrichum.py —— 端到端重算物种级对账表。

三块输入：
  左（文献侧）  PDF 文本层的 `_T` 图标签        ← labels.extract_T_labels
  右（NCBI 侧） NCBI Taxonomy 的 type material ← sources.TaxonomyAPI（联网，落盘缓存）
  裁定表       已人工确认的错值清单             ← 副表1（rule ①，最高优先）

用法：
    python tools/reconcile_colletotrichum.py --sim-pdf <全文.txt> \
        --errors <副表1.tsv> --limit 40 --out <目录> [--no-manual]

--no-manual 时不加载人工裁定表，用于测量**机械规则本身**的一致率。
"""
import argparse
import csv
import io
import os
import re
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
from typeaudit import labels as LB, reconcile as RC, sources as SRC   # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8")
    sys.stderr.reconfigure(encoding="utf-8")
except Exception:
    pass


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--sim-pdf", required=True, help="SiM 全文（PDF 文本层抽出的 .txt）")
    ap.add_argument("--errors", help="副表1（人工裁定表，rule ①）")
    ap.add_argument("--ref", help="参考用 副表3（对照，不参与判定）")
    ap.add_argument("--registry", help="Biocollections CSV（判机构类型）")
    ap.add_argument("--cache", default=".cache")
    ap.add_argument("--limit", type=int, default=40)
    ap.add_argument("--out", required=True)
    ap.add_argument("--no-manual", action="store_true")
    ap.add_argument("--budget", type=float, default=0,
                    help="时间预算（秒）；到点即停并把已完成部分落盘，便于多次调用累积")
    a = ap.parse_args(argv)

    txt = io.open(a.sim_pdf, encoding="utf-8", errors="replace").read()
    left = LB.labels_by_species(txt)
    print("   左侧（PDF 图标签）：{} 个种".format(len(left)))

    manual = set()
    if a.errors and not a.no_manual and os.path.exists(a.errors):
        # ⚠️ rule ① 只对**物种级**错值生效。
        #    实测教训：E11（alatae）的错值所在源是「GenBank/RefSeq 记录」= **记录级**，
        #    副表3 的 alatae 行物种级判定是 ✅ 同株双保藏，错的只是某条记录挂错了号。
        #    把记录级错值也当物种级，会把 alatae 误判成 🔴。
        RECORD_LEVEL = "GenBank/RefSeq 记录"
        for r in csv.DictReader(io.open(a.errors, encoding="utf-8-sig", newline=""), delimiter="\t"):
            if (r.get("错值所在源") or "").strip() == RECORD_LEVEL:
                continue
            for x in re.split(r"[;；]", r.get("涉及物种") or ""):
                if x.strip():
                    manual.add(x.strip())
    print("   人工裁定表（仅物种级）：{} 个种{}".format(len(manual), "（已禁用）" if a.no_manual else ""))

    inst = SRC.Biocollections(a.registry) if a.registry else None
    tax = SRC.TaxonomyAPI(cache_dir=a.cache, email=os.environ.get("NCBI_EMAIL"))

    ref = {}
    if a.ref and os.path.exists(a.ref):
        for r in csv.DictReader(io.open(a.ref, encoding="utf-8-sig", newline=""), delimiter="\t"):
            ref[(r["species"] or "").strip()] = r

    # 物种清单：左侧 ∪ 参考表（限量，跨分布抽样）
    species = sorted(ref) if ref else sorted(left)
    # ⚠️ 取样要**顺序**取（`species[:N]`）而不是等距取——这样 limit=40 的样本是 limit=100 的超集，
    #    配合落盘缓存即可**渐进累积**（本次超时了，下次接着跑，不必重来）。
    if a.limit:
        species = species[:a.limit]

    rows, agree, total_cmp = [], 0, [0]
    t_start = time.time()
    done = 0
    for sp in species:
        if a.budget and time.time() - t_start > a.budget:
            print('   ⏱ 预算用尽，已处理 {} / {}（缓存已落盘，可再次运行续跑）'.format(done, len(species)))
            break
        done += 1
        ep = LB.canon(sp.split()[-1].lower())      # ⚠️ 查表侧也要归一：PDF 写 axonopodi，
                                                   #    副表3 写 axonopodis，不归一就对不上
        lf = left.get(ep, "")
        try:
            _tid, tms = tax.type_material_of(sp, retmax=1)
        except Exception as e:
            tms = []
        rf = " | ".join(tms)
        v, why = RC.verdict(lf, rf, species=sp, manual_hits=manual, inst=inst)
        rv = (ref.get(sp, {}) or {}).get("verdict", "")
        rows.append([sp, lf, rf, v, why, rv])
        if rv and not a.no_manual:
            # 度量口径：**排除副表3 判 🔴 的行**——
            # 那些是人工裁定的结果（rule ① 的输入），机械规则无从得知，
            # 把它们算进分母会系统性低估机械规则的一致率。
            def cls(x):
                if "有确认错值" in x or x.startswith("🔴"): return "err"
                if "一致" in x and "标本" not in x and "同株" not in x: return "ok"
                if "仅" in x: return "one"
                return "other"
            if cls(rv) != "err":
                total_cmp[0] += 1
                agree += cls(v) == cls(rv)

    os.makedirs(a.out, exist_ok=True)
    p = os.path.join(a.out, "重算_物种级对账.tsv")
    with io.open(p, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["物种", "左侧(SiM 图标签)", "右侧(NCBI Taxonomy)", "重算判定", "依据", "副表3 原判定"])
        w.writerows(rows)
    print("   写出 {}（{} 行）".format(p, len(rows)))
    print()
    print("   与副表3 的判定类别一致：{} / {}（{:.0%}）".format(
        agree, total_cmp[0], agree / total_cmp[0] if total_cmp[0] else 0))
    print("   （分母 = 副表3 中**非 🔴** 的行数；🔴 行属人工裁定，机械规则无从得知）")
    print()
    print("   样例：")
    for r in rows[:10]:
        print("      {:<32}{:<18}{:<18}{:<16}{}".format(r[0][:30], r[1][:16], r[2][:16], r[3][:14], r[5][:16]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
