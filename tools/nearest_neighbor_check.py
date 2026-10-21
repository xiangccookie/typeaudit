#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
nearest_neighbor_check.py —— 按【株】做最近邻检验（并名/错配候选），不建树。

【为什么不是全对建树】
    现有 type_synonym_verdict.py 是按【对】跑 8 基因树 —— 只能比两个号，
    全属两两组合是 O(n²)，不可行。
    本工具改成按【株】做 best-close-match 检验，复杂度 O(n)：
        对每株 X（声明属于种 G）：
            d_in  = X 到 G 内【其他株】的最大距离
            d_out = X 到任意【非 G】株的最小距离
            若 d_out < d_in ⇒ X 落到别的种里 ⇒ 并名 / 错配候选
    只对挑出来的少数候选，才值得跑 8 基因建树确认（把最贵的工具用在刀刃上）。

【结论强度不对称 —— 必须写进任何引用本结果的报告】
    单基因判「不一致」= 强信号（可信）
    单基因判「一致」= 弱信号（近缘复合体内 ITS 常分不开）
    所以本工具只用来【挑候选】，**不能作为「并名正确」的证据**。

用法
    python tools/nearest_neighbor_check.py --aln <已比对的 fasta> [--out out/nn.tsv]
    python tools/nearest_neighbor_check.py --aln X.fa --locus ITS   # 单基因也支持
输入约定
    已比对（等长）的 FASTA；头形如 `Csiam_18578`（组名_号）。
    组名 = 第一个下划线之前的部分（种缩写或种名）。`D.acericola_MFLUCC_17-0956` 也吃得下。
"""
import sys
try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

import argparse, io, os, re
from collections import defaultdict

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GAP = set("-~.?")


def read_aln(path):
    """读已比对 FASTA。返回 [(header, seq_upper)]；不是等长就报错（不猜、不自动比对）。"""
    recs, name = [], None
    with io.open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n").rstrip("\r")
            if line.startswith(">"):
                if name is not None:
                    recs.append((name, "".join(buf)))   # ⚠️ 必须 join —— buf 是行列表
                name, buf = line[1:].strip(), []
            elif line:
                buf.append(line.strip())
        if name is not None:
            recs.append((name, "".join(buf)))
    recs = [(n, s.upper()) for n, s in recs if s]
    lens = {len(s) for _, s in recs}
    if len(lens) > 1:
        raise SystemExit("序列不等长（{} 种长度）—— 本工具要求【已比对】的输入；"
                         "请先用 mafft 比对。不做隐式比对，免得给出不可比的数。".format(len(lens)))
    return recs


def group_of(header):
    h = re.sub(r"^\(strain\)\s*", "", header.strip())
    return h.split("_")[0] if "_" in h else h


def pid(s):
    return re.sub(r"[^A-Za-z0-9\.\-]", "", s)


def pdist(a, b):
    """p-距离：只比双方都不是 gap 的位置。可比位点为 0 时返回 None。"""
    same = diff = 0
    for x, y in zip(a, b):
        if x in GAP or y in GAP:
            continue
        if x == y:
            same += 1
        else:
            diff += 1
    n = same + diff
    return (diff / float(n)) if n else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--aln", required=True, help="已比对的 FASTA")
    ap.add_argument("--out", default=None, help="结果 TSV（默认 out/nearest_neighbor.tsv）")
    ap.add_argument("--locus", default="", help="位点名（只写进报告，单基因时填）")
    a = ap.parse_args()

    recs = read_aln(a.aln)
    labels = [pid(n) for n, _ in recs]
    groups = [group_of(n) for n, _ in recs]
    seqs = [s for _, s in recs]
    n = len(recs)
    gl = sorted(set(groups))
    print("=" * 96)
    print("最近邻检验（按株，不建树）")
    print("=" * 96)
    print("  输入：{}".format(os.path.relpath(a.aln, REPO)))
    print("  株数 {:,}   比对长度 {:,}   组（种）数 {:,}".format(n, len(seqs[0]), len(gl)))
    if len(gl) < 2:
        raise SystemExit("只有一个组 —— 无法做「组内 vs 组间」比较。")
    print("  位点：{}".format(a.locus or "(拼接/未标注)"))
    print("  计算 {} 对距离……".format(n * (n - 1) // 2), flush=True)

    d = [[None] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            v = pdist(seqs[i], seqs[j])
            d[i][j] = d[j][i] = v

    rows, flags = [], []
    for i in range(n):
        gin = [d[i][j] for j in range(n) if j != i and groups[j] == groups[i] and d[i][j] is not None]
        gout = [(d[i][j], j) for j in range(n) if groups[j] != groups[i] and d[i][j] is not None]
        if not gout:
            continue
        d_in_max = max(gin) if gin else None
        d_in_mean = (sum(gin) / len(gin)) if gin else None
        d_out_min, jo = min(gout)
        rows.append({"strain": labels[i], "group": groups[i],
                     "n_in_group": len(gin) + 1,
                     "d_in_max": d_in_max, "d_in_mean": d_in_mean,
                     "d_out_min": d_out_min, "nearest": labels[jo], "nearest_group": groups[jo]})
        if d_in_max is not None and d_out_min < d_in_max:
            flags.append(rows[-1])

    # 条形码间隙：每组「组内最大」vs「到最近他组的最小」
    gap = defaultdict(lambda: [0.0, None])
    for r in rows:
        g = r["group"]
        if r["d_out_min"] is not None:
            if gap[g][1] is None or r["d_out_min"] < gap[g][1]:
                gap[g][1] = r["d_out_min"]
        if r["d_in_max"] is not None:
            gap[g][0] = max(gap[g][0], r["d_in_max"])

    print()
    print("=" * 96)
    print("结果")
    print("=" * 96)
    print("  并名/错配候选（d_out < d_in）：{:,} / {:,} 株".format(len(flags), len(rows)))
    for r in flags[:15]:
        print("     {:<18} 声明组 {:<10} d_in_max={:.5f}  d_out_min={:.5f}  最近={}（{}）".format(
            r["strain"][:18], r["group"], r["d_in_max"], r["d_out_min"],
            r["nearest"][:14], r["nearest_group"]))
    print()
    print("  条形码间隙（组内最大 vs 到他组最近）：")
    no_gap = [g for g, (imax, omin) in gap.items() if omin is not None and imax >= omin]
    for g in sorted(gap)[:14]:
        imax, omin = gap[g]
        mark = "  ← 无间隙（组内比组间还远）" if (omin is not None and imax >= omin) else ""
        print("     {:<14} 组内max={:.5f}  他组min={}{}".format(
            g, imax, "{:.5f}".format(omin) if omin is not None else "—", mark))
    if len(gap) > 14:
        print("     …（共 {} 组）".format(len(gap)))
    print()
    print("  无间隙的组 {:,} / {:,} —— 这些种的序列与近缘种重叠，单基因无法分辨。".format(
        len(no_gap), len(gap)))

    outp = a.out or os.path.join(REPO, "out", "nearest_neighbor.tsv")
    if not os.path.isdir(os.path.dirname(outp)):
        os.makedirs(os.path.dirname(outp))
    with io.open(outp, "w", encoding="utf-8", newline="") as f:
        f.write("strain\tgroup\tn_in_group\td_in_max\td_in_mean\td_out_min\tnearest\tnearest_group\tflag\n")
        for r in rows:
            f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                r["strain"], r["group"], r["n_in_group"],
                "" if r["d_in_max"] is None else "{:.5f}".format(r["d_in_max"]),
                "" if r["d_in_mean"] is None else "{:.5f}".format(r["d_in_mean"]),
                "{:.5f}".format(r["d_out_min"]), r["nearest"], r["nearest_group"],
                "并名/错配候选" if r in flags else ""))
    print()
    print("  已写 {}".format(outp))
    print()
    print("  ⚠️ 结论强度不对称：单基因判「不一致」可信；判「一致」不可信（近缘复合体内分不开）。")
    print("     本工具只挑候选，不作为「并名正确」的证据。")
    return 0


if __name__ == "__main__":
    sys.exit(main())
