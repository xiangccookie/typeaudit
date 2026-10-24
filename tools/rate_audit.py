#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""rate_audit.py —— 从「证据链」重算可核号与可核号内错值（号级口径）。

为什么需要它
------------
论文的核心比率 12.6% 一直是硬编码在 tools/make_figs.py 第 35 行：

    SAME, DIFF, REAL, FALSEP = 2634, 558, 402, 156

其中 558（可核号内的号级不一致）只在一份说明文档里有记录
（交付_三属对照_20261022/06_反向审计/反向审计_说明.md），**没有可重跑的脚本**。
对数据论文来说这是硬伤：读者无法复现这个分母与分子。

本脚本用流水线自己落盘的证据链（out_<genus>/cNN/<genus>_证据.jsonl，
每个种一行、每个号一条记录，含 flags）重算，口径写死在这里，可逐项核对。

口径（与说明文档逐条对应）
--------------------------
1. 可核号 = 号级单元格中，机构属于 CBS / ICMP / CGMCC 的那些
   （这三家有公开机读目录，能核对；其余机构没有，不进分母）。
2. 计数按【号级单元格的出现次数】：同一种下同一号若在「菌株号」与「模式菌」
   两列各出现一次，算两个单元格（与三列表的展示一致）。
3. 错值 = 该单元格的任一标记被 pipeline_web.classify_flag 判为「错值」。
4. 修正后（反向审计）= 错值中再扣除「号不存在但 GenBank 有引用」的，
   这部分在 06_反向审计 里逐条核过（156 条），本脚本不改判，只报未扣除值。

用法
----
    python tools/rate_audit.py                      # 打印三属对照
    python tools/rate_audit.py --json out.json      # 同时落盘

已知差异（诚实声明）
--------------------
本脚本结果与说明文档的记录值不完全一致（可核号 3,036 vs 3,192；错值 572 vs 558）。
差异来源已定位但不逐条解释：
  · 原口径可能还并入了「序列轴」与「四方核验」里出现、但没进三列表的号；
  · 同一号在两列重复出现时，原口径的选择规则未记录。
因此本脚本的价值不是「证明 558 错了」，而是【给出一个可复现的口径】，
让读者能自己算、也能看出与旧值的差距有多大（约 2–5%）。
"""
import argparse
import glob
import io
import json
import os
import re
import sys
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, HERE)

GENERA = [("Colletotrichum", "colletotrichum", "out_v2"),
          ("Trichoderma", "trichoderma", "out_trichoderma"),
          ("Diaporthe", "diaporthe", "out_diaporthe")]
VERIFIABLE = {"CBS", "ICMP", "CGMCC"}

# 说明文档里记录的原值，用于对照
LEGACY = {"Colletotrichum": (1553, 302), "Trichoderma": (970, 174), "Diaporthe": (669, 82)}


def load_classify():
    """复用网页侧的标记归类函数，保证与论文口径同源。"""
    import importlib.util as ilu
    spec = ilu.spec_from_file_location("_pw", os.path.join(HERE, "pipeline_web.py"))
    mod = ilu.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.classify_flag


def code_of(numkey):
    """号键 'CBS149616' → 'CBS'"""
    m = re.match(r"([A-Za-z]+)", numkey or "")
    return m.group(1).upper() if m else ""


PREFIX = re.compile(r"^\((?:strain|isolate|culture|specimen|other)\)\s*", re.I)
TOK = re.compile(r"^([A-Za-z][\w\-\.]*)\s*[:#]?\s*([\d\.]+[A-Za-z]?)\s*(?:\[(.*?)\])?\s*$")


def audit_genus(genus, slug, classify):
    """主口径：三列表里 CBS/ICMP/CGMCC 号的出现次数（与展示一致，两列各计一次）。

    为什么用三列表而不用证据链：证据链只落盘「有标记」的号（实测 1,871），
    没有问题的号不入账；三列表把「[—]」也写出来，因此更接近文档的 3,192。
    """
    p = os.path.join(REPO, "results", genus, "{}_三列表.tsv".format(slug))
    occ = wrong = 0
    by_code = Counter()
    if not os.path.exists(p):
        return {"可核号": 0, "错值": 0, "按标记": {}}
    for row in io.open(p, encoding="utf-8-sig").read().splitlines()[1:]:
        c = row.split("\t")
        if len(c) < 3:
            continue
        for cell in c[1:3]:
            for seg in cell.split("｜"):
                s = PREFIX.sub("", seg.strip())
                if not s or s == "—":
                    continue
                m = TOK.match(s)
                if not m:
                    continue
                code = m.group(1).upper()
                if code not in VERIFIABLE:
                    continue
                occ += 1
                by_code[code] += 1
                fl = (m.group(3) or "").strip()
                if fl and fl != "—" and classify(fl) == "错值":
                    wrong += 1
    return {"可核号": occ, "错值": wrong, "按标记": dict(by_code)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", default=None, help="把结果写成 JSON")
    a = ap.parse_args()
    classify = load_classify()

    print("从证据链复算（口径见文件头说明）")
    print("{:<16}{:>10}{:>10}{:>12}".format("属", "可核号", "错值", "文档记录"))
    t_ver = t_err = 0
    detail = {}
    for g, gl, od in GENERA:
        r = audit_genus(g, gl, classify)
        detail[g] = r
        t_ver += r["可核号"]
        t_err += r["错值"]
        e_ver, e_err = LEGACY[g]
        print("{:<16}{:>10}{:>10}{:>12}".format(
            g, r["可核号"], r["错值"], "{} / {}".format(e_ver, e_err)))
        # 机构构成分项 —— 图 2B 需要，与「可核号」同一口径
        bc = r.get("按标记") or {}
        if bc:
            tot = sum(bc.values())
            print("    └ 构成：CBS {:>5} ｜ ICMP {:>4} ｜ CGMCC {:>4}   （ICMP {:.1f}%）".format(
                bc.get("CBS", 0), bc.get("ICMP", 0), bc.get("CGMCC", 0),
                100.0 * bc.get("ICMP", 0) / tot if tot else 0))
    print("{:<16}{:>10}{:>10}{:>12}".format("合计", t_ver, t_err, "3192 / 558"))
    print("\n与文档的差：可核号 {:+.0f}（{:.1f}%），错值 {:+.0f}（{:.1f}%）".format(
        t_ver - 3192, 100.0 * (t_ver - 3192) / 3192, t_err - 558, 100.0 * (t_err - 558) / 558))
    print("错值率（未扣反向审计）：{:.1f}%（文档 17.5%）".format(100.0 * t_err / t_ver if t_ver else 0))
    for lbl, mv in (("保守", 156), ("乐观", 176)):
        n = t_err - mv
        print("  扣反向审计（{}，移出 {} 条）：{}/{} = {:.1f}%".format(
            lbl, mv, n, t_ver, 100.0 * n / t_ver if t_ver else 0))
    if a.json:
        payload = {"各属": detail, "合计": {"可核号": t_ver, "错值": t_err},
                   "文档记录": {"可核号": 3192, "错值": 558}}
        with io.open(a.json, "w", encoding="utf-8", newline="\n") as f:
            f.write(json.dumps(payload, ensure_ascii=False, indent=2))
        print("已写出:", a.json)
    return 0


if __name__ == "__main__":
    sys.exit(main())
