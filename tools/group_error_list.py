#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""group_error_list.py —— 把错误清单按**形态**分组输出。

为什么按形态分组（而不是按编号）：**形态决定发现手段**。
  · 查无型   值指向**不存在**的东西        → 机械检查即可（J-02 代码表、J-08 位数）
  · 安静型   值指向**另一个真实存在**的对象 → 必须靠**跨源对账**（J-01/J-13/J-14）；
             存在性检查全过——这正是它"安静"的原因
  · 归属型   值本身可能对，问题在**归属关系**（漏登／异名混入／一号两主）
  · 形式错配 形式错但**后果为零**（用于说明"发现的错值数 ≠ 危害"）

用法：python tools/group_error_list.py --src <副表1_错误清单.tsv> --out-dir <目录>
"""
import argparse
import csv
import io
import os
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8")
    sys.stderr.reconfigure(encoding="utf-8")
except Exception:
    pass

ORDER = ["查无型", "安静型", "归属型", "形式错配"]
GROUP = {
    "E1":  ("安静型",   "数字转置 → 落到另一个**真实存在**的号上（属对、号错）"),
    "E2":  ("安静型",   "相邻编号的**非模式菌**被当模式菌（指向真实存在的非模式菌株）"),
    "E6":  ("安静型",   "数字转置 → 落到**另一个属**的真实模式号上（最危险）"),
    "E9":  ("安静型",   "相邻编号的非模式菌被当模式菌；由**同源内部矛盾**暴露（图 vs 正文）"),
    "E11": ("形式错配", "号挂错行（形式错配成立，但**后果为零**：仅差 4 bp / 4136 bp）"),
    "E3":  ("查无型",   "保藏代码错字 → 该号**不存在**"),
    "E8":  ("查无型",   "保藏号丢字 → 该号**不存在**"),
    "E10": ("查无型",   "保藏号丢字（前导字符）→ 该号**不存在**"),
    "E4":  ("归属型",   "接受名的 ex-holotype **漏登**（值本身对，缺在另一侧）"),
    "E5":  ("归属型",   "异名模式号**混入**接受名（值属于另一个名）"),
    "E7":  ("归属型",   "**一号两主**（同号指派给两个种；方向已定为 tongrenense 侧错）"),
}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", required=True)
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args(argv)
    with io.open(a.src, encoding="utf-8-sig", newline="") as f:
        rd = csv.DictReader(f, delimiter="\t")
        fn, rows = list(rd.fieldnames), list(rd)
    for c in ["形态分类", "形态说明"]:
        if c not in fn:
            fn = fn[:1] + ["形态分类", "形态说明"] + fn[1:]
    for r in rows:
        g, d = GROUP.get(r["编号"], ("未分类", ""))
        r["形态分类"], r["形态说明"] = g, d
    outs = [(os.path.join(a.out_dir, os.path.basename(a.src)), rows),
            (os.path.join(a.out_dir, "副表1b_错误清单_按形态分组.tsv"),
             sorted(rows, key=lambda r: (ORDER.index(r["形态分类"]) if r["形态分类"] in ORDER else 99, r["编号"])))]
    for path, data in outs:
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
        with io.open(path, "w", encoding="utf-8-sig", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fn, delimiter="\t", lineterminator="\n")
            w.writeheader(); w.writerows(data)
        print("   写出 {}（{} 行）".format(path, len(data)))
    from collections import Counter
    print()
    print("   分组小计：{}".format({g: sum(1 for r in rows if r["形态分类"] == g) for g in ORDER}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
