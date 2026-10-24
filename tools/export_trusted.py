# -*- coding: utf-8 -*-
"""export_trusted.py —— 把网页里「受信菌株表」的实时计算结果落盘成文件。

为什么要落盘
------------
网页的 /api/trusted 是**每次请求现算**的，刷新就没了。
论文若要用它作补充表，必须先有稳定文件（可引用、可核对、可提交）。

⚠️ 「受信」的含义**仅限于「号与外部目录/登记表一致」**——
   不等于这株菌真是模式菌株，也不等于序列属于该种。
   只覆盖 CBS / ICMP / CGMCC 三家；其余机构一律不进表（把「没检查」当「受信」是错的）。

用法
----
    python tools/export_trusted.py                 # 三个属全导
    python tools/export_trusted.py --genus colletotrichum
"""
import argparse
import importlib.util
import io
import json
import os
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GEN = [("colletotrichum", "out_colleto_web"), ("trichoderma", "out_tricho_web"),
       ("diaporthe", "out_dia_web")]
COLS = ["种名", "模式菌号", "菌株号", "模式号说明", "菌株号说明"]


def load_web():
    spec = importlib.util.spec_from_file_location("pw", os.path.join(REPO, "tools", "pipeline_web.py"))
    mod = importlib.util.module_from_spec(spec)
    sys.modules["pw"] = mod
    spec.loader.exec_module(mod)
    return mod


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genus", default="all")
    ap.add_argument("--out", default="out_supp")
    a = ap.parse_args()
    pw = load_web()
    targets = GEN if a.genus == "all" else [g for g in GEN if g[0] == a.genus]
    od = os.path.join(REPO, a.out)
    os.makedirs(od, exist_ok=True)
    summ = []
    for slug, outdir in targets:
        _, rows = pw.read_trusted(os.path.join(REPO, outdir), slug)
        nums = set()
        for r in rows:
            for k in ("strain", "type"):
                for x in (r.get(k) or "").split(";"):
                    if x.strip():
                        nums.add(x.strip())
        p = os.path.join(od, "{}_受信菌株表.tsv".format(slug))
        with io.open(p, "w", encoding="utf-8-sig", newline="") as f:
            f.write("\t".join(COLS) + "\n")
            for r in rows:
                f.write("\t".join([r.get("species", ""), r.get("type", ""), r.get("strain", ""),
                                   r.get("why_type", ""), r.get("why_strain", "")]) + "\n")
        summ.append({"属": slug, "有受信号的种": len(rows), "受信号（去重）": len(nums),
                     "受信号（人次）": sum(1 for r in rows for k in ("strain", "type")
                                          for x in (r.get(k) or "").split(";") if x.strip()),
                     "文件": p})
        print("{:<16} 种 {:>4} ｜ 去重号 {:>4} ｜ {}".format(slug, len(rows), len(nums),
                                                           os.path.relpath(p, REPO)))
    json.dump(summ, io.open(os.path.join(od, "受信菌株表_汇总.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=2)
    print("\n汇总 → {}".format(os.path.join(od, "受信菌株表_汇总.json")))


if __name__ == "__main__":
    main()
