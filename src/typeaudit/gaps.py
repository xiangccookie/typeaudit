# -*- coding: utf-8 -*-
"""缺环分类器：把一条条「判不了」的挂账项，分类成**可交接的任务单**。

设计立场：**只出清单，不裁定。** 这些缺环的补齐动作需要外部资源或人工判断，
脚本越界去猜就会制造假证据。列表本身就是可交付物——它把「还有 N 项没结」
变成「这 N 项各自缺什么、找谁要、影响多大」。

配置外部化：见 config/gaps_example.json。
"""
import argparse
import csv
import io
import json
import os
import re
import sys
from collections import Counter


def load_tsv(path):
    with io.open(path, encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def classify(text, classes, default="no_local_evidence"):
    """返回命中的缺环类（按配置顺序、去重）。配置顺序即优先级。"""
    out = []
    for c in classes:
        if re.search(c["regex"], text) and c["name"] not in out:
            out.append(c["name"])
    return out or [default]


def collect(cfg, data_dir):
    tasks = []
    for src in cfg["sources"]:
        for r in load_tsv(os.path.join(data_dir, src["tsv"])):
            blob = " ".join(r.get(c) or "" for c in src.get("text_cols", []))
            if cfg.get("blocker_regex") and not re.search(cfg["blocker_regex"], blob):
                continue
            tasks.append({
                "来源": src.get("label", src["tsv"]),
                "对象": " · ".join((r.get(c) or "").strip() for c in src["id_cols"]).strip(" ·"),
                "现判定": (r.get(src.get("status_col", "")) or "")[:60],
                "缺环类型": classify(blob, cfg["classes"]),
                "证据": (r.get(src.get("evidence_col", "")) or "")[:160] if src.get("evidence_col") else "",
            })
    return tasks


def main(argv=None):
    ap = argparse.ArgumentParser(description="缺环分类器：挂账项 -> 可交接任务单")
    ap.add_argument("--config", required=True, help="配置 JSON 路径")
    ap.add_argument("--data-dir", default=".", help="数据目录")
    ap.add_argument("--out", help="任务单 TSV 输出路径")
    a = ap.parse_args(argv)
    cfg = json.load(io.open(a.config, encoding="utf-8"))
    tasks = collect(cfg, a.data_dir)
    needs, impact, order = cfg["needs"], cfg["impact"], cfg.get("impact_order", [])

    print("=" * 96)
    print("【缺环任务单】共 {} 项挂账".format(len(tasks)))
    print("=" * 96)
    rows = []
    for i, t in enumerate(tasks, 1):
        cls = t["缺环类型"]
        tid = "G{:02d}".format(i)
        key = sorted(cls, key=lambda c: order.index(c) if c in order else 99)[0]
        rows.append([tid, t["来源"], t["对象"], t["现判定"], " ＋ ".join(cls),
                     " ＋ ".join(needs.get(c, c) for c in cls), impact.get(key, ""), t["证据"]])
        print()
        print("   {}  {}".format(tid, t["对象"]))
        print("       缺环类型  {}".format(" ＋ ".join(cls)))
        print("       需要什么  {}".format(" ＋ ".join(needs.get(c, c) for c in cls)))
        print("       影响      {}".format(impact.get(key, "")))
    print()
    print("【分类汇总】")
    for k, v in Counter(c for t in tasks for c in t["缺环类型"]).most_common():
        print("   {:<22}{:>3} 项   {}".format(k, v, needs.get(k, "")))
    if a.out:
        d = os.path.dirname(os.path.abspath(a.out))
        os.makedirs(d, exist_ok=True)
        with io.open(a.out, "w", encoding="utf-8-sig", newline="") as f:
            w = csv.writer(f, delimiter="\t", lineterminator="\n")
            w.writerow(["任务ID", "来源", "对象", "现判定", "缺环类型", "需要什么资源", "对结论的影响", "证据"])
            w.writerows(rows)
        print()
        print("   任务单落盘：{}".format(a.out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
