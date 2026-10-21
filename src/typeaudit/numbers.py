# -*- coding: utf-8 -*-
"""核数器：文档里「声明的数字」↔ 源表现场重算，报不一致。

为什么需要（实测代价）：人工做这件事时一次逮到 4 处旧值残留 + 22 行快照漂移，
全是「改了权威值、下游没跟」的级联漏改。漏检代价高——读者会把旧值当结论。

**规则全部外部化**（JSON），故本工具对任何项目都通用：
    python -m typeaudit.numbers --rules rules/example_rules.json

两条实现约束（踩过坑）：
  约束一·分区：改文档必须区分**正文**与**历史更正记录**。同一串可能在两处出现，
      历史处必须保留。故支持 body_marker：扫描时裁掉标记之后的段落。
  约束二·先确认源表：报「算不出来」之前，先确认自己读的是不是原句声明的那张表。
      实测曾因在 A 表里复算某比例而误判"算不出"，实际它的出处是 B 表。
"""
import argparse
import csv
import io
import json
import os
import re
import sys

OPS = {
    "count_rows":   lambda T, a: str(len(T[a["tsv"]])),
    "sum_col":      lambda T, a: str(sum(int(r[a["col"]]) for r in T[a["tsv"]]
                                          if str(r.get(a["col"], "")).strip().isdigit())),
    "value_counts": lambda T, a: str(sum(1 for r in T[a["tsv"]]
                                          if (r.get(a["col"]) or "").strip() == a["key"])),
    "lookup":       lambda T, a: str(next((r.get(a["col"], "") for r in T[a["tsv"]]
                                            if (r.get(a["key_col"]) or "").strip() == a["key"]), "")),
    "const":        lambda T, a: str(a["value"]),
}


def load_tsv(path):
    if not os.path.exists(path):
        return []
    with io.open(path, encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def read_doc(path, body_marker=None):
    if not os.path.exists(path):
        raise FileNotFoundError("文档找不到：{}".format(path))
    s = io.open(path, encoding="utf-8-sig", errors="replace").read()
    if body_marker:
        i = s.find(body_marker)
        if i > 0:
            s = s[:i]
    return s


def run(rules_path, docs_dir=None, data_dir=None, report=None):
    cfg = json.load(io.open(rules_path, encoding="utf-8"))
    docs_dir = docs_dir or cfg.get("docs_dir", ".")
    data_dir = data_dir or cfg.get("data_dir", ".")
    marker = cfg.get("body_marker")
    T = {k: load_tsv(os.path.join(data_dir, v)) for k, v in cfg.get("tables", {}).items()}

    rows, ok = [], 0
    print("=" * 96)
    print("【核数器：文档声明数字 ↔ 源表重算】")
    print("=" * 96)
    for r in cfg["rules"]:
        try:
            doc = read_doc(os.path.join(docs_dir, r["doc"]), marker)
        except FileNotFoundError as e:
            print("   [MISS] {:<30}{}".format(r["name"], e))
            rows.append([r["name"], r["doc"], "MISS", "", ""])
            continue
        m = re.search(r["regex"], doc)
        try:
            want = OPS[r["op"]](T, r.get("args", {}))
        except Exception as e:
            print("   [ERR ] {:<30}源表计算失败：{}".format(r["name"], str(e)[:50]))
            rows.append([r["name"], r["doc"], "ERR", "", str(e)[:60]])
            continue
        if not m:
            status, got = "MISS", "(该串不存在)"
        else:
            got_g = [m.group(i) for i in range(1, m.lastindex + 1)] if m.lastindex else [m.group(0)]
            want_g = [want] if len(got_g) == 1 else want.split(" / ")
            if len(got_g) != len(want_g):
                status = "SHAPE"
            else:
                status = "PASS" if got_g == want_g else "FAIL"
            got = " / ".join(got_g)
            want = " / ".join(want_g)
        ok += status == "PASS"
        rows.append([r["name"], r["doc"], status, got, want])
        print("   [{}] {:<30}{:<20}{}".format(
            status, r["name"][:28], r["doc"][:18], "" if status == "PASS" else "找到={} 应为={}".format(got, want)))
    print()
    print("   PASS {}/{}".format(ok, len(cfg["rules"])))
    if report:
        d = os.path.dirname(os.path.abspath(report))
        os.makedirs(d, exist_ok=True)
        with io.open(report, "w", encoding="utf-8-sig", newline="") as f:
            w = csv.writer(f, delimiter="\t", lineterminator="\n")
            w.writerow(["检查项", "文档", "结果", "文档中找到", "源表重算"])
            w.writerows(rows)
        print("   报告落盘：{}".format(report))
    return ok, len(cfg["rules"])


def main(argv=None):
    ap = argparse.ArgumentParser(description="核数器：文档声明数字 ↔ 源表重算")
    ap.add_argument("--rules", required=True, help="规则 JSON 路径")
    ap.add_argument("--docs-dir", help="文档目录（覆盖规则文件里的设置）")
    ap.add_argument("--data-dir", help="数据目录（覆盖规则文件里的设置）")
    ap.add_argument("--report", help="报告 TSV 输出路径")
    a = ap.parse_args(argv)
    ok, n = run(a.rules, a.docs_dir, a.data_dir, a.report)
    return 0 if ok == n else 1


if __name__ == "__main__":
    sys.exit(main())
