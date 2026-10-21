# -*- coding: utf-8 -*-
"""读 NCBI 结果表：给人看的视图。

表里「值」列是 JSON（因为不同方法的返回形状不同，JSON 是保证精确重放的最简载体），
本工具把它翻成人能读的形式。

用法：
    python tools\ncbi_table_info.py                      # 找 data\ncbi_*.tsv 全部列出
    python tools\ncbi_table_info.py data\ncbi_Colletotrichum.tsv
    python tools\ncbi_table_info.py --genus Colletotrichum --show 30
"""
import argparse
import collections
import glob
import io
import json
import os
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "src"))


def _decode(v):
    """把表里的内部标记解码：{"__tuple__": [...]} → 元组。
    （不解码的话，界面上会看到 __tuple__ 这种实现细节。）"""
    if isinstance(v, dict):
        if "__tuple__" in v and len(v) == 1:
            return tuple(_decode(x) for x in v["__tuple__"])
        return {k: _decode(x) for k, x in v.items()}
    if isinstance(v, list):
        return [_decode(x) for x in v]
    return v


def human(value, width=64):
    """把方法的返回值翻成一行人话。"""
    value = _decode(value)
    if isinstance(value, tuple):
        return "({})".format("、".join(str(_decode(x))[:40] for x in value))
    if isinstance(value, list):
        if not value:
            return "（空列表）"
        if all(isinstance(x, str) for x in value):
            s = "、".join(value[:6])
            return "{}{}  共{}项".format(s[:width], "…" if len(s) > width else "", len(value))
        return json.dumps(value[:3], ensure_ascii=False)[:width] + " 共{}项".format(len(value))
    if isinstance(value, dict):
        return json.dumps(value, ensure_ascii=False)[:width]
    return str(value)[:width]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("path", nargs="?", default=None)
    ap.add_argument("--genus", default=None)
    ap.add_argument("--show", type=int, default=15, help="每个方法展示几条")
    a = ap.parse_args()

    if a.path:
        paths = [a.path if os.path.isabs(a.path) else os.path.join(REPO, a.path)]
    else:
        pat = "ncbi_{}.tsv".format(a.genus) if a.genus else "ncbi_*.tsv"
        paths = sorted(glob.glob(os.path.join(REPO, "data", pat)))
    if not paths:
        print("找不到表。先跑一次流水线（表会随跑自动生成到 data/ncbi_<属名>.tsv）。")
        return 1

    for p in paths:
        if not os.path.isfile(p):
            print("❌ 不存在：{}".format(p))
            continue
        rows = [l.rstrip("\n").split("\t")
                for l in io.open(p, encoding="utf-8") if l.strip()]
        print("=" * 90)
        print("表：{}".format(os.path.relpath(p, REPO)))
        print("    {:,} 字节 / {} 条记录".format(os.path.getsize(p), max(0, len(rows) - 1)))
        mp = os.path.splitext(p)[0] + ".meta.json"
        if os.path.isfile(mp):
            m = json.loads(io.open(mp, encoding="utf-8").read())
            print("    生成/更新：{}".format(m.get("updated")))
            print("    代码指纹：{}".format(m.get("code_fingerprint")))
            print("    累计命中 {hits} / 联网新增 {misses} / 不可序列化跳过 {skipped_unserializable}".format(**m))
        print("=" * 90)

        by = collections.OrderedDict()
        for r in rows[1:]:
            if len(r) < 3:
                continue
            by.setdefault(r[0], []).append((r[1], r[2]))

        for method, items in sorted(by.items(), key=lambda x: -len(x[1])):
            print()
            print("【{}】{} 条".format(method, len(items)))
            for argjson, valjson in items[:a.show]:
                try:
                    args = json.loads(argjson)
                except Exception:
                    args = {"?": argjson[:40]}
                # 参数里最像「对象」的那一项，让人一眼看出是哪个种/哪个号
                key = (args.get("organism") or args.get("species") or args.get("code")
                       or args.get("number") or "")
                try:
                    val = json.loads(valjson)
                except Exception:
                    val = valjson
                extra = ", ".join("{}={}".format(k, v) for k, v in args.items()
                                  if k not in ("organism", "species", "code", "number"))
                print("    {:<38} {}  {}".format(
                    str(key)[:38], human(val), ("（{}）".format(extra) if extra else "")))
            if len(items) > a.show:
                print("    … 还有 {} 条（用 --show N 展开）".format(len(items) - a.show))
        print()
        print("=" * 90)
        print("说明：本表按【源类名.方法名 + 完整实参】为键、方法返回值(JSON)为值。")
        print("      重放时命中即零请求；未命中会联网并把结果补写进来。")
        print("      要强制重抓：删掉本表 + 同名 .meta.json 再跑。")
    return 0


if __name__ == "__main__":
    sys.exit(main())
