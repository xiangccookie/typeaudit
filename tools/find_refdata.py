#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""自动找「参考数据」的位置，写进 registry——省得每次手填那三个长路径。

管理哪三个文件
--------------
  cbs            cbs_catalogue.tsv            CBS 保藏目录（判「号不存在 / 目录记为别种」）
  mycobank       mycobank_full.tsv            MycoBank 全库（判异名关系）
  biocollections NCBI_biocollections_full.csv 机构代码登记表（判「未知机构代码」）

用法
----
  python tools/find_refdata.py              # 探测并写 registry，打印结果
  python tools/find_refdata.py --show       # 只看当前 registry
  python tools/find_refdata.py --add cbs <path>   # 手动登记一条

registry 位置：~/.dsh/typeaudit_refs.json

为什么只做「探测」不做「下载」
------------------------------
三个文件的获取成本差三个数量级（本项目实测）：
  · 机构登记表  835 KB，有公开下载                → 可自动下（见 fetch_refdata.py）
  · MycoBank    61 MB / 40 万条，只能爬           → 小时级
  · CBS 目录    6.3 万个号、**官方无批量接口**
                只能逐号请求，1.3 s/次 ⇒ 约 22 小时
所以本脚本只负责「在本地已有的副本里挑最合适的那个」；
真要重新抓，用 tools/fetch_refdata.py（会如实告知耗时）。
"""
import argparse
import io
import json
import os
import sys
import time

try:
    sys.stdout.reconfigure(encoding="utf-8")
    sys.stderr.reconfigure(encoding="utf-8")
except Exception:
    pass

FILES = {
    "cbs": ("cbs_catalogue.tsv", "CBS 保藏目录"),
    "mycobank": ("mycobank_full.tsv", "MycoBank 全库"),
    "biocollections": ("NCBI_biocollections_full.csv", "机构登记表"),
}
# 搜索根目录 + 最大深度（太深会很慢）
ROOTS = [
    (os.path.join(os.path.expanduser("~"), "Desktop"), 4),
    (os.path.join(os.path.expanduser("~"), ".dsh"), 5),
    (os.path.join(os.path.expanduser("~"), "Documents"), 3),
]


def registry_path():
    d = os.path.join(os.path.expanduser("~"), ".dsh")
    os.makedirs(d, exist_ok=True)
    return os.path.join(d, "typeaudit_refs.json")


def load():
    p = registry_path()
    if os.path.exists(p):
        try:
            return json.load(io.open(p, encoding="utf-8"))
        except Exception:
            pass
    return {}


def save(d):
    io.open(registry_path(), "w", encoding="utf-8").write(
        json.dumps(d, ensure_ascii=False, indent=1))


def probe(verbose=True):
    """在常见根目录下找三个文件；同名多份时取**最新修改**的那份。

    为什么取最新：本项目里 CBS 目录曾有多份副本（含 _SUPERSEDED 下的旧版），
    取到旧版会让核验结果整体失真。
    """
    out = {}
    for key, (fn, label) in FILES.items():
        hits = []
        for root, depth in ROOTS:
            if not os.path.isdir(root):
                continue
            base_depth = root.rstrip(os.sep).count(os.sep)
            for dirpath, dirnames, filenames in os.walk(root):
                if dirpath.count(os.sep) - base_depth >= depth:
                    dirnames[:] = []
                if fn in filenames:
                    p = os.path.join(dirpath, fn)
                    try:
                        hits.append((os.path.getmtime(p), os.path.getsize(p), p))
                    except OSError:
                        pass
        if hits:
            hits.sort(reverse=True)          # 最新优先
            out[key] = {"path": hits[0][2], "mtime": hits[0][0], "size": hits[0][1],
                        "n_copies": len(hits),
                        "others": [h[2] for h in hits[1:4]]}
            if verbose:
                print("  [{}] {}  ({:,} B, {} 份副本)".format(
                    key, hits[0][2], hits[0][1], len(hits)))
                for _, sz, p in hits[1:4]:
                    print("        备选（较旧）：{}  ({:,} B)".format(p, sz))
        else:
            if verbose:
                print("  [{}] ❌ 未找到 {}（{}）".format(key, fn, label))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--show", action="store_true")
    ap.add_argument("--add", nargs=2, metavar=("KEY", "PATH"))
    a = ap.parse_args()

    if a.show:
        d = load()
        print("registry：{}".format(registry_path()))
        for k, v in d.items():
            print("  {:<16}{}".format(k, v.get("path")))
            print("      {:,} B  改动于 {}".format(
                v.get("size", 0),
                time.strftime("%Y-%m-%d %H:%M", time.localtime(v.get("mtime", 0)))))
        if not d:
            print("  （空；先跑一次不带参数的 find_refdata.py）")
        return 0

    if a.add:
        key, path = a.add
        if key not in FILES:
            print("❌ key 必须是 {}".format("/".join(FILES)))
            return 2
        if not os.path.exists(path):
            print("❌ 路径不存在：{}".format(path))
            return 2
        d = load()
        d[key] = {"path": os.path.abspath(path), "mtime": os.path.getmtime(path),
                  "size": os.path.getsize(path), "n_copies": 1, "others": [],
                  "manual": True}
        save(d)
        print("✅ 已登记 {} → {}".format(key, path))
        return 0

    print("探测参考数据（根目录：{}）".format("、".join(r for r, _ in ROOTS)))
    found = probe()
    d = load()
    for k, v in found.items():
        d[k] = v
    save(d)
    print("\n已写入 registry：{}".format(registry_path()))
    miss = [k for k in FILES if k not in found]
    if miss:
        print("未找到：{}（这三个文件仓库不带，须自备；见 README）".format(
            "、".join(FILES[m][1] for m in miss)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
