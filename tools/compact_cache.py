# -*- coding: utf-8 -*-
# --- console utf8 bootstrap（Windows 控制台默认 GBK，输出 ↔ ⚠️ 🔴 会崩；详见 src/typeaudit/_console.py）
import sys as _sys
for _s in (_sys.stdout, _sys.stderr):
    try:
        _s.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass
r"""缓存瘦身：一次性的磁盘回收。

做两件事（默认预演，--apply 才落笔）：

  一、**就地剥掉旧缓存里的碱基序列**（ORIGIN 段）。
      旧缓存是加 _strip_origin 之前下载的，实测 abscissum 两份分别 90.9 MB / 88.1 MB，
      其中 74% 是审计根本用不到的碱基。剥掉后解析结果不变（见 sources._strip_origin）。

  二、**收掉重复的 gbx_ 副本**。
      `gbx_<哈希>.gb` 是按请求内容寻址的新名字，旧缓存是 cc_/sv_/tm_ 这种调用方名字。
      两者内容相同时，早先的实现会复制一份 ⇒ 凭空多出约 179 MB。
      现在改用 gbx_index.json 记「新键 → 已有文件名」，副本可以删掉。
      ⚠️ 删之前逐字节比对（sha256），并且把别名写进索引，确保还找得到。

    python tools\compact_cache.py            # 预演
    python tools\compact_cache.py --apply    # 落笔
"""
import argparse
import collections
import hashlib
import io
import json
import os
import re
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
CACHE = os.path.join(REPO, ".cache")
INDEX = os.path.join(CACHE, "gbx_index.json")
ORIGIN_RE = re.compile(r"^ORIGIN\b.*?(?=^//)", re.S | re.M)


def sha(p):
    h = hashlib.sha256()
    with io.open(p, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--apply", action="store_true")
    a = ap.parse_args()

    files = [f for f in os.listdir(CACHE)
             if os.path.isfile(os.path.join(CACHE, f)) and f.endswith(".gb")]
    print("缓存目录：{}".format(CACHE))
    print("  .gb 文件 {} 个，合计 {:,} 字节".format(
        len(files), sum(os.path.getsize(os.path.join(CACHE, f)) for f in files)))
    print()

    # ---- 一、就地剥序列 ----
    print("=" * 74)
    print("一、就地剥掉 ORIGIN（碱基序列）")
    print("=" * 74)
    stripped = 0
    saved1 = 0
    for f in sorted(files):
        p = os.path.join(CACHE, f)
        try:
            s = io.open(p, encoding="utf-8", errors="replace").read()
        except Exception:
            continue
        if "\nORIGIN" not in s:
            continue
        s2 = ORIGIN_RE.sub("", s)
        b1, b2 = len(s.encode()), len(s2.encode())
        if b2 >= b1:
            continue
        stripped += 1
        saved1 += b1 - b2
        print("   {:>12,} → {:>12,} B  （省 {:>5.1f}%）  {}".format(
            b1, b2, 100.0 * (b1 - b2) / b1, f[:52]))
        if a.apply:
            try:
                io.open(p, "w", encoding="utf-8", errors="replace").write(s2)
            except Exception as e:
                print("      ⚠️ 写回失败：{}".format(e))
    print("   合计：{} 个文件，省 {:,} 字节 = {:.1f} MB".format(
        stripped, saved1, saved1 / 1e6))

    # ---- 二、收掉重复的 gbx_ ----
    print()
    print("=" * 74)
    print("二、收掉重复的 gbx_ 副本（改为写索引）")
    print("=" * 74)
    gbx = [f for f in files if f.startswith("gbx_")]
    others = [f for f in files if not f.startswith("gbx_")]
    # 按大小分组，只在同尺寸内算 sha，省时间
    bysize = collections.defaultdict(list)
    for f in others:
        bysize[os.path.getsize(os.path.join(CACHE, f))].append(f)
    idx = {}
    if os.path.isfile(INDEX):
        try:
            idx = json.load(io.open(INDEX, encoding="utf-8")) or {}
        except Exception:
            idx = {}
    removed = saved2 = 0
    for f in sorted(gbx):
        p = os.path.join(CACHE, f)
        sz = os.path.getsize(p)
        cand = bysize.get(sz) or []
        if not cand:
            print("   [留] {:<44} {:>12,} B（无同尺寸对象，保留）".format(f[:44], sz))
            continue
        h = sha(p)
        twin = None
        for g in cand:
            if sha(os.path.join(CACHE, g)) == h:
                twin = g
                break
        if not twin:
            print("   [留] {:<44} {:>12,} B（同尺寸但内容不同）".format(f[:44], sz))
            continue
        removed += 1
        saved2 += sz
        print("   [删] {:<44} {:>12,} B  与 {} 完全相同".format(f[:44], sz, twin[:30]))
        idx[f] = twin
        if a.apply:
            try:
                os.remove(p)
            except Exception as e:
                print("      ⚠️ 删除失败：{}".format(e))
                idx.pop(f, None)
                removed -= 1
                saved2 -= sz
    print("   合计：删 {} 个副本，省 {:,} 字节 = {:.1f} MB".format(
        removed, saved2, saved2 / 1e6))
    if a.apply and idx:
        io.open(INDEX, "w", encoding="utf-8").write(
            json.dumps(idx, ensure_ascii=False, indent=0))
        print("   已写索引：{}（{} 条别名）".format(
            os.path.relpath(INDEX, REPO), len(idx)))

    print()
    print("=" * 74)
    print("合计可回收 {:,} 字节 = {:.1f} MB".format(saved1 + saved2, (saved1 + saved2) / 1e6))
    if not a.apply:
        print()
        print("这是预演（没有落笔）。确认后：")
        print("    python tools\\compact_cache.py --apply")
    else:
        tot = 0
        for root, dirs, fs in os.walk(REPO):
            dirs[:] = [d for d in dirs if d not in (".git", "__pycache__")]
            for f in fs:
                try:
                    tot += os.path.getsize(os.path.join(root, f))
                except Exception:
                    pass
        print("   仓库现在 {:,} 字节 = {:.1f} MiB".format(tot, tot / 1024 / 1024))
    print("=" * 74)
    return 0


if __name__ == "__main__":
    sys.exit(main())
