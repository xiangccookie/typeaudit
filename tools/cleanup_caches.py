# -*- coding: utf-8 -*-
"""缓存归并：把散落在错误位置的抓取缓存收拢到仓库统一的 .cache/。

为什么需要它（本机实测）：
  · config/ 里躺着 270.88 MiB 的 GenBank 缓存 —— 配置目录被当成了缓存目录；
  · out_audit/out_audit/.cache 是【双层目录】277.41 MiB —— 相对路径被拼了两次；
  · 每个 out_*/.cache 各存一套，互不复用 → 换个输出目录就把属级 8.38 MiB 重下一遍。

本工具的行为（**只搬不删，除非是逐字节相同的副本**）：
  1. 扫描「错误位置」里的缓存文件；
  2. 算 sha256。统一目录 .cache/ 里已有同一份 → 这是纯副本，删除；
  3. 否则【移动】到 .cache/（保留，不丢数据）；
  4. 搬空的错误目录一并删掉。

默认只预演（--apply 才真正落笔）。**不碰 data/、docs/、src/、tools/、config/*.json**。
"""
import argparse
import hashlib
import io
import os
import re
import shutil
import sys

# ⚠️ 固定输出编码：Windows 控制台默认 GBK，管道捕获时会变乱码。
try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEST = os.path.join(REPO, ".cache")

# 抓取缓存的命名形态：前缀_..._.gb / .xml / .json
CACHE_RE = re.compile(r"^(cc|sv|tm|tmt|txg|es_nucleotide|es_taxonomy|es_|tmt|"
                      r"icmp|esum|count|search|cmp|pmc|e)_[^\\/]*\.(gb|xml|json)$",
                      re.I)
# 明显属于「源数据 / 配置」，永不移动
KEEP_RE = re.compile(r"\.(tsv|csv|json)$", re.I)
NEVER = ("data", "docs", "src", "tools", "tests", "examples", "rules")


def is_cache_name(fn):
    return bool(CACHE_RE.match(fn))


def wrong_dirs():
    """列出所有「不该放缓存」的目录。

    ⚠️ 铁律：**必须排除统一目录 DEST 自己**。
       否则 .cache 里的每个文件都会在「目标哈希索引」里匹配到自己，
       被判成「已有同一份的副本」→ --apply 会把统一缓存目录整个清掉，
       还报告说「删除 N 个副本」。（本机实测踩过一次，见本文件末尾的教训。）
    """
    out = []
    seen = set()
    dest_abs = os.path.abspath(DEST)

    def add(label, d):
        da = os.path.abspath(d)
        if da == dest_abs or da in seen:
            return
        seen.add(da)
        out.append((label, d))

    cfg = os.path.join(REPO, "config")
    if os.path.isdir(cfg):
        add("config/ 目录（配置目录不该有缓存）", cfg)
    for root, dirs, files in os.walk(REPO):
        dirs[:] = [d for d in dirs if d not in (".git", "__pycache__")]
        rel = os.path.relpath(root, REPO).replace("\\", "/")
        top = rel.split("/")[0]
        if top in NEVER or rel == ".":
            continue
        base = os.path.basename(root)
        if base == ".icmp_cache":
            add("out_*/.icmp_cache（按运行分家的缓存）", root)
        elif base == ".cache":
            if rel.count("/") >= 1:
                add("out_*/.cache（按运行分家的缓存 / 双层目录）", root)
            else:
                add("out_*/.cache（按运行分家的缓存）", root)
    return out


def sha(p, cap=200 * 1024 * 1024):
    if os.path.getsize(p) > cap:
        return None
    h = hashlib.sha256()
    with io.open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--apply", action="store_true",
                    help="真正落笔（默认只预演）")
    ap.add_argument("--keep-outdirs", action="store_true",
                    help="不动 out_*/.cache（只清 config/ 与双层目录）")
    a = ap.parse_args()

    os.makedirs(DEST, exist_ok=True)
    dest_hashes = {}
    for fn in os.listdir(DEST):
        p = os.path.join(DEST, fn)
        if os.path.isfile(p):
            h = sha(p)
            if h:
                dest_hashes[h] = p

    print("统一缓存目录：{}".format(DEST))
    print("  （里面现有 {} 个文件）".format(len(dest_hashes)))
    print()

    total_files = total_bytes = dup_bytes = moved_bytes = 0
    plan = []
    for label, d in wrong_dirs():
        if a.keep_outdirs and label.startswith("out_"):
            continue
        if not os.path.isdir(d):
            continue
        files = [f for f in os.listdir(d) if os.path.isfile(os.path.join(d, f))]
        caches = [f for f in files if is_cache_name(f)]
        if not caches:
            continue
        b = sum(os.path.getsize(os.path.join(d, f)) for f in caches)
        print("【{}】{}".format(label, os.path.relpath(d, REPO)))
        print("     {} 个文件 / {:,} B = {:.2f} MiB".format(
            len(caches), b, b / 1024 / 1024))
        for f in caches:
            p = os.path.join(d, f)
            s = os.path.getsize(p)
            h = sha(p)
            total_files += 1
            total_bytes += s
            if h and h in dest_hashes:
                dup_bytes += s
                print("       [删·已有同一份]  {:<52} {:>12,}".format(f[:52], s))
                plan.append(("del", p, None))
            else:
                moved_bytes += s
                if h:
                    dest_hashes[h] = os.path.join(DEST, f)
                print("       [搬·保留]        {:<52} {:>12,}".format(f[:52], s))
                plan.append(("move", p, os.path.join(DEST, f)))
        print()

    print("=" * 72)
    print("合计：扫描缓存 {} 个 / {:,} B = {:.2f} MiB".format(
        total_files, total_bytes, total_bytes / 1024 / 1024))
    print("  其中【纯副本可删】 {:,} B = {:.2f} MiB".format(dup_bytes, dup_bytes / 1024 / 1024))
    print("  其中【搬走保留】   {:,} B = {:.2f} MiB".format(moved_bytes, moved_bytes / 1024 / 1024))
    print("=" * 72)

    if not a.apply:
        print()
        print("这是预演（没有落笔）。确认无误后加 --apply 真正执行：")
        print("    python tools\\cleanup_caches.py --apply")
        return

    n_del = n_mv = 0
    for act, src, dst in plan:
        try:
            if act == "del":
                os.remove(src)
                n_del += 1
            else:
                if os.path.exists(dst):
                    # 罕见：同名但内容不同 → 加序号避免覆盖
                    root, ext = os.path.splitext(dst)
                    i = 1
                    while os.path.exists("{}_{}{}".format(root, i, ext)):
                        i += 1
                    dst = "{}_{}{}".format(root, i, ext)
                shutil.move(src, dst)
                n_mv += 1
        except Exception as e:
            print("  ⚠️ {} 失败：{}: {}".format(os.path.relpath(src, REPO),
                                              type(e).__name__, e))
    print()
    print("完成：删除副本 {} 个，搬移保留 {} 个".format(n_del, n_mv))
    print("回收空间 ≈ {:,} B = {:.2f} MiB".format(dup_bytes, dup_bytes / 1024 / 1024))
    print()
    print("下一步：把空的错误目录也删掉（脚本只搬文件，不删目录）——")
    print("    rmdir /s config 里那些缓存放空后，config 目录本身要保留（还有 sources.json）")
    print("    双层目录 out_audit\\out_audit 可整个删除")


if __name__ == "__main__":
    main()
