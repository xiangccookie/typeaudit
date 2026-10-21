#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""参考数据的「获取与维护」——不假装能一键下载，如实分工。

三个文件的现实（都实测过）
--------------------------
  机构登记表（9,672 行 / 152 KB）   ✅ **已打包进本仓库** data/biocollections_codes.csv，开箱即用
  MycoBank 全库（40 万条）          ⚠️ **官方提供 MBlist 压缩包，你下、我解析**
                                       → tools/import_mycobank.py
  CBS 目录（6.3 万号）              ❌ **抓不了**：wi.knaw.nl 已改 Angular SPA
                                       （实测与 MycoBank 共用同一份 main.*.js；HTML 只是 7,506 字符空壳，
                                        去标签后仅剩 "Loading..."）。你现有那份是改版前抓的，够用。

用法
----
  python tools/fetch_refdata.py guide                        # 打印三个文件的获取步骤（含下载地址）
  python tools/fetch_refdata.py cbs --check                  # 报告现有 CBS 目录覆盖到哪些号段
  python tools/fetch_refdata.py cbs --add CBS 130343 "holotype of X"   # 手工补录一条（从网站查到后）
  python tools/fetch_refdata.py mycobank --import <MBlist.zip>         # 转调导入器

为什么保留「手工补录」
----------------------
既然自动抓不了，**增量补充的现实路径就是手工**：你在 CBS 网站上查到某个号，
用 `--add` 补进目录。脚本负责把它归成规范键（去点/去空格）并与现有目录合并——
**这一步才是最容易出错的**（norm_key 规则不一致会让号"存在却查不到"）。
"""
import argparse
import csv
import io
import os
import re
import sys
import time

try:
    sys.stdout.reconfigure(encoding="utf-8")
    sys.stderr.reconfigure(encoding="utf-8")
except Exception:
    pass

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(REPO, "src"))

MYCOBANK_URL = "https://www.mycobank.org/Basic%20names%20search"
CBS_URL_PAT = "https://wi.knaw.nl/page/fungal_display/name/CBS%20{num}"


def norm_key(code, number):
    """必须与 sources.norm_key 完全一致的归一键（含去点！）。
    目录生成时把点也去了（实测 CBS 1546.1 → CBS15461），不一致会让号"存在却查不到"。"""
    try:
        from typeaudit.sources import norm_key as _nk
        return _nk(code, number)
    except Exception:
        return re.sub(r"[\s:：.\-]", "", str(code) + str(number)).upper()


def cbs_path():
    """找 CBS 目录：先 registry，再仓库 config。"""
    rp = os.path.join(os.path.expanduser("~"), ".dsh", "typeaudit_refs.json")
    try:
        if os.path.exists(rp):
            import json
            d = json.load(io.open(rp, encoding="utf-8"))
            p = (d.get("cbs") or {}).get("path")
            if p and os.path.exists(p):
                return p
    except Exception:
        pass
    p = os.path.join(REPO, "data", "cbs_catalogue.tsv")
    return p if os.path.exists(p) else None


def load_cbs(path):
    rows = []
    if path and os.path.exists(path):
        with io.open(path, encoding="utf-8-sig", newline="", errors="replace") as f:
            rows = list(csv.DictReader(f, delimiter="\t"))
    return rows


def cmd_guide():
    print("=" * 88)
    print("三个参考数据的获取步骤")
    print("=" * 88)
    print("""
① 机构登记表 ★ 已随仓库提供，无需任何操作
     data/biocollections_codes.csv   （9,672 行 / 152 KB，由 NCBI Biocollections 精简而来）
     若要与最新完整版同步：NCBI Biocollections（原 GrBio）下载完整 CSV，替换 config 里的路径。

② MycoBank 异名表  ← 官方打包下载 + 本工具解析
     1) 打开 {url}
     2) 找到「MBList」下载入口，下载压缩包（通常几十 MB）
     3) 解析成本工具格式：
          python tools/import_mycobank.py <下载的.zip> -o data/mycobank_full.tsv
        （脚本会自动找表头、兼容逗号/制表符、汇总多文件、按 Taxon name 去重）
     4) 把网页版「MycoBank 异名表」指向 data/mycobank_full.tsv

③ CBS 目录  ← ❌ 无法自动抓，用你手上那份
     wi.knaw.nl 已改 Angular SPA（与 MycoBank 共用同一份构建产物），
     HTML 抓不到数据。你现有那份（改版前抓的）继续用即可。
     需要补新号时：
       · 在网站查到该号的状态
       · python tools/fetch_refdata.py cbs --add CBS 130343 "holotype of X"
""".format(url=MYCOBANK_URL))
    p = cbs_path()
    print("当前 CBS 目录：{}".format(p or "❌ 未找到"))
    return 0


def cmd_cbs_check():
    p = cbs_path()
    if not p:
        print("❌ 未找到 CBS 目录（先跑 find_refdata.py 或 --add 手工登记）")
        return 2
    rows = load_cbs(p)
    print("CBS 目录：{}".format(p))
    print("  行数：{:,}".format(len(rows)))
    nums = []
    for r in rows:
        n = (r.get("cbs_number") or "").replace("CBS", "").strip()
        m = re.match(r"^(\d+)", n)
        if m:
            nums.append(int(m.group(1)))
    if nums:
        nums.sort()
        print("  号段范围：CBS {} – {}（中位 {}）".format(nums[0], nums[-1], nums[len(nums)//2]))
        # 直方图（每 1 万一个桶）
        from collections import Counter
        c = Counter((x // 10000) * 10000 for x in nums)
        print("  分布（每 1 万号一桶，前 12 桶）：")
        for k in sorted(c)[:12]:
            print("      CBS {:>6}–{:<6}{:>6,}  {}".format(
                k, k + 9999, c[k], "█" * min(60, c[k] // 200)))
    st = sum(1 for r in rows if (r.get("status") or "").strip())
    print("  有状态标注：{:,} 行（{:.1f}%）".format(st, 100.0 * st / max(1, len(rows))))
    print("\n  ⚠️ 没有标注的行查出来会是「未核」（不是「号不存在」）——这是正确的：没标注就不能判。")
    return 0


def cmd_cbs_add(code, num, status):
    p = cbs_path()
    if not p:
        p = os.path.join(REPO, "data", "cbs_catalogue.tsv")
        os.makedirs(os.path.dirname(p), exist_ok=True)
    rows = load_cbs(p)
    key = norm_key(code, num)
    hit = next((r for r in rows if (r.get("norm_key") or "") == key), None)
    if hit:
        print("ℹ️ 已存在，更新状态：{} {}  {!r} → {!r}".format(
            code, num, (hit.get("status") or "")[:40], status[:40]))
        hit["status"] = status
        hit["has_status"] = "1" if status.strip() else ""
        hit["cbs_number"] = "{} {}".format(code, num)
    else:
        print("➕ 新增：{} {}  {!r}".format(code, num, status[:50]))
        rows.append({"cbs_number": "{} {}".format(code, num), "norm_key": key,
                     "status": status, "has_status": "1" if status.strip() else ""})
    rows.sort(key=lambda r: (r.get("norm_key") or ""))
    with io.open(p, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["cbs_number", "norm_key", "status", "has_status"],
                           delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print("   已写回 {}（{:,} 行）".format(p, len(rows)))
    print("   归一键 = {}（与 sources.norm_key 一致：**去点**）".format(key))
    return 0


def cmd_mycobank_import(zip_path):
    imp = os.path.join(HERE, "import_mycobank.py")
    if not os.path.exists(imp):
        print("❌ 找不到 import_mycobank.py")
        return 2
    import subprocess
    out = os.path.join(REPO, "data", "mycobank_full.tsv")
    r = subprocess.run([sys.executable, imp, zip_path, "-o", out], cwd=REPO)
    return r.returncode


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd")
    sub.add_parser("guide")
    c = sub.add_parser("cbs")
    c.add_argument("--check", action="store_true")
    c.add_argument("--add", nargs=3, metavar=("CODE", "NUM", "STATUS"))
    m = sub.add_parser("mycobank")
    m.add_argument("--import", dest="imp", metavar="MBLIST_ZIP")
    a = ap.parse_args()
    if a.cmd == "guide" or not a.cmd:
        return cmd_guide()
    if a.cmd == "cbs":
        if a.add:
            return cmd_cbs_add(*a.add)
        return cmd_cbs_check()
    if a.cmd == "mycobank":
        if a.imp:
            return cmd_mycobank_import(a.imp)
        return cmd_guide()
    return 0


if __name__ == "__main__":
    sys.exit(main())
