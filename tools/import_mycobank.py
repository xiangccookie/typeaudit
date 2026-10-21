#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""把 MycoBank 下载的 MBlist 压缩包 → 本工具能用的 mycobank_full.tsv

为什么是"你下载、我解析"
------------------------
MycoBank 是单页应用（SPA），**所有路由返回同一份空壳**（实测 `/Basic names search`
与首页字节数完全相同，7,506 B，HTML 里没有任何导出链接）⇒ 抓不到。
官方提供的是**打包下载**：去 MycoBank 的 Basic names search 页下载 **MBlist 压缩包**，
然后用本脚本解析成统一格式。

用法
----
  python tools/import_mycobank.py <MBlist.zip> [-o data/mycobank_full.tsv]
  python tools/import_mycobank.py <某个.txt>    [-o ...]     # 直接给解压后的文件也行

它做什么
--------
  ① 解压（zip / rar 不支持 → 请先手动解压）或直接读文件
  ② 在文件里找**表头**：必须同时含「Taxon name」与「Current name」两列
     （MBlist 的列名可能带引号、可能用逗号或制表符分隔、可能带 BOM —— 都做兼容）
  ③ 汇总所有命中的文件（MBlist 常按字母拆成多个文件）
  ④ 去重（同一 Taxon name 只留一条）后写出统一 TSV：
        Taxon name | Current name | Current name (bare) | MyCoBank # | Name status | Rank

⚠️ 若表头对不上，脚本**不猜**——会把每个文件的前 3 行打出来让你看，
   并提示需要改哪个常量。
"""
import argparse
import csv
import io
import os
import re
import sys
import zipfile

try:
    sys.stdout.reconfigure(encoding="utf-8")
    sys.stderr.reconfigure(encoding="utf-8")
except Exception:
    pass

# 本工具（sources.MycoBank）认的字段名 → 在 MBlist 里可能的写法（大小写/空格不敏感）
WANT = {
    "Taxon name":        ["taxonname", "name", "taxon"],
    "Current name":      ["currentname", "current", "acceptname", "acceptedname"],
    "Current name (bare)": ["currentnamebare", "currentbare"],
    "MycoBank #":        ["mycobank", "mycobankno", "mycobanknumber", "id"],
    "Name status":       ["namestatus", "status"],
    "Rank":              ["rank", "taxonrank"],
}
REQUIRED = ("Taxon name", "Current name")


def _norm(s):
    return re.sub(r"[\s_\-\.#]", "", (s or "").strip().strip('"')).lower()


def _sniff_delim(head):
    return max(["\t", ",", ";"], key=lambda d: head.count(d))


def _map_header(fields):
    """把文件表头映射到本工具的字段名；缺必需列则返回 None。"""
    got = {}
    for i, f in enumerate(fields or []):
        k = _norm(f)
        for canon, alts in WANT.items():
            if canon in got:
                continue
            if k == _norm(canon) or k in [_norm(a) for a in alts]:
                got[canon] = i
                break
    if all(r in got for r in REQUIRED):
        return got
    return None


def parse_stream(text, name, verbose=True):
    """解析一段文本；返回 (rows, header_map)。失败返回 ([], None)。"""
    lines = text.splitlines()
    if not lines:
        return [], None
    # 找表头行（可能前几行是说明）
    for start in range(min(8, len(lines))):
        head = lines[start]
        d = _sniff_delim(head)
        if d not in head:
            continue
        fields = next(csv.reader([head], delimiter=d))
        hm = _map_header(fields)
        if hm:
            rd = csv.reader(lines[start + 1:], delimiter=d)
            rows = []
            for r in rd:
                if not r or all(not x.strip() for x in r):
                    continue
                def g(canon):
                    i = hm.get(canon)
                    return (r[i].strip().strip('"') if (i is not None and i < len(r)) else "")
                if not g("Taxon name"):
                    continue
                rows.append([g("Taxon name"), g("Current name"), g("Current name (bare)"),
                             g("MycoBank #"), g("Name status"), g("Rank")])
            if verbose:
                print("   [{}] 表头在第 {} 行，分隔符 {!r}，{} 条".format(
                    name, start + 1, d, len(rows)))
            return rows, hm
    if verbose:
        print("   [{}] ⚠️ 找不到表头（需同时含 {}）。前 3 行：".format(
            name, " 与 ".join(REQUIRED)))
        for ln in lines[:3]:
            print("        " + ln[:140])
    return [], None


def iter_members(src):
    """产出 (成员名, 文本)。zip 就遍历成员；普通文件直接读。"""
    if zipfile.is_zipfile(src):
        with zipfile.ZipFile(src) as z:
            for info in z.infolist():
                if info.is_dir():
                    continue
                low = info.filename.lower()
                if low.endswith((".txt", ".csv", ".tsv", ".dat")) or "." not in os.path.basename(low):
                    data = z.read(info)
                    yield info.filename, data.decode("utf-8-sig", "replace")
    else:
        yield os.path.basename(src), io.open(src, encoding="utf-8-sig", errors="replace").read()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("src", help="MBlist 压缩包（.zip）或解压后的 .txt/.csv")
    ap.add_argument("-o", "--out", default=os.path.join("data", "mycobank_full.tsv"))
    a = ap.parse_args()

    if not os.path.exists(a.src):
        print("❌ 找不到 {}".format(a.src)); return 2
    print("读入：{}（{:,} B）".format(a.src, os.path.getsize(a.src)))

    all_rows, hm0, ok_files = [], None, 0
    for name, text in iter_members(a.src):
        rows, hm = parse_stream(text, name)
        if hm:
            ok_files += 1
            hm0 = hm0 or hm
            all_rows.extend(rows)
    print("\n命中文件 {} 个，累计 {} 条".format(ok_files, len(all_rows)))
    if not all_rows:
        print("❌ 一条都没解析出来。请把上面每个文件的「前 3 行」贴给我，我按真实表头调 WANT 映射。")
        return 1

    # 去重（同一 Taxon name 只留一条，保留首次出现）
    seen, uniq = set(), []
    for r in all_rows:
        k = r[0].lower()
        if k in seen:
            continue
        seen.add(k); uniq.append(r)
    print("去重后 {} 条（按 Taxon name）".format(len(uniq)))

    os.makedirs(os.path.dirname(os.path.abspath(a.out)) or ".", exist_ok=True)
    with io.open(a.out, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["Taxon name", "Current name", "Current name (bare)",
                    "MycoBank #", "Name status", "Rank"])
        w.writerows(uniq)
    print("\n✅ 已写出：{}（{:,} B）".format(a.out, os.path.getsize(a.out)))

    # 自检：用本工具的 MycoBank 适配器读回来，并抽查一条异名
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
    try:
        from typeaudit.sources import MycoBank
        mb = MycoBank(a.out)
        print("   适配器读回：{} 条".format(mb.n_rows()))
        for nm in ["Colletotrichum hymenocallidis", "Colletotrichum siamense"]:
            cur = mb.current_name(nm)
            print("   {:<34}→ {}".format(nm, cur or "（查不到）"))
        print("\n   把网页版的「MycoBank 异名表」指向这个文件即可。")
    except Exception as e:
        print("   （自检跳过：{}）".format(str(e)[:80]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
