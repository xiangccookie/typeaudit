#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
setup_data.py —— 一条命令装齐外部数据，并留下溯源。

【为什么做成「安装器」而不是「把数据打包进仓库」】
    CBS 目录 2.0 MB、MycoBank 全表 58.3 MB，塞进仓库毫无技术障碍。但不打包，三个理由：
      ① 许可：两者都是第三方数据库（Westerdijk Institute / CBS-KNAW），再分发须先确认条款；
      ② 时效：**本机实测 CBS 会变** —— status 会被改、号会加、编号还有两套写法
         （CBS 100.06 与 CBS 10006 是两株不同的菌）。打包一份快照，用户会把它
         当成现况用，**这比让他自己去下更危险**；
      ③ 取数工具已经在仓库里，缺的只是「一条命令把顺序跑完 + 记下溯源」。
    所以这里做的是：体检 → 给指引 → 能自动的自动 → 记溯源。

用法
    python tools/setup_data.py                              # 体检：三个数据集在不在、多大、从哪来
    python tools/setup_data.py --fetch-cbs                  # 拉 CBS 全量（约 30 次请求 / 10–40 分钟）
    python tools/setup_data.py --fetch-cbs --genus Colletotrichum   # 只拉一个属
    python tools/setup_data.py --import-mycobank <MBlist.zip>       # 解析官方 MycoBank 压缩包
    python tools/setup_data.py --provenance                 # 只看溯源记录
"""
import sys
try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

import argparse, hashlib, io, json, os, subprocess, time

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROV = os.path.join(REPO, "data", "_provenance.json")

# 三个数据集：显示名 → 候选路径（与 sources._FALLBACK_PATHS 保持一致）
DATASETS = [
    ("CBS 目录", "CBSCatalogue",
     ["data/cbs_catalogue.tsv",
      r"C:\Users\Administrator\Desktop\DeepSeek\bio原理\out\cbs_catalogue.tsv"],
     "CBS 官方后端 BioAware / BioloMICS REST API（tools/fetch_cbs_live.py 拉取）"),
    ("MycoBank 全表", "MycoBank",
     ["data/mycobank_full.tsv",
      r"C:\Users\Administrator\Desktop\Tool\脚本\调取mycobank数据，分析sp_curret,编码，是否为复合体\mycobank_full.tsv"],
     "MycoBank 官方 MBlist 压缩包（tools/import_mycobank.py 解析）"),
    ("机构登记表", "Biocollections",
     ["data/biocollections_codes.csv"],
     "NCBI Biocollections 登记表（仓库自带精简版，开箱即用）"),
]


def resolve(cands):
    for rel in cands:
        p = rel if os.path.isabs(rel) else os.path.join(REPO, rel)
        if os.path.exists(p):
            return p
    return None


def sha256(path, limit=None):
    """算指纹；limit 非空时只算前 N 字节（大文件全量算太慢，够做变更检测）。"""
    h = hashlib.sha256()
    with io.open(path, "rb") as f:
        if limit:
            h.update(f.read(limit))
        else:
            for chunk in iter(lambda: f.read(1 << 20), b""):
                h.update(chunk)
    return h.hexdigest()


def count_rows(path):
    try:
        with io.open(path, encoding="utf-8", errors="replace") as f:
            return sum(1 for _ in f)
    except Exception:
        return None


def load_prov():
    if os.path.exists(PROV):
        try:
            return json.loads(io.open(PROV, encoding="utf-8").read())
        except Exception:
            return {}
    return {}


def save_prov(d):
    if not os.path.isdir(os.path.dirname(PROV)):
        os.makedirs(os.path.dirname(PROV))
    with io.open(PROV, "w", encoding="utf-8") as f:
        json.dump(d, f, ensure_ascii=False, indent=1)


def record(name, path, source, how):
    """把「这份数据从哪来、什么时候来的」写进溯源文件。"""
    prov = load_prov()
    prov[name] = {
        "path": path,
        "size": os.path.getsize(path),
        "mtime": time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(os.path.getmtime(path))),
        "rows": count_rows(path),
        "sha256_first1MB": sha256(path, 1 << 20),
        "source": source,
        "how": how,
        "recorded_at": time.strftime("%Y-%m-%d %H:%M:%S"),
    }
    save_prov(prov)
    print("  已记溯源：{}".format(name))
    return prov[name]


def show():
    prov = load_prov()
    print("=" * 100)
    print("外部数据体检")
    print("=" * 100)
    print("  {:<16} {:<8} {:>12} {:<18} {}".format("数据集", "状态", "大小", "修改时间", "路径"))
    print("  " + "-" * 94)
    n_ok = 0
    for label, adapter, cands, source in DATASETS:
        p = resolve(cands)
        if p:
            n_ok += 1
            st = "✅"
            size = "{:,} B".format(os.path.getsize(p))
            mt = time.strftime("%Y-%m-%d %H:%M", time.localtime(os.path.getmtime(p)))
        else:
            st, size, mt = "❌ 缺", "-", "-"
        print("  {:<16} {:<8} {:>12} {:<18} {}".format(
            label, st, size, mt, (p.replace(REPO + os.sep, "") if p else "（没找到）")))
    print()
    if prov:
        print("  溯源记录（data/_provenance.json）：")
        for k, v in prov.items():
            print("     {:<16} 来源 {}　{}".format(k, v.get("source", "?")[:44], v.get("recorded_at", "")))
            print("                      sha256(前1MB)={}…".format((v.get("sha256_first1MB") or "")[:16]))
    else:
        print("  还没有溯源记录。跑一次 --fetch-cbs 或 --import-mycobank 就会生成。")
    print()
    print("  缺失项怎么补：")
    print("     CBS 目录      python tools/setup_data.py --fetch-cbs [--genus 属名]")
    print("     MycoBank 全表 去官网下 MBlist 压缩包，然后：")
    print("                   python tools/setup_data.py --import-mycobank <下载的.zip>")
    print("       官网：https://www.mycobank.org/Basic names search")
    print("     机构登记表     仓库自带，无需处理")
    print()
    print("  提示：路径不存在时，工具会按上表候选位置 + ~/.dsh/typeaudit_refs.json 自动回退；")
    print("        所以「放在这些位置」就等于装好了。")
    return 0 if n_ok == len(DATASETS) else 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fetch-cbs", action="store_true", help="拉 CBS 目录（走 fetch_cbs_live.py）")
    ap.add_argument("--genus", default=None, help="配合 --fetch-cbs：只拉该属")
    ap.add_argument("--import-mycobank", default=None, help="解析 MycoBank MBlist 压缩包")
    ap.add_argument("--out", default=None, help="产物路径（默认 data/ 下）")
    ap.add_argument("--email", default=os.environ.get("NCBI_EMAIL") or "", help="写进 User-Agent 的联系邮箱")
    ap.add_argument("--provenance", action="store_true", help="只看溯源记录")
    a = ap.parse_args()

    if a.provenance:
        show()
        return 0

    if a.fetch_cbs:
        cmd = [sys.executable, "-u", os.path.join(REPO, "tools", "fetch_cbs_live.py")]
        if a.genus:
            cmd += ["--genus", a.genus]
        if a.out:
            cmd += ["--out-dir", a.out]
        if a.email:
            cmd += ["--email", a.email]
        print("$ " + " ".join(cmd[2:]))
        rc = subprocess.call(cmd, cwd=REPO)
        # 抓到产物路径就记溯源
        od = a.out or os.path.join(REPO, "out")
        cands = []
        if os.path.isdir(od):
            cands = [os.path.join(od, f) for f in sorted(os.listdir(od))
                     if f.startswith("cbs_") and f.endswith(".tsv") and "diff" not in f]
        if cands:
            newest = max(cands, key=os.path.getmtime)
            record("CBS 目录（本次拉取）", newest,
                   "CBS 官方后端 BioAware REST API",
                   "tools/fetch_cbs_live.py {}".format("--genus " + a.genus if a.genus else "（全量）"))
        print()
        show()
        return rc

    if a.import_mycobank:
        zip_p = a.import_mycobank
        if not os.path.exists(zip_p):
            print("找不到压缩包：{}".format(zip_p))
            return 1
        dst = a.out or os.path.join(REPO, "data", "mycobank_full.tsv")
        cmd = [sys.executable, "-u", os.path.join(REPO, "tools", "import_mycobank.py"),
               zip_p, "-o", dst]
        print("$ " + " ".join(cmd[2:]))
        rc = subprocess.call(cmd, cwd=REPO)
        if rc == 0 and os.path.exists(dst):
            record("MycoBank 全表", dst, "MycoBank 官方 MBlist",
                   "tools/import_mycobank.py " + os.path.basename(zip_p))
        print()
        show()
        return rc

    return show()


if __name__ == "__main__":
    sys.exit(main())
