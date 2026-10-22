#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""verify_missing_generic.py —— 对「号不存在」做 GenBank 反查（任意属通用）。

【为什么要重写】
    旧版 verify_missing.py 用 '"CBS 126514"'（带引号的空格形式）检索。
    实测这种形式会被 NCBI 分词成 "CBS" AND "126514"，命中大量无关记录
    （本次会话实测：检索 "CBS 112996" 返回 20,561 条假命中）。
    改为 '"CBS:126514"' 的【冒号形式】（不引号）后精确，再逐条核对记录原文。

【判定】
    有记录且原文含该号、且 /organism 种加词相符 → 号可能真实，机构未收录
    有记录但种名不符                              → 需人工看
    零记录                                        → 支持「号不存在」
"""
import argparse, io, os, re, sys, time

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CACHE = os.path.join(REPO, ".cache", "verify_missing2")


def norm_num(s):
    return re.sub(r"[^A-Za-z0-9]", "", (s or "")).upper()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True, help="如 out_diaporthe")
    ap.add_argument("--genus", required=True, help="如 diaporthe")
    ap.add_argument("--email", default=os.environ.get("NCBI_EMAIL") or "audit@example.org")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--start", type=int, default=0, help="从第几条开始（配合缓存断续跑）")
    a = ap.parse_args()

    from Bio import Entrez
    Entrez.email = a.email
    os.makedirs(CACHE, exist_ok=True)

    src = os.path.join(REPO, a.out, "{}_发现汇总.tsv".format(a.genus))
    L = [x.split("\t") for x in io.open(src, encoding="utf-8-sig").read().splitlines() if x.strip()]
    h = L[0]; ci, oi = h.index("类别"), h.index("对象")
    items = []
    for r in L[1:]:
        if not r[ci].startswith("号不存在"):
            continue
        m = re.match(r"^(.*?)\s*/\s*(.+)$", r[oi].strip())
        if m:
            items.append((m.group(1).strip(), m.group(2).strip()))
    items = items[a.start:]
    if a.limit:
        items = items[:a.limit]
    print("待查 {} 条（从第 {} 条起）".format(len(items), a.start), flush=True)

    # ⚠️ 逐条落盘：单条查询可能拖很久（NCBI 限流），
    #    原来只在末尾一次性写文件，中途被超时杀掉就全丢（实测丢过一次 43 条）。
    op = os.path.join(REPO, a.out, "反向审计_号不存在.tsv")
    first = not os.path.exists(op) or a.start == 0
    fh = io.open(op, "w" if first else "a", encoding="utf-8-sig", newline="")
    if first:
        fh.write("号\t声明种\tGenBank命中数\tGenBank里的种名\t判定\n")
        fh.flush()

    out = []
    for i, (sp, num) in enumerate(items, 1):
        t_item = time.time()
        key = norm_num(num)[:64] or "x"
        cp = os.path.join(CACHE, key + ".txt")
        if os.path.exists(cp):
            cnt, orgs = io.open(cp, encoding="utf-8").read().split("\n", 1)
        else:
            cnt, orgs = "0", ""
            code = re.match(r"^([A-Za-z][A-Za-z\.\-]*)", num.strip())
            digit = re.search(r"([\d][\w\.\-]*)$", num.strip())
            term = "{}:{}".format(code.group(1), digit.group(1)) if (code and digit) else '"{}"'.format(num)
            try:
                hh = Entrez.esearch(db="nuccore", term=term, retmax=20)
                ids = Entrez.read(hh)["IdList"]; hh.close(); time.sleep(0.35)
                cnt = str(len(ids))
                if ids:
                    # ⚠️ 必须先用 esummary 看记录长度：实测有些号命中的是【基因组级】记录
                    #    （CBS 102804 单条 66 MB、CBS 120623 36 MB），efetch 会 IncompleteRead，
                    #    一条就拖十几分钟，把整批拖超时。
                    #    这里只取 ≤200 kb 的记录做逐条核对；其余记下"体量过大"。
                    nm = norm_num(num)
                    small, big = [], 0
                    try:
                        ss = Entrez.esummary(db="nuccore", id=",".join(ids))
                        sums = Entrez.read(ss); ss.close(); time.sleep(0.35)
                        for s_ in sums:
                            try:
                                ln = int(s_.get("Length") or 0)
                            except Exception:
                                ln = 0
                            if 0 < ln <= 200000:
                                small.append(str(s_["Id"]))
                            else:
                                big += 1
                    except Exception:
                        small = list(ids)
                    keep = []
                    if small:
                        nmax = min(len(small), 8)
                        hh = Entrez.efetch(db="nuccore", id=",".join(small[:nmax]),
                                           rettype="gb", retmode="text")
                        txt = hh.read(); hh.close(); time.sleep(0.35)
                        for rec in txt.split("\n//\n"):
                            if not rec.strip():
                                continue
                            if nm not in norm_num(re.sub(r"[^A-Za-z0-9]", "", rec)[:20000]):
                                continue
                            mo = re.search(r'/organism="([^"]+)"', rec)
                            if mo:
                                keep.append(mo.group(1))
                    cnt = str(len(keep))
                    orgs = ",".join(sorted(set(keep)))
                    if not keep and big:
                        cnt = "-2"        # -2 = 有记录但都是大记录，未逐条核对
                        orgs = "（命中 {} 条体量过大记录，未核对）".format(big)
            except Exception as e:
                cnt, orgs = "-1", "查询失败:" + str(e)[:40]
            io.open(cp, "w", encoding="utf-8").write("{}\n{}".format(cnt, orgs))

        n = int(cnt) if str(cnt).lstrip("-").isdigit() else -1
        epi = sp.split()[-1].lower() if sp else ""
        ok = epi and epi in (orgs or "").lower()
        if n == -2:
            v = "GenBank 有引用但均为基因组级大记录，未逐条核对 → 需人工看"
        elif n < 0:
            v = "查询失败"
        elif n == 0:
            v = "GenBank 亦无引用 → 支持「号不存在」"
        elif ok:
            v = "GenBank 有 {} 条记录引用该号且种名相符 → 号可能真实，机构未收录".format(n)
        else:
            v = "GenBank 有 {} 条记录引用该号但种名不符 → 需人工看".format(n)
        out.append({"号": num, "声明种": sp, "GenBank命中数": cnt,
                    "GenBank里的种名": (orgs or "")[:90], "判定": v})
        fh.write("{}\t{}\t{}\t{}\t{}\n".format(num, sp, cnt, (orgs or "")[:90], v))
        fh.flush()
        print("  [{}/{}] {:<16} {:>3} 条 ｜ {:.0f}s".format(
            i, len(items), num[:14], cnt if str(cnt).lstrip("-").isdigit() else -2,
            time.time() - t_item), flush=True)

    fh.close()
    print("完成 → {} （本批 {} 条）".format(op, len(out)), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
