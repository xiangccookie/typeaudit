#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""verify_missing.py —— 对「号不存在」做 GenBank 反查。

【为什么这样验证】
    机构目录查无 ≠ 这个号是编造的。有两种可能：
      ① 号是真的，只是该机构没把它公开上线（新种、私人保藏常见）
      ② 号写错了/编造
    区分办法：去 GenBank 搜这个号。有人引用过 → 说明这号在实物世界存在；
              零引用 → 支持"号不存在"的判断（仍非铁证，但强得多）。

输出：out_final/核对1_号不存在_GenBank反查.csv
"""
import csv, io, os, re, sys, time

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CACHE = os.path.join(REPO, ".cache", "verify_missing")


def main():
    from Bio import Entrez
    email = os.environ.get("NCBI_EMAIL") or "audit@example.org"
    Entrez.email = email
    p = os.path.join(REPO, "out_final", "核对1_号不存在.csv")
    rows = list(csv.DictReader(io.open(p, encoding="utf-8-sig")))
    print("待查 {} 条".format(len(rows)))
    os.makedirs(CACHE, exist_ok=True)
    out = []
    for i, r in enumerate(rows, 1):
        num = (r.get("号") or "").strip()
        sp = (r.get("种名") or "").strip()
        key = re.sub(r"[^A-Za-z0-9]", "_", num)
        cp = os.path.join(CACHE, key + ".txt")
        if os.path.exists(cp):
            c, orgs = io.open(cp, encoding="utf-8").read().split("\n", 1)
        else:
            c, orgs = "0", ""
            try:
                h = Entrez.esearch(db="nuccore", term='"{}"'.format(num), retmax=20)
                c = Entrez.read(h)["Count"]; h.close()
                time.sleep(0.4)
                if int(c) > 0:
                    h = Entrez.esearch(db="nuccore", term='"{}"'.format(num), retmax=20)
                    ids = Entrez.read(h)["IdList"]; h.close()
                    time.sleep(0.4)
                    if ids:
                        h = Entrez.efetch(db="nuccore", id=",".join(ids), rettype="gb", retmode="text")
                        txt = h.read(); h.close()
                        orgs = ",".join(sorted({m.group(1) for m in
                                                re.finditer(r'/organism="([^"]+)"', txt)}))
                        time.sleep(0.4)
            except Exception as e:
                c, orgs = -1, "查询失败:" + str(e)[:40]
            io.open(cp, "w", encoding="utf-8").write("{}\n{}".format(c, orgs))
        n = int(c) if str(c).lstrip("-").isdigit() else -1
        sp_bare = sp.split()[-1] if sp else ""
        org_has_sp = sp_bare and sp_bare.lower() in (orgs or "").lower()
        if n < 0:
            verdict = "查询失败"
        elif n == 0:
            verdict = "GenBank 也查无 → 支持「号不存在」"
        elif org_has_sp:
            verdict = "GenBank 有 {} 条记录引用该号，且种名相符 → 号可能真实，机构未收录".format(n)
        else:
            verdict = "GenBank 有 {} 条记录引用该号，但种名不符 → 需人工看".format(n)
        out.append({"号": num, "声明种": sp, "机构": (r.get("查无于") or ""),
                    "GenBank命中数": c, "GenBank里的种名": (orgs or "")[:90],
                    "判定": verdict})
        if i % 10 == 0:
            print("  [{}/{}]".format(i, len(rows)), flush=True)

    op = os.path.join(REPO, "out_final", "核对1_号不存在_GenBank反查.csv")
    with io.open(op, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["号", "声明种", "机构", "GenBank命中数",
                                          "GenBank里的种名", "判定"])
        w.writeheader(); w.writerows(out)
    from collections import Counter
    print("\n结果 → {}".format(os.path.relpath(op, REPO)))
    for k, v in Counter(x["判定"] for x in out).most_common():
        print("  {:<60} {}".format(k[:58], v))
    return 0


if __name__ == "__main__":
    sys.exit(main())
