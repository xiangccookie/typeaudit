#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""verify_typematerial.py —— 抓 GenBank 记录的 /type_material 原文，判定类型等级。

【为什么这一步能自动做】
    「机构记为次要模式」的争议点是：记录自己说的是哪一级？
    如果记录原文写 "culture from epitype of ..." → 记录准确，非错；
    如果写 "type material" 或 "holotype" → 与机构记的次要模式冲突。
    这些文字就在 GenBank 记录里，不用查论文。

输出：out_final/核对4b_类型材料原文.csv
"""
import csv, io, os, re, sys, time

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CACHE = os.path.join(REPO, ".cache", "typematerial")


def main():
    from Bio import Entrez
    email = os.environ.get("NCBI_EMAIL") or "audit@example.org"
    Entrez.email = email
    p = os.path.join(REPO, "out_final", "核对4_需查原始发表.csv")
    rows = [r for r in csv.DictReader(io.open(p, encoding="utf-8-sig"))
            if r["类别"].startswith("模式地位")]
    print("待查 {} 条".format(len(rows)))
    os.makedirs(CACHE, exist_ok=True)
    out = []
    for i, r in enumerate(rows, 1):
        num, sp = r["号"].strip(), r["声明种"].strip()
        key = re.sub(r"[^A-Za-z0-9]", "_", num)
        cp = os.path.join(CACHE, key + ".txt")
        if os.path.exists(cp):
            body = io.open(cp, encoding="utf-8").read()
        else:
            body = ""
            try:
                h = Entrez.esearch(db="nuccore", term='"{}"'.format(num), retmax=30)
                ids = Entrez.read(h)["IdList"]; h.close()
                time.sleep(0.4)
                if ids:
                    h = Entrez.efetch(db="nuccore", id=",".join(ids), rettype="gb", retmode="text")
                    body = h.read(); h.close()
                    time.sleep(0.4)
            except Exception as e:
                body = "ERROR:" + str(e)[:60]
            io.open(cp, "w", encoding="utf-8").write(body)
        # 抽所有 /type_material 原文 + 该记录的 organism
        recs = re.split(r"\n//\n", body) if body and not body.startswith("ERROR") else []
        hits = []
        for rec in recs:
            org = re.search(r'/organism="([^"]+)"', rec)
            tm = re.findall(r'/type_material="([^"]+)"', rec)
            if tm and org:
                hits.append((org.group(1), "; ".join(tm)))
        sp_bare = sp.split()[-1].lower()
        mine = [h for h in hits if sp_bare in h[0].lower()]
        texts = " | ".join(t for _, t in (mine or hits))[:160]
        low = texts.lower()
        if not hits:
            verdict = "该号的记录里没有 /type_material 字段"
        elif "epitype" in low:
            verdict = "记录原文也写 epitype → 与机构一致，措辞准确（非错）"
        elif re.search(r"\bholotype\b|\bneotype\b|\blectotype\b", low):
            verdict = "记录称正模/新模式/选模式，机构记为次要模式 → 冲突"
        elif "type material" in low:
            verdict = "记录只笼统写 type material，未标等级 → 不精确"
        else:
            verdict = "有类型字段但等级不明"
        out.append({"号": num, "声明种": sp, "机构记的类型": "次要模式",
                    "记录的/type_material原文": texts or "（无）",
                    "记录里出现的种名": " | ".join(sorted({o for o, _ in hits}))[:70],
                    "判定": verdict})
        if i % 5 == 0:
            print("  [{}/{}]".format(i, len(rows)), flush=True)

    op = os.path.join(REPO, "out_final", "核对4b_类型材料原文.csv")
    with io.open(op, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["号", "声明种", "机构记的类型",
                                          "记录的/type_material原文", "记录里出现的种名", "判定"])
        w.writeheader(); w.writerows(out)
    from collections import Counter
    print("\n结果 → {}".format(os.path.relpath(op, REPO)))
    for k, v in Counter(x["判定"] for x in out).most_common():
        print("  {:<52} {}".format(k[:50], v))
    return 0


if __name__ == "__main__":
    sys.exit(main())
