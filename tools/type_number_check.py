# -*- coding: utf-8 -*-
"""type_number_check.py —— 反查：一个种声明的「模式材料号」，在 GenBank 里被标成了哪些物种。

为什么要这一层
--------------
正向审计是「从记录里的号 → 查机构目录」；这一层反过来：
**拿该种权威侧声明的模式材料号，去 GenBank 里看它出现在哪些 organism 下。**
若同一个号在 GenBank 里被写成两个物种，就是数据库自己写下的矛盾（直接证据，不依赖判定规则）。

与 tools/seq_axis.py 的区别
---------------------------
seq_axis 的 same_number 只检查「三列表里已有的号」，而三列表的号受 /strain 扫描上限约束；
本脚本直接以**类型号清单**为输入，不受该上限影响。

用法
----
    python tools/type_number_check.py --out out_colleto_web --genus colletotrichum --limit 0
"""
import argparse
import io
import json
import os
import re
import time
from collections import Counter

from Bio import Entrez

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CACHE = os.path.join(REPO, ".cache", "typecheck")
os.makedirs(CACHE, exist_ok=True)


def log(m):
    print(m, flush=True)


def key_of(display):
    """'CBS 128530' → ('CBS','128530')；'COAD 1877' → ('COAD','1877')"""
    m = re.match(r"^\s*([A-Za-z][A-Za-z\.\-]*)\s*[:#]?\s*([\w\.\-]+)\s*$", display)
    if not m:
        return None
    return m.group(1).upper(), m.group(2)


def binomials_of(ids):
    """批量 esummary → 每条记录 Title 的前两个词（= 物种名）"""
    names = Counter()
    for i in range(0, len(ids), 100):
        chunk = ids[i:i + 100]
        try:
            h = Entrez.esummary(db="nucleotide", id=",".join(chunk))
            s = Entrez.read(h)
            h.close()
            for d in s:
                title = str(d.get("Title", ""))
                words = title.split()
                if len(words) >= 2:
                    names["{} {}".format(words[0], words[1].rstrip(",.;"))] += 1
        except Exception as e:
            log("   ! esummary 失败: {}".format(str(e)[:60]))
        time.sleep(0.36)
    return names


def check_one(code, num):
    """返回 (记录数, 物种集合 Counter)"""
    cp = os.path.join(CACHE, "{}_{}.json".format(code, num))
    if os.path.exists(cp):
        d = json.load(io.open(cp, encoding="utf-8"))
        return d["n"], Counter(d["names"])
    ids = []
    for term in ('"{}:{}"'.format(code, num), '"{} {}"'.format(code, num)):
        try:
            h = Entrez.esearch(db="nucleotide", term=term, retmax=100)
            r = Entrez.read(h)
            h.close()
            ids = r.get("IdList", [])
            time.sleep(0.36)
        except Exception:
            ids = []
        if ids:
            break
    names = binomials_of(ids) if ids else Counter()
    json.dump({"n": len(ids), "names": dict(names)},
              io.open(cp, "w", encoding="utf-8"), ensure_ascii=False)
    return len(ids), names


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--genus", required=True)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--email", default=os.environ.get("NCBI_EMAIL", "audit@example.org"))
    a = ap.parse_args()
    Entrez.email = a.email

    od = os.path.join(REPO, a.out)
    ev = os.path.join(od, "{}_证据.jsonl".format(a.genus))
    items, seen = [], set()
    for line in io.open(ev, encoding="utf-8", errors="replace"):
        d = json.loads(line)
        for k, v in (d.get("numbers") or {}).items():
            if not v.get("claimed_type"):
                continue
            disp = v.get("display") or k
            kn = key_of(disp)
            if not kn:
                continue
            sig = (d["species"], kn[0] + kn[1])
            if sig in seen:
                continue
            seen.add(sig)
            items.append((d["species"], disp, kn))
    if a.limit:
        items = items[:a.limit]
    log("待查类型号 {} 个（覆盖 {} 个种）".format(len(items), len({i[0] for i in items})))

    rows = []
    for n, (sp, disp, kn) in enumerate(items, 1):
        cnt, names = check_one(kn[0], kn[1])
        claim = sp.split()[-1].lower()
        homonyms = {x: c for x, c in names.items() if x.split()[-1].lower() != claim}
        verdict = ("未检出记录" if cnt == 0
                   else ("两个及以上物种" if len(names) > 1
                         else ("记为别的种" if homonyms else "一致")))
        rows.append({
            "种": sp, "类型号": disp, "记录数": cnt,
            "GenBank 里的物种": " / ".join("{}×{}".format(k, v) for k, v in names.most_common(4)),
            "判定": verdict,
        })
        if n % 50 == 0 or n == len(items):
            log("   [{}/{}] {}".format(n, len(items), disp))

    tsv = os.path.join(od, "{}_类型号反查.tsv".format(a.genus))
    cols = ["种", "类型号", "记录数", "GenBank 里的物种", "判定"]
    with io.open(tsv, "w", encoding="utf-8-sig", newline="") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")

    c = Counter(r["判定"] for r in rows)
    log("\n== 汇总 ==")
    for k, v in c.most_common():
        log("   {:<14} {}".format(k, v))
    log("   明细 → {}".format(tsv))
    json.dump({"属": a.genus, "类型号数": len(rows), "判定分布": dict(c), "明细": tsv},
              io.open(os.path.join(od, "{}_类型号反查_汇总.json".format(a.genus)), "w",
                      encoding="utf-8"), ensure_ascii=False, indent=2)


if __name__ == "__main__":
    main()
