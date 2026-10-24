# -*- coding: utf-8 -*-
"""genome_scan.py —— 基因组层扫描（三属）：模式菌株有没有基因组、多少组装被官方判为归属错误。

v2 修正（第一版把 0 判成"无模式声明"，根因见下）
------------------------------------------------
⚠️ `type-material` **不在** `assembly_info.type_material` 里（该字段恒为 None），
   它在 `assembly_info.biosample.attributes` 里，键名为 `type-material`。
   第一版只看前者，于是 97 个有基因组的种全部被判成「无模式声明」——错得干净利落。

数据来源
--------
1. esearch db=taxonomy  `"Genus species"`  → taxid
2. Datasets API v2alpha `genome/taxon/{taxid}/dataset_report` → 该种**全部**组装
      assembly_info.biosample.attributes['type-material']  → 记录自带的模式材料声明
      assembly_info.genome_notes                          → 含 'unverified source organism' 即官方判归属错误

判据
----
有模式菌株基因组：该种**任一**组装的 type-material 非空
官方判错：genome_notes 里出现 'unverified source organism'

用法
----
    python tools/genome_scan.py --genus Colletotrichum --out out_genome_colleto
    python tools/genome_scan.py --genus all
"""
import argparse
import io
import json
import os
import time
import urllib.request

from Bio import Entrez

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CACHE = os.path.join(REPO, ".cache", "genome_scan2")
os.makedirs(CACHE, exist_ok=True)
API = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/taxon/{}/dataset_report?page_size=1000"
SLUG = {"Colletotrichum": "colletotrichum", "Trichoderma": "trichoderma", "Diaporthe": "diaporthe"}


def log(m):
    print(m, flush=True)


def species_of(g):
    p = os.path.join(REPO, "data", "species_{}_all.txt".format(SLUG[g]))
    return [l.strip() for l in io.open(p, encoding="utf-8") if l.strip()]


def taxid_of(sp):
    cp = os.path.join(CACHE, "tax_{}.json".format(sp.replace(" ", "_")))
    if os.path.exists(cp):
        return json.load(io.open(cp, encoding="utf-8"))
    tid = None
    try:
        h = Entrez.esearch(db="taxonomy", term='"{}"'.format(sp), retmax=1)
        r = Entrez.read(h)
        h.close()
        if r.get("IdList"):
            tid = r["IdList"][0]
        time.sleep(0.36)
    except Exception as e:
        log("   ! taxonomy 失败 {}: {}".format(sp, str(e)[:50]))
    json.dump(tid, io.open(cp, "w", encoding="utf-8"))
    return tid


def reports_of(tid):
    cp = os.path.join(CACHE, "rep_{}.json".format(tid))
    if os.path.exists(cp):
        return json.load(io.open(cp, encoding="utf-8"))
    got = []
    token = None
    for _ in range(6):
        url = API.format(tid) + ("&page_token=" + token if token else "")
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                d = json.loads(r.read().decode("utf-8"))
        except Exception as e:
            log("   ! dataset_report 失败 taxid={}: {}".format(tid, str(e)[:60]))
            break
        got += d.get("reports", [])
        token = d.get("next_page_token")
        if not token:
            break
        time.sleep(0.3)
    json.dump(got, io.open(cp, "w", encoding="utf-8"))
    time.sleep(0.25)
    return got


def run(genus, out):
    sps = species_of(genus)
    log("== {} ｜ {} 个种 ==".format(genus, len(sps)))
    rows, per = [], {}
    for n, sp in enumerate(sps, 1):
        tid = taxid_of(sp)
        reps = reports_of(tid) if tid else []
        has_type, n_unv, seen = False, 0, set()
        for rep in reps:
            ai = rep.get("assembly_info", {}) or {}
            bs = ai.get("biosample", {}) or {}
            key = bs.get("accession") or ai.get("assembly_name") or ""
            if key in seen:          # GCA/GCF 成对出现，按 biosample 去重
                continue
            seen.add(key)
            attr = {a.get("name", ""): a.get("value", "") for a in (bs.get("attributes") or [])}
            tm = (attr.get("type-material") or attr.get("type_material") or "").strip()
            notes = ai.get("genome_notes") or []
            unv = any("unverified source organism" in str(x).lower() for x in notes)
            if tm:
                has_type = True
            if unv:
                n_unv += 1
            rows.append({
                "属": genus, "种": sp,
                "记录物种名": (rep.get("organism", {}) or {}).get("organism_name", ""),
                "assembly": ai.get("assembly_name", ""),
                "biosample": bs.get("accession", ""),
                "strain": attr.get("strain", ""),
                "culture_collection": attr.get("culture_collection", attr.get("culture-collection", "")),
                "type_material": tm[:200],
                "genome_notes": ";".join(str(x) for x in notes)[:120],
            })
        per[sp] = {"n": len(seen), "has_type": has_type, "unv": n_unv}
        if n % 25 == 0 or n == len(sps):
            log("   [{}/{}] {} ｜ 组装 {} ｜ 有模式声明 {}".format(
                n, len(sps), sp, len(seen), "是" if has_type else "否"))

    od = os.path.join(REPO, out)
    os.makedirs(od, exist_ok=True)
    cols = ["属", "种", "记录物种名", "assembly", "biosample", "strain",
            "culture_collection", "type_material", "genome_notes"]
    tsv = os.path.join(od, "{}_基因组层.tsv".format(SLUG[genus]))
    with io.open(tsv, "w", encoding="utf-8-sig", newline="") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")

    withg = [s for s, v in per.items() if v["n"] > 0]
    hast = [s for s in withg if per[s]["has_type"]]
    summ = {
        "属": genus, "种总数": len(sps), "有基因组的种": len(withg),
        "有模式菌株基因组的种": len(hast),
        "基因组全部来自非模式菌株的种": len(withg) - len(hast),
        "非模式占比%": round(100.0 * (len(withg) - len(hast)) / max(1, len(withg)), 1),
        "组装数": sum(v["n"] for v in per.values()),
        "官方判归属错误的组装": sum(v["unv"] for v in per.values()),
        "明细": tsv,
    }
    json.dump(summ, io.open(os.path.join(od, "{}_基因组层_汇总.json".format(SLUG[genus])),
                            "w", encoding="utf-8"), ensure_ascii=False, indent=2)
    log("\n== {} 汇总 ==".format(genus))
    for k, v in summ.items():
        if k != "明细":
            log("   {:<26} {}".format(k, v))
    return summ


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--genus", default="all")
    ap.add_argument("--out", default="out_genome")
    a = ap.parse_args()
    Entrez.email = os.environ.get("NCBI_EMAIL", "audit@example.org")
    tgt = list(SLUG) if a.genus.lower() == "all" else [a.genus]
    res = [run(g, a.out) for g in tgt]
    json.dump(res, io.open(os.path.join(REPO, a.out, "三属_基因组层_汇总.json"), "w",
                           encoding="utf-8"), ensure_ascii=False, indent=2)
    log("\n全部完成 → {}".format(os.path.join(REPO, a.out)))
