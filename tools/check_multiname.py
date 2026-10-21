#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# --- console utf8 bootstrap（Windows 控制台默认 GBK，输出 ↔ ⚠️ 🔴 会崩；详见 src/typeaudit/_console.py）
import sys as _sys
for _s in (_sys.stdout, _sys.stderr):
    try:
        _s.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass
"""一株多名检测 + 三档判定（grep 友好文字标记）

用法
----
  python check_multiname.py --gb <genbank.txt>  --out <out.tsv> [--mycobank <mycobank_full.tsv>]
  python check_multiname.py --json <raw.json>  --out <out.tsv> [--mycobank <mycobank_full.tsv>]

输入
----
  --gb    GenBank flat file（`efetch -rettype gb` 的输出，可含多条记录）
  --json  已解析记录（{accession: {defn, cc, strain, host}}），用于免网络复跑
  --mycobank  MycoBank 全库 TSV（列含 Taxon name / Current name / Name status / Rank）。

输出
----
  TSV，每行一个「同一保藏号下物种名不唯一」的菌株，含三档判定。

三档判定（grep 文字，勿改成 emoji——项目里 emoji 已被其它语义占用）
----
  多名-3 号重复使用（真错）        ← 优先判：同号下 /strain 或其它保藏号不同 ⇒ 实为两株
  多名-1 异名（非错）              ← 各名字的 Current name 相同 ⇒ 只是换了名字
  多名-2 两皆接受名（需查出处）    ← 各名字都是接受名且不同 ⇒ 分类学修订/位阶变化，异名表判不了

判据来源
--------
  分组键一律用 /culture_collection（保藏号），**绝不用物种名**——物种名会随分类学修订漂移。
  保藏号写法归并：去空格、大写。优先 ICMP，其次 CBS/IMI/DAR/BRIP，最后 /strain。
"""
import argparse, csv, io, json, os, re, sys

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

# ---- 1. 输入解析 ----
def parse_gb(text):
    """把 GenBank flat file 解析成 {version: {defn, cc, strain, host}}。"""
    out = {}
    for b in re.split(r"(?m)^LOCUS\s+", text)[1:]:        # 注意 \s+：吃掉对齐空格
        m = re.match(r"(\S+)\s+(\d+)\s+bp", b)
        if not m:
            continue
        base = m.group(1)
        vm = re.search(r"(?m)^VERSION\s+(\S+)", b)
        ver = vm.group(1) if vm else base
        defn = (re.search(r"DEFINITION\s+(.+?)(?=\n\s*ACCESSION)", b, re.S) or [None, ""])[1]
        get = lambda k: re.findall(r'/' + k + r'="((?:[^"]|"[^"]*")*?)"', b)
        out[ver] = dict(defn=re.sub(r"\s+", " ", defn).strip(),
                        cc=get("culture_collection"), strain=get("strain"),
                        host=get("host"), tm=get("type_material"))
    return out

def norm_code(s):
    """保藏号/菌株号归一：去空格、冒号、点号、连字符，统一大写。
    项目实测：CBS 124.22 / CBS124.22 / CBS:124.22 是同一号的不同写法。"""
    s = re.sub(r"<[^>]*>", "", s or "")          # 剥掉国家标签：MFLUCC<THA>090551 → MFLUCC090551
    return re.sub(r"[\s:：.\-]", "", s).upper()

def strain_key(rec):
    """分组键：保藏号归一化（优先 ICMP）。"""
    pref = ["ICMP", "CBS", "IMI", "DAR", "BRIP", "CGMCC", "MFLUCC", "MAFF", "NBRC"]
    for p in pref:
        for c in rec["cc"]:
            if norm_code(c).startswith(p):
                return norm_code(c)
    return norm_code(rec["strain"][0]) if rec["strain"] else "?"

def other_cc(rec, key):
    """同记录里除分组键以外的其它保藏号（用于判「号重复使用」）。"""
    out = []
    for c in rec["cc"]:
        n = norm_code(c)
        if n != key and re.match(r"^[A-Z]{2,8}", n):
            out.append(n)
    return out

# ---- 2. 物种名提取（剥掉 f.sp. / f. / subsp. 等位阶后缀）----
RE_NAME = re.compile(r"\b((?:Colletotrichum|Glomerella|Gloeosporium)\s+[a-z][a-z\-]+)")

def names_of(defn):
    return sorted(set(RE_NAME.findall(defn or "")))

# ---- 3. MycoBank 索引（按需惰性构建）----
class MycoBank:
    def __init__(self, path):
        self.idx = {}
        if not path or not os.path.exists(path):
            self.ok = False
            return
        self.ok = True
        with io.open(path, encoding="utf-8-sig", errors="replace", newline="") as f:
            rd = csv.reader(f, delimiter="\t")
            hdr = next(rd)
            i_cur = hdr.index("Current name") if "Current name" in hdr else 1
            i_st = hdr.index("Name status") if "Name status" in hdr else None
            i_rk = hdr.index("Rank") if "Rank" in hdr else None
            for row in rd:
                if not row:
                    continue
                key = row[0].strip().lower()
                if key in self.idx:
                    continue
                self.idx[key] = (row[i_cur], row[i_st] if i_st is not None else "",
                                 row[i_rk] if i_rk is not None else "")

    def current(self, name):
        """返回 (Current name 归一名, status, rank)；查不到返回 (None, '', '')。"""
        if not self.ok:
            return None, "", ""
        v = self.idx.get(name.strip().lower())
        if not v:
            return None, "", ""
        cur, st, rk = v
        # 归一到「属 + 种」（去掉命名人、发表处）
        m = re.match(r"\s*([A-Z][a-z]+)\s+([a-z][a-z\-]+)", cur)
        return (m.group(1) + " " + m.group(2) if m else cur), st, rk

# ---- 4. 主流程 ----
def main():
    ap = argparse.ArgumentParser()
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--gb")
    src.add_argument("--json")
    ap.add_argument("--out", required=True)
    ap.add_argument("--mycobank", default=None)
    a = ap.parse_args()

    if a.gb:
        recs = parse_gb(io.open(a.gb, encoding="utf-8-sig", errors="replace").read())
    else:
        recs = json.load(io.open(a.json, encoding="utf-8"))
    print("读入记录：{} 条".format(len(recs)))

    mb = MycoBank(a.mycobank)
    print("MycoBank 索引：{}".format("{} 条".format(len(mb.idx)) if mb.ok else "未提供（将只做多名检测）"))

    # 按保藏号归组
    grp = {}
    for ver, r in recs.items():
        r.setdefault("defn", ""); r.setdefault("cc", []); r.setdefault("strain", []); r.setdefault("host", [])
        grp.setdefault(strain_key(r), []).append((ver, r))

    rows = []
    for key, lst in sorted(grp.items()):
        names = sorted({n for _, r in lst for n in names_of(r["defn"])})
        if len(names) < 2:
            continue
        # —— 判据用归一化值，显示用原始值（两套必须分开）——
        def _raw_codes(rs):
            return sorted({c.strip() for _, r in rs for c in (r["cc"] or [])})
        raw_cc = _raw_codes(lst)
        norm_cc = sorted({norm_code(c) for c in raw_cc})
        raw_strain = sorted({" ".join(r["strain"]).strip() for _, r in lst if r["strain"]})
        norm_strain = sorted({norm_code(" ".join(r["strain"])) for _, r in lst if r["strain"]})
        hosts = sorted({" ".join(r["host"]).strip() for _, r in lst if r["host"]})

        # 判据：分组键之外，还有没有【归一化后仍不同】的保藏号 / strain
        oth_norm = sorted({x for x in norm_cc if x != norm_code(key)})
        n_strain_distinct = len(norm_strain)

        verdict, why = "多名-2 两皆接受名（需查该菌株原始出处）", ""
        if len(oth_norm) > 1 or n_strain_distinct > 1:
            verdict = "多名-3 号重复使用（真错）"
            why = "同号下归一化后仍不唯一：strain_norm={} cc_norm={}".format(
                "|".join(norm_strain)[:60], "|".join(oth_norm)[:60])
        elif mb.ok:
            cur = {}
            for n in names:
                c, st, rk = mb.current(n)
                cur[n] = c
            uniq = {v for v in cur.values() if v}
            unknown = [k for k, v in cur.items() if not v]
            if len(uniq) == 1 and len(cur) > 1:
                verdict = "多名-1 异名（非错）"
                why = "各名 Current name 相同：{}".format(list(uniq)[0])
            else:
                why = "Current name 各异：{}".format(
                    "; ".join("{}→{}".format(k.split()[-1], v or "?") for k, v in cur.items())[:90])
            if unknown:
                why += "  ｜未收录(非有效名/未定名)：{}".format(",".join(x.split()[-1] for x in unknown))

        # 第一列显示【分组键的原始写法】（不是排序后的第一个 cc）
        disp = key
        for _, _r in lst:
            for _c in (_r["cc"] or []):
                if norm_code(_c) == norm_code(key):
                    disp = _c.strip(); break
            if disp != key:
                break
        rows.append([disp, str(len(lst)), str(len(names)),
                     ";".join(n.split()[-1] for n in names[:4]),
                     "|".join(raw_strain)[:40], "|".join([c for c in raw_cc if norm_code(c) != norm_code(key)])[:40],
                     "|".join(hosts)[:40], verdict, why])

    hdr = ["保藏号", "记录数", "物种名数", "物种名", "strain", "其它保藏号", "寄主", "判定", "依据"]
    with io.open(a.out, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(hdr); w.writerows(rows)
    print("\n检测到「一株多名」：{} 个保藏号".format(len(rows)))
    from collections import Counter
    for k, v in Counter(r[7] for r in rows).most_common():
        print("   {:<40}{}".format(k, v))
    print("\n清单 → {}".format(a.out))

if __name__ == "__main__":
    main()
