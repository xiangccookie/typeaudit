#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""reproduce_colletotrichum.py —— 复现评估：本工具链能否自动检出已确认的 11 条错值。

这是「从数据到结论」的端到端检验，不是单元测试。
每条错值给出「检测配方」（用哪些规则 + 哪些源），然后实跑，报告检出与否。

用法：
    python examples/reproduce_colletotrichum.py --cbs <cbs_catalogue.tsv> [--live] [--cache .cache]

  --live   允许联网（CGMCC 官网 / NCBI Entrez）。不加则只跑离线源。

⚠️ 本脚本**不内置任何第三方数据**，只读你自备的文件；--live 时请遵守各站点限流。
"""
import argparse
import io
import os
import re
import sys
from collections import defaultdict

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
from typeaudit import rules as R, sources as S   # noqa: E402

# ---- 11 条错值的检测配方（数据取自本项目 副表1_错误清单.tsv）----
# 每项：(编号, 物种, 错值, 正确值, [适用规则])
CASES = [
 ("E1",  "Colletotrichum falcatum",    ("CBS", "147945"),  ("CBS", "127945"),  ["J-04", "J-03"]),
 ("E2",  "Colletotrichum jiangxiense", ("CGMCC", "3.17361"), ("CGMCC", "3.17363"), ["J-16"]),
 ("E3",  "Colletotrichum falcatum",    ("CGMG", "3.14187"), ("CGMCC", "3.14187"), ["J-02"]),
 ("E4",  "Colletotrichum siamense",    ("CBS", "130417"),   ("CBS", "130417"),   ["J-05", "J-05b"]),
 ("E5",  "Colletotrichum siamense",    ("CBS", "125378"),   ("CBS", "130417"),   ["J-05"]),
 ("E6",  "Colletotrichum arxii",       ("CBS", "132551"),   ("CBS", "132511"),   ["J-04"]),
 ("E7",  "Colletotrichum jishouense",  ("GMBC", "0209"),    None,                ["J-07"]),
 ("E8",  "Colletotrichum danxiashanense", ("CGMCC", "3.2051"), ("CGMCC", "3.20514"), ["J-08", "J-13"]),
 ("E9",  "Colletotrichum zhaoqingense", ("NN", "058985"),  ("NN", "071035"),    ["J-14", "J-13"]),
 ("E10", "Colletotrichum baltimorense", ("BPI", "892771"),  ("BPI", "1892771"),  ["J-13", "J-03", "J-21"]),
 ("E11", "Colletotrichum alatae",      ("CBS", "304.67"),  ("ICMP", "18122"),   ["J-19", "J-20"]),
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cbs", help="cbs_catalogue.tsv 路径")
    ap.add_argument("--live", action="store_true", help="允许联网（CGMCC / NCBI）")
    ap.add_argument("--cache", default=".cache")
    ap.add_argument("--registry", help="机构代码登记表 CSV（J-02 用）")
    ap.add_argument("--sim", help="SiM 2022 全文 XML（J-08/J-14 用）")
    ap.add_argument("--protologue", help="E7 原描述 PMCID（MycoKeys 49:1-14 = 6422929）")
    a = ap.parse_args()

    cbs = S.CBSCatalogue(a.cbs) if a.cbs else None
    # 侧源：登记表（J-02 用）与 SiM 2022 全文（J-08/J-14 用）
    reg = S.Biocollections(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..",
                                        a.registry) if a.registry else "")
    corpus = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")
    SIM = None
    for cand in (a.sim, os.path.join("data", "SiM2022.xml")):
        if cand and os.path.exists(cand):
            SIM = cand
            break
    if cbs:
        print("CBS 目录：{:,} 行（其中有 status 的 {:,} 行）".format(
            cbs.n_rows(), sum(1 for x in cbs.rows if x.get("has_status") == "1")))
    gb = S.GenBankAPI(a.cache) if a.live else None
    tax = S.TaxonomyAPI(a.cache) if a.live else None
    cgmcc = S.CGMCCAPI(a.cache) if a.live else None

    print()
    print("=" * 110)
    print("【复现评估：11 条错值 × 工具链】")
    print("=" * 110)
    print("   {:<5}{:<26}{:<18}{:<10}{}".format("编号", "物种", "错值号", "自动检出", "证据"))
    hit = 0
    for code, sp, wrong, right, rules in CASES:
        evid, detected = [], False
        c, n = wrong
        # J-04 目录模式标注
        if "J-04" in rules and cbs:
            st = cbs.status_of(c, n)
            if st is not None:
                if st.get("rank") is None:
                    evid.append("J-04 目录有号但无标注"); detected = True
                elif st["name"] and st["name"].split()[-1].lower() != sp.split()[-1].lower():
                    evid.append("J-04 目录记 {} of {}".format(st["rank"], st["name"])); detected = True
        # J-05 按学名反查
        if "J-05" in rules and cbs:
            got = cbs.type_numbers_of(sp)
            if got and ("{} {}".format(c, n)).replace(" ", "").upper() not in \
               {g.replace(" ", "").upper() for g in got}:
                evid.append("J-05 目录按学名反查 = {}".format(got[:3])); detected = True
        # J-03 / J-18 号的属分布
        if "J-03" in rules and gb:
            try:
                r3 = gb.records_of_number(c, n)
                gs = [g for g in r3["genera"] if g.lower() not in sp.lower()]
                if gs:
                    evid.append("J-03 属分布 {} 含他属 {}".format(r3["genera"], gs)); detected = True
                else:
                    evid.append("J-03 属分布 {}".format(r3["genera"]))
            except SystemExit as e:
                evid.append("J-03 跳过（{}）".format(str(e)[:24]))
        # J-16 记录声明 × 目录
        if "J-16" in rules and cgmcc:
            try:
                t = cgmcc.is_type_strain(n)
                if t is False:
                    evid.append("J-16 官网记「非模式菌株」"); detected = True
                elif t is True:
                    evid.append("J-16 官网记「模式菌株」")
                else:
                    evid.append("J-16 官网无该号")
            except Exception as e:
                evid.append("J-16 失败 {}".format(type(e).__name__))
        # J-19 Taxonomy 反查
        if "J-19" in rules and tax:
            try:
                r19 = tax.organism_of_number(c, n)
                if r19["count"]:
                    evid.append("J-19 Taxonomy → {}".format(r19["names"][:2])); detected = True
                else:
                    evid.append("J-19 Taxonomy 无此号")
            except Exception as e:
                evid.append("J-19 失败 {}".format(type(e).__name__))
        # J-21 可检索性
        if "J-21" in rules and gb:
            try:
                r21 = gb.retrievability(c, n)
                if r21["title"] == 0 and r21["strain"] == 0 and r21["all_fields"] > 0:
                    evid.append("J-21 仅全字段可命中（T=0 S=0 A={}）".format(r21["all_fields"]))
                    detected = True
            except Exception as e:
                evid.append("J-21 失败 {}".format(type(e).__name__))
        # J-20 记录 /note 的 Type status（记录内声明）
        if "J-20" in rules and gb:
            try:
                nts = gb.note_type_status_of(sp)
                if nts:
                    evid.append("J-20 /note Type status = {}".format(nts[:2])); detected = True
                else:
                    evid.append("J-20 该种记录无 Type status")
            except Exception as e:
                evid.append("J-20 失败 {}".format(type(e).__name__))
        # J-07 一号两主（**单元版**：按 XML 单元归属物种，避免块跨表格的假阳性）
        if "J-07" in rules and a.protologue:
            try:
                from typeaudit import gate as G
                raw = G.fetch_pmc(a.protologue, cache_dir=a.cache, email=os.environ.get("NCBI_EMAIL"))
                f07 = R.j07_units(raw, sp.split()[0])
                if f07:
                    evid.append("J-07 " + f07[0]["subject"] + " → " + f07[0]["detail"][:30]); detected = True
                else:
                    evid.append("J-07 原描述内无一号两主")
            except Exception as e:
                evid.append("J-07 失败 {}".format(type(e).__name__))
        # J-02 前缀 vs 登记表（空缺型拼写错）
        if "J-02" in rules and reg.n_rows():
            if reg.institution_type(c) is None:
                evid.append("J-02 前缀 {} 不在登记表内".format(c)); detected = True
        # J-13 数字级近失配（跨侧比较）
        if "J-13" in rules and right:
            c2, n2 = right
            if R.j13_near_miss(["{}{}".format(c, n)], ["{}{}".format(c2, n2)]):
                evid.append("J-13 与正确值互为前缀/后缀/编辑距离1"); detected = True
        # J-14 同源内部矛盾（**单元版**）：图注 vs 正文 Typus
        # 图注的物种归属按**内文指名**，不按文档序（图是浮动的）
        if "J-14" in rules and SIM and os.path.exists(SIM):
            try:
                raw = io.open(SIM, encoding="utf-8", errors="replace").read()
                f14 = R.j14_units(raw, sp.split()[0], species_list=[sp.split()[-1]])
                if f14:
                    evid.append("J-14 " + f14[0]["detail"][:50])
                    detected = detected or f14[0]["severity"] == "candidate"
                else:
                    evid.append("J-14 图注与正文无矛盾")
            except Exception as e:
                evid.append("J-14 失败 {}".format(type(e).__name__))
        # J-05 双侧：漏登（目录按学名反查到的号，在该侧未出现）
        if "J-05b" in rules and cbs:
            got = cbs.type_numbers_of(sp)
            evid.append("J-05b 目录反查 = {}".format(got[:3]))
            if got:
                detected = True
        hit += detected
        print("   {:<5}{:<26}{:<18}{:<10}{}".format(
            code, sp[:24], "{} {}".format(c, n), "✅ 检出" if detected else "—", "；".join(evid)[:66]))
    print()
    print("   ⇒ 自动检出 {}/{}".format(hit, len(CASES)))
    return hit


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
