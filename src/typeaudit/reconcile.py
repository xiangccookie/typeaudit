# -*- coding: utf-8 -*-
"""物种级对账：把「两侧的号」按判定规则 ①–⑥ 归到判定档。

判定规则（抄自《口径说明.md》第二节，**顺序即优先级**）：
    ① 人工裁定表命中                        → 🔴 有确认错值（最高优先，机械规则不得覆盖）
    ② 两侧有共有号                          → ✅ 一致
    ③ 一侧全为{标本馆, 未登记}、另一侧含培养物 → 🔶 标本/活培养物配对
    ④ 两侧均含培养物保藏代码                → 🔶 疑似同株双保藏
    ⑤ 两侧都有号但类型不明                  → ⚠️ 需人工
    ⑥ 仅一侧有号                            → 仅 NCBI 有 / 仅 SiM 有  （**≠ 不一致**）

⚠️ 两条铁律（都是踩过坑的）：
  · **rule ① 必须最高优先**，否则已人工确认的错值会被机械规则冲掉（实测发生过一次）。
  · **rule ⑥ 不可与 ⑤ 合并**：「空」与「不同」在集合运算里同形（空集与任何集合的交集都是空集），
    合并会把不可检验的行误当可疑（本项目曾误算 31 行）。
"""
import csv
import io
import os
import re

# 机构类型：培养物保藏 / 标本馆。优先查 Biocollections 适配器，这里是兜底表（未登记的常用代码）。
CULTURE_CODES = set("""CBS CGMCC NRRL ATCC MAFF IFO JCM NBRC ICMP IMI BRIP MUCL LC NN GMBC ACCC CNUCC
BCC CPC CBS-KNAW CUFC DAOM UAMH NRRL DSM""".split())
HERBARIUM_CODES = set("""BPI HMAS HMJAU K M NY C G S W PC P UPS ZT FH FLAS MICH MO MSC OSC PH PRM UC WIS
LE LPS MEL TUR CM E F L LD O UME CFMR KRAM HAL WA NIFS HAST GDGM HALL TENN VPI NYS HLA GACP
MFLU MFLUCC HASTT""".split())

# ⚠️ 分隔符要允许「连字符」：实测目录里同一号有 `CGMCC-3.18952` 与 `CGMCC 3.18952`
#    两种写法，只允许空格/冒号会整条解析失败（静默丢号）。
RE_NUM = re.compile(r"\b([A-Z][A-Za-z]{1,11})\s*[\s:#\-]?\s*(\d[\w.\-]*)")


def split_nums(field):
    """把一侧的号字符串拆成 [(code, num)]，去重保序。"""
    out, seen = [], set()
    for part in re.split(r"[|;]", field or ""):
        for m in RE_NUM.finditer(part):
            k = (m.group(1).upper(), re.sub(r"[\s:\-]", "", m.group(2)))
            if k not in seen:
                seen.add(k); out.append(k)
    return out


def kind_of(code, inst=None):
    """机构类型：culture / herbarium / unknown。inst 为可选的适配器（有 institution_type）。"""
    if inst is not None:
        try:
            t = inst.institution_type(code)
            if t in ("culture", "herbarium"):
                return t
        except Exception:
            pass
    c = (code or "").upper()
    if c in CULTURE_CODES: return "culture"
    if c in HERBARIUM_CODES: return "herbarium"
    return "unknown"


def verdict(left_field, right_field, species=None, manual_hits=(), inst=None):
    """按 ①–⑥ 给出判定。返回 (verdict, 说明)。"""
    L, R = split_nums(left_field), split_nums(right_field)
    # ① 人工裁定优先（机械规则不得覆盖）
    if species and species in manual_hits:
        return "🔴 有确认错值", "rule ① 人工裁定表命中（最高优先）"
    if not L and not R:
        return "⚠️ 需人工", "两侧都无号"
    if not L:
        return "仅 SiM 无号 / 仅 NCBI 有", "rule ⑥ 仅一侧有号（≠ 不一致）"
    if not R:
        return "仅 NCBI 无号 / 仅 SiM 有", "rule ⑥ 仅一侧有号（≠ 不一致）"
    # ② 共有号（规范化后比较）
    ls = {c + n for c, n in L}
    rs = {c + n for c, n in R}
    if ls & rs:
        return "✅ 一致", "rule ② 两侧有共有号"
    # 类型分档
    lk = {kind_of(c, inst) for c, _ in L}
    rk = {kind_of(c, inst) for c, _ in R}
    l_cult, r_cult = "culture" in lk, "culture" in rk
    l_all_spec = lk <= {"herbarium", "unknown"}
    # ③ 一侧全为{标本馆,未登记}、另一侧含培养物
    if (l_all_spec and r_cult) or (rk <= {"herbarium", "unknown"} and l_cult):
        return "🔶 标本/活培养物配对", "rule ③ 一侧全为{标本馆,未登记}，另一侧含培养物"
    # ④ 两侧均含培养物保藏代码
    if l_cult and r_cult:
        return "🔶 疑似同株双保藏", "rule ④ 两侧均含培养物保藏代码（同株双保藏本来就该有两套号）"
    # ⑤ 都有号但类型不明
    return "⚠️ 需人工", "rule ⑤ 两侧都有号但机构类型不明"

# ==================== 实测：本模块的判别力与**固有盲区** ====================
# 全量重跑炭疽菌 293 行，与人工裁定表对照（284 行可比）：
#   交叉表**完美对角**——240 一致 / 31 仅一侧 / 5 ✅双保藏 / 5 ✅标本配对 / 2 🔶疑似同株 / 1 🔶标本配对
#   全 284 行逐行一致，零假阳性。
#
# 但 9 行真错值里只标出 7 行（**召回 78%**）。漏掉的 2 行是：
#   · C. danxiashanense：SiM=NN055218，NCBI 侧**同时**含 NN055218（共有）与错的 CGMCC 3.2051
#   · C. siamense     ：SiM=ICMP 18578，NCBI 侧**同时**含 ICMP 18578（共有）与混入的 CBS 125378
#
# ⚠️ **rule ② 的固有盲区：只看「交集非空」，不检查「有无多余/错误的值」。**
#    上面两条都是**"多发了一个错号"**（多值型），而不是"把唯一的号给错了"（替换型）。
#    只要两侧存在一对共有号，rule ② 就判一致，多出来的错号不被审视。
#
# ⇒ 分工（这才是正确用法）：
#    **对账层（本模块）负责缩范围**：293 → 53 行需人工（18.1%），命中 7/9 真错值；
#    **J 规则层负责定罪**：漏掉的 2 行由 J-05（按学名反查目录，抓异名混入/漏登）
#      与 J-13（数字级近失配，3.2051 vs 3.20514）捕获——实测两者都能抓到。
#
# 若要让对账层自己覆盖多值型，需给 rule ② 补一条子检查：
#    「两侧有共有号时，对**非共有号**逐个做归属核查（J-05 反查 / J-13 近失配）」。
#    ⚠️ 别把它做成"有非共有号就判可疑"——那会把所有双保藏/标本配对全部误报。

# ==================== rule ② 的子检查：非共有号的「归属核查」 ====================
from collections import Counter as _Counter


def sub_attribution(species, side_sim, side_ncbi, source=None, batch=None, synonyms=()):
    """rule ② 命中（两侧有共有号）之后，逐个审「一侧独有的号」是否有归属问题。

    为什么需要：rule ② 只看「交集非空」，**不检查有没有多余/错误的号**。
    实测漏掉的两条正错值都是"多发了一个错号"型：
        C. danxiashanense：左 NN055218 ↔ 右 NN055218（共有）+ 错的 CGMCC 3.2051
        C. siamense     ：左 ICMP 18578 ↔ 右 ICMP 18578（共有）+ 混入的 CBS 125378（异名模式号）

    ⚠️ **不可做成「有非共有号就判可疑」**——那会把所有双保藏/标本配对（本项目 13 行）全部误报。
       本函数要求**正面证据**，只有两条路线之一命中才报：
         路线 A（目录归属）：该号在保藏目录里登记的是**别的种**（source.status_of）
         路线 B（同批位数孤例）：该号在**同机构号批**里位数偏离众数（source 可缺）

    返回 [ {number, route, detail, severity} ]
    """
    def toks(field):
        return split_nums(field)
    L, R = toks(side_sim), toks(side_ncbi)
    lset = {c + n for c, n in L}
    rset = {c + n for c, n in R}
    shared = lset & rset
    if not shared:
        return []                       # 只在 rule ② 命中的场景下做
    epi = (species or "").split()[-1].lower()
    ok_names = {epi} | {str(x).lower() for x in synonyms}
    out = []
    for side, nums, own in (("SiM", L, lset), ("NCBI", R, rset)):
        for code, num in nums:
            if code + num in shared:
                continue
            # 路线 A：目录归属 —— 该号登记的是别的种
            if source is not None:
                try:
                    st = source.status_of(code, num)
                except Exception:
                    st = None
                nm = (st or {}).get("name") or ""
                if nm and nm.split()[-1].lower() not in ok_names:
                    out.append({"number": "{} {}".format(code, num), "route": "目录归属",
                                "detail": "目录记「{} of {}」，非本种".format(st.get("rank") or "?", nm),
                                "severity": "candidate", "side": side})
                    continue
            # 路线 B：**同子系列断号**（不是"同机构位数孤例"）
            # ⚠️ 判据演进（实测教训）：最初用「该机构全局位数众数」，结果 5 条假阳性。
            #    根因：`CGMCC 3.2051` 与 `CGMCC 3.9617` **都是 6 位**，光看位数分不开——
            #      3.2051 与 3.20509/3.20510/…/3.20536 同属 `3.205` 系列（28 个兄弟全是 7 位）⇒ 断号
            #      3.9617 是 CGMCC 旧编号，前缀 `3.961` **没有任何兄弟** ⇒ 孤立老号，不判
            #    ⇒ 正确判据是「**同前缀（前 len-1 位）的兄弟号**里位数不一致」，不是机构全局。
            if batch and len(num) >= 3:
                pref = num[:len(num) - 1]
                sibs = [n for c, n in batch if c == code and n.startswith(pref) and n != num]
                if len(sibs) >= 2:
                    mode, mn = _Counter(len(s) for s in sibs).most_common(1)[0]
                    if len(num) != mode and mn / len(sibs) >= 0.8:
                        out.append({"number": "{} {}".format(code, num), "route": "同子系列断号",
                                    "detail": "该号 {} 位，而同前缀 {} 的 {} 个兄弟号均为 {} 位（疑为断号）".format(
                                        len(num), pref + "…", len(sibs), mode),
                                    "severity": "candidate", "side": side})
    return out


def verdict_deep(side_sim, side_ncbi, species=None, manual_hits=(), inst=None,
                 source=None, batch=None, synonyms=()):
    """`verdict()` + rule ② 子检查。返回 (verdict, 说明, 子检查明细)。"""
    v, why = verdict(side_sim, side_ncbi, species=species, manual_hits=manual_hits, inst=inst)
    subs = sub_attribution(species, side_sim, side_ncbi, source=source, batch=batch, synonyms=synonyms)
    if v == "✅ 一致" and subs:
        first = subs[0]
        return ("🔶 疑似错值（非共有号归属异常）",
                "rule ② 命中，但子检查发现：{}（{}）".format(first["number"], first["detail"]),
                subs)
    return v, why, subs

# ==================== 实测：rule ② 子检查的召回/精确率权衡 ====================
# 在 293 行全量实测（基础规则 + 子检查）：
#
#   配置                       召回       假阳性(/284 干净行)
#   基础规则(③④⑥)              7/9 = 78%   0
#   + 子检查（本函数）           9/9 = 100%  5
#     ├ 只路线 A（目录归属）       +1         1
#     └ 只路线 B（同批位数）       +2         4
#
# 两条路线的假阳性根因**都不是算法问题，而是缺一张知识表**：
#
#   路线 B 的 4 条假阳性：**「同批等宽」假设在这些机构不成立**——
#     CGMCC 有两代编号（`3.9617` 旧 4 位 / `3.20514` 新 5 位，同一批共存）；
#     ICMP / BRIP / PDD 的号位数也随年代变化。
#     ⇒ 需要一张「机构编号是否定宽」表。
#
#   路线 A 的 1 条假阳性：*C. godetiae* 行的 `CBS 130251` 在目录里记
#     「Holotype of *Colletotrichum clavatum*」——而 `clavatum` 与 `godetiae`
#     各有自己的模式号（目录里 `godetiae → CBS 133.44`）。
#     **NCBI 常把异名的类型材料列在接纳名下**，故这条很可能是合法的。
#     ⇒ 需要一张「异名关系」表（规范 L-08 已点名）。
#
# 结论：**子检查要真正好用，缺的不是算法，是两张知识表**。
#   在此之前，建议只开路线 A（+1 召回、+1 假阳性），把路线 B 的结果标为「低置信候选」。

# ==================== L 表：合法解释（放行，但**留痕**） ====================
def synonym_explains(species, code, num, source, syn_table):
    """该「非共有号」是否可用**异名**解释掉（规范 L 表：合法解释，放行）。

    ⚠️ 为什么只是「降级留痕」而不是「清除」：
       MycoBank 的 `Current name` **不可当真值**——实测（Colletotrichum 1,252 行）：
         · 393 条被并名里，**131 条（33%）种加词完全相同**，那只是正字/性别变体
           （`C. acidae → C. acidi`、`C. aecidiicola → C. aecidiicolum`），**不是异名**；
         · `Colletotrichum godetiae` 在表里出现 **3 次**（三个不同 MyCoBank 号，
           其中一个 Name status 为 **Invalid**）。
       ⇒ 本函数只在「目录说该号属于 X 种」且「X 是 species 的异名」时给出解释，
         并**原样返回 MycoBank 号**供人工复核；命中者降级为 `🔷 合法（异名）`，仍留在表里。
    """
    if source is None or not syn_table:
        return None
    st = source.status_of(code, num)
    nm = (st or {}).get("name")
    if not nm:
        return None
    epi = (species or "").split()[-1].lower()
    # 该号在目录里登记的种名，是否在本种的异名表里
    for p in syn_table:
        if p["current"].split()[-1].lower() == epi and p["taxon"].split()[-1].lower() == nm.split()[-1].lower():
            return {"number": "{} {}".format(code, num), "route": "合法异名（L 表）",
                    "detail": "目录记它是 `{}` 的模式，而 `{}` 在 MycoBank 里现名为 `{}`（#{}）".format(
                        nm, nm, p["current"], p["mb"]),
                    "mb": p["mb"], "severity": "legitimate"}
    return None


def load_synonyms(path):
    """读异名表 TSV（taxon/current/mb/status）。"""
    out = []
    if path and os.path.exists(path):
        with io.open(path, encoding="utf-8-sig", newline="") as f:
            out = list(csv.DictReader(f, delimiter="\t"))
    return out
