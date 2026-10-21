# -*- coding: utf-8 -*-
"""J 规则库：把《类型材料元数据审计规范 v1》的 J-01…J-21 编码为可执行检查。

规范的 J 表原本是散文，本模块把它变成代码。**每条规则都标注了自动化程度与依赖**：

  offline      纯本地计算，不需要任何外部源
  needs_table  需要一个外部来源表（用户提供 TSV/JSON）
  needs_net    需要联网或外部数据库适配器

设计立场（抄自规范零节）：**脚本产出候选，人产出结论**。
  本模块只做「并置 + 判据 + 分级」，不做终审。凡涉及「谁对」的判断，
  一律输出为候选，并在 finding 里带上证据与档位。

关键实现约束（规范原文里点名的坑，逐条落到代码）：

  ⚠️ J-01 的机械陷阱：**「对方没给值」与「对方给了别的值」在集合运算里同形**
     （空集与任何集合的交集都是空集），但含义相反。故必须单列分支——
     见 `j01_cross_source` 里 status ∈ {agree, disagree, only_a, only_b} 四态。
     把 only_* 混进 disagree 会把"不可检验"的行误当可疑（本项目曾误算 31 行）。

  ⚠️ J-15 是元规则且最易漏：`CBS 147945` **不在** NCBI 给 *C. falcatum* 的
     type_material 字段里，所以用「类型材料字段」枚举时它从未被检查过。
     见 `j15_enumeration_scope`——它比较**两份枚举**，而不是检查单个值。

  ⚠️ J-17 是口径修正，缺了它会系统性高估传播：光数「字符串出现在多少篇文献」
     会把「该号的正主在正常使用它」算成「错误在传播」。
     故 `j17_citation_with_attribution` 要求每个引用都带**归属判据**。

  ⚠️ J-18 的适用范围必须声明：只适用于「错号指向另一个属」的跨属子类；
     同属错误没有「另一方」，指标不适用。本项目 10 条错值中适用 2 条。
"""
import os
import re
from collections import Counter, defaultdict

# ---------------------------------------------------------------- 工具

def edit_distance_1(a, b):
    """判断两串是否编辑距离 ≤1（替换/插入/删除各一）。"""
    if a == b:
        return True
    la, lb = len(a), len(b)
    if abs(la - lb) > 1:
        return False
    if la == lb:
        return sum(1 for x, y in zip(a, b) if x != y) == 1
    if la > lb:
        a, b = b, a
    i = j = 0
    diff = 0
    while i < len(a) and j < len(b):
        if a[i] != b[j]:
            diff += 1
            j += 1
            if diff > 1:
                return False
        else:
            i += 1
            j += 1
    return True


def split_code(num):
    """'CBS 147945' -> ('CBS', '147945')；容忍冒号/下划线/无空格。"""
    m = re.match(r"^\s*([A-Za-z][A-Za-z\-\.]{0,11})\s*[:：_\-]?\s*(.+?)\s*$", num or "")
    return (m.group(1).upper(), m.group(2)) if m else ("", (num or "").strip())


def findings(rule, subject, detail, severity="candidate", evidence=""):
    return {"rule": rule, "subject": subject, "detail": detail, "severity": severity, "evidence": evidence}


# ---------------------------------------------------------------- J-01

def j01_cross_source(claims_a, claims_b):
    """跨来源对账（规范评为『最高产』：9 条错值中 7 条由它冒出）。

    claims_x: {物种: 号 或 [号...]}
    返回四态之一，**四态不可合并**：
      agree     两边给出同一个号
      disagree  两边都给了值且不同        <-- 候选错值
      only_a / only_b  只有一侧给了值      <-- 不可检验，**不是**可疑

    ⚠️ **输出约定**：只有 `disagree` 是候选（severity=candidate）；
       `agree` 记 ok、`only_*` 记 untestable——**后两者不得写进 findings 表**，
       否则会把"不可检验"当成"有问题"（本项目曾误算 31 行）。
    """
    out = []
    for sp in sorted(set(claims_a) | set(claims_b)):
        A = claims_a.get(sp) or []
        B = claims_b.get(sp) or []
        A = [A] if isinstance(A, str) else list(A)
        B = [B] if isinstance(B, str) else list(B)
        if A and B:
            if set(A) & set(B):
                out.append(findings("J-01", sp, "agree", "ok", "{} ∩ {}".format(A, B)))
            else:
                out.append(findings("J-01", sp, "disagree: A={} vs B={}".format(A, B),
                                    "candidate", "两侧都有值且无交集"))
        elif A and not B:
            out.append(findings("J-01", sp, "only_a", "untestable", "A={}".format(A)))
        elif B and not A:
            out.append(findings("J-01", sp, "only_b", "untestable", "B={}".format(B)))
    return out


# ---------------------------------------------------------------- J-07

def j07_one_number_two_species(number_to_species):
    """一号两主：**同一侧内部**，同一个号出现在两个物种下。

    number_to_species: {号: set(物种)}
    规范示例：`GMBC 0209` 同时指给两个同文发表的新种（本项目 E7）。
    """
    out = []
    for num, sps in sorted(number_to_species.items()):
        if len(sps) > 1:
            out.append(findings("J-07", num, "一号两主：{}".format(sorted(sps)),
                                "candidate", "同侧内部同号多主，至少一方有误或该号对应不唯一保藏物"))
    return out


# ---------------------------------------------------------------- J-08

def j08_batch_regularity(records, group_key=None, digit_field="num",
                         min_n=5, dominance=0.7):
    """批量规律异常：同批数据里位数/连号的孤例。

    records: [{'group':..., 'num':...}] 或 [{'num':...}]

    ⚠️ **假阳性驱动**（实测：12 篇炭疽菌论文上出 14 条候选，抽检几乎全是假阳性）：
      本条要求「同批」是**同构编号体系**。`LC 0951` 与 `LC 13850` 只是流水号，
      LC 本就不定宽 ⇒ 位数不同**不是异常**。真异常（E8：`CGMCC 3.2051` vs `3.20514`）
      之所以成立，是因为 CGMCC 的 `3.xxxxx` **固定 5 位**。
      ⇒ 故加了两个门槛：**同批样本数 ≥ min_n** 且 **众数位数占比 ≥ dominance**；
        不够门槛的批直接跳过，**不猜**。使用者应按机构补一张「是否定宽」表来提高精度。
    """
    out = []
    groups = defaultdict(list)
    for r in records:
        groups[r.get(group_key) if group_key else "_all"].append(r)
    for g, rs in groups.items():
        lens = Counter(len(str(r[digit_field])) for r in rs if r.get(digit_field))
        if len(lens) < 2:
            continue
        major, major_n = lens.most_common(1)[0]
        if len(rs) < min_n or major_n / len(rs) < dominance:
            continue                      # 非同构批，不做位数判断
        for r in rs:
            L = len(str(r.get(digit_field) or ""))
            if L and L != major:
                out.append(findings("J-08", "{}/{}".format(g, r[digit_field]),
                                    "位数孤例：{} 位（同批众数 {} 位，n={}）".format(L, major, len(rs)),
                                    "candidate", "批量规律异常，规范示例 E8"))
    return out


# ---------------------------------------------------------------- J-13

def j13_near_miss(numbers_a, numbers_b):
    """数字级近失配（规范评为『最高产』的数字级判据）。

    ✅ 判据：互为前缀／互为后缀／编辑距离 1。
    规范示例命中 4 个物种：E1、E2、E6、E10。
    ⚠️ 与 J-08 互补：若一侧压根没给该机构的号，J-13 抓不到，只能靠 J-08。
    """
    out = []
    for a in sorted(set(numbers_a)):
        ca, na = split_code(a)
        for b in sorted(set(numbers_b)):
            cb, nb = split_code(b)
            if not na or not nb or na == nb:
                continue
            if ca and cb and ca != cb:
                continue          # 不同机构前缀不比较
            if na.startswith(nb) or nb.startswith(na) or edit_distance_1(na, nb):
                out.append(findings("J-13", "{} vs {}".format(a, b),
                                    "近失配：前缀相同、数字互为前缀/后缀或编辑距离 1",
                                    "candidate", "规范示例 E1/E2/E6/E10"))
    return out


# ---------------------------------------------------------------- J-14

RE_T_MARK = re.compile(r"([A-Z][a-z\-]{3,})\s+([A-Z][A-Za-z]{1,10}\s*\d[\w\-]*)_T\b")

def j14_figure_vs_body(figure_text, body_text, genus):
    """同一文献：图标签 vs 正文模式指定段（同源内部矛盾）。

    规范实测：图标签是 **5** 条错值的所在源，而正文裁决了 **3** 条；
    **图的可靠性显著低于正文**（图是转引层的转引）。
    故判定优先级：**正文有 Typus 段时，图标签让位于正文**。
    """
    fig_nums = defaultdict(set)
    for m in RE_T_MARK.finditer(figure_text or ""):
        fig_nums[m.group(1).lower()].add(m.group(2).replace(" ", ""))
    body_nums = defaultdict(set)
    for m in re.finditer(r"ex-type culture\s+([^);]{5,200})", body_text or "", re.I):
        for code, num in re.findall(r"\b([A-Z][A-Za-z]{1,11})\s*[:：]?\s*(\d+[\w.\-]*)", m.group(1)):
            body_nums[None].add(code + num)          # 归属由调用方按物种块细化
    out = []
    for ep, nums in sorted(fig_nums.items()):
        if not body_nums.get(None):
            continue
        if not (nums & body_nums[None]):
            out.append(findings("J-14", "{}. {}".format(genus, ep),
                                "图标签 {} 与正文的 ex-type 号无交集".format(sorted(nums)),
                                "candidate", "图 vs 正文；按规范正文优先"))
    return out


# ---------------------------------------------------------------- J-15

def j15_enumeration_scope(enum_a, enum_b, label_a="枚举A", label_b="枚举B"):
    """⚙️ 元规则：**枚举范围本身决定检出**。

    规范实测：层一用「类型材料字段里的号」枚举（212 个）→ 检出 **0**；
    改用「全部记录标识符」枚举 → 立刻抓到 E1。
    `CBS 147945` 不在 type_material 字段里，所以前者**从来没检查过它**。
    """
    A, B = set(enum_a), set(enum_b)
    only_b = sorted(B - A)
    return {
        "label_a": label_a, "label_b": label_b,
        "n_a": len(A), "n_b": len(B),
        "only_in_b": only_b[:50], "n_only_in_b": len(only_b),
        "note": "⚠️ 若 only_in_b 非空，说明枚举 A 把一批对象排除在检查之外——"
                "不是检查无效，是枚举范围漏了。凡报『检出 0』必先跑本条。",
    }


# ---------------------------------------------------------------- J-17

def j17_citation_with_attribution(citations):
    """传播计数**必须带归属判据**（规范：缺了它会系统性高估传播）。

    citations: [{'ref':..., 'number':..., 'attributed_to':...,
                 'expected_organism':...}]
    ✅ 只有「该文献把这个号归给了**被审分类单元**」才算传播；
       文献把它当**它的正主**用 = 正确使用，**零传播**。
    规范示例：`CBS 147945` 4 篇全归给炭疽菌（真传播）；
             `CBS 132551` 1 篇归给 *Arthroderma*（**零传播**）。
    """
    propagated, correct_use = [], []
    for c in citations:
        if c.get("attributed_to") and c.get("expected_organism") \
           and c["attributed_to"] == c["expected_organism"]:
            propagated.append(c)
        else:
            correct_use.append(c)
    return {
        "n_total_hits": len(citations),
        "n_propagated": len(propagated),
        "n_correct_use_excluded": len(correct_use),
        "propagated_refs": [c.get("ref") for c in propagated],
        "excluded": [{"ref": c.get("ref"), "attributed_to": c.get("attributed_to")} for c in correct_use],
        "note": "光数字符串命中会把「正主在正常用它」算成传播；故必须带归属判据。",
    }


# ---------------------------------------------------------------- J-18

def j18_conspicuity(records_of_number, records_of_wrong_side):
    """显眼度 = 该号在 GenBank 里**被错误那一方**占据的比例（可计算）。

    ⚠️ **适用范围必须声明**：只适用于「错号指向**另一个属**」的跨属子类；
       同属错误（如 E2、E9）没有「另一方」，指标不适用。
       本项目 10 条错值中适用 = **2 条**。
    规范示例：E1 = 5/6 = 0.8333 → 传播 4 篇；E6 = 0/9 = 0.0000 → 传播 0 篇。
    """
    n = records_of_number
    if not n:
        return None
    return {"conspicuity": round(records_of_wrong_side / n, 4),
            "n_records": n, "n_wrong_side": records_of_wrong_side,
            "applicable": "仅跨属子类；同属错误无『另一方』，不适用"}


# ---------------------------------------------------------------- J-20

def j20_note_type_status(genbank_text):
    """记录的 `/note` 自带 `Type status:` —— 可机读的机构/项目级指认。

    规范示例：`MH858975` 的 `/note` 自带
      `Type status: Isotype of C. gloeosporioides f. alata …; Holotype of C. alatae …`
    ⚠️ FASTA 只给 `>header` + 碱基，**不带限定符** ⇒ 必须用 GenBank 格式取。
    """
    out = []
    for m in re.finditer(r'/note\s*=\s*"([^"]*Type status:[^"]*)"', genbank_text or "", re.I):
        body = m.group(1)
        for st in re.findall(r"(Isotype|Holotype|Epitype|Neotype|Paratype)\s+of\s+([^;.]+)", body):
            out.append(findings("J-20", st[1].strip(), "Type status: {} of {}".format(*st[:1] + (st[1].strip(),)),
                                "evidence", "记录内声明，可机读"))
    return out


# ---------------------------------------------------------------- 注册表



def coverage():
    """规则库覆盖度自检：哪些能离线跑、哪些缺源。"""
    c = Counter(r["auto"] for r in all_rules())
    return {"total": len(all_rules()), "by_auto": dict(c),
            "offline_ready": [r["id"] for r in all_rules() if r["auto"] == "offline"],
            "needs_source": [r["id"] for r in all_rules() if r["auto"] != "offline"]}

# ---------------------------------------------------------------- J-04 / J-05 / J-09（需外部源）

def j04_catalogue_status(code, number, name_claimed, source):
    """目录模式标注：号 → 目录里记的「某某的模式」。

    判据：目录 status 里的**学名**与论文声称的物种不符 ⇒ 被占型候选。
    规范示例：`CBS 127945` = "Neotype of *Colletotrichum falcatum*"；
             `CBS 147945` 无标注 ⇒ 错值。
    实测（本地 CBS 数据 62,761 行）：E1/E4/E5/E6/incanum 的关键号全部可查。
    """
    st = source.status_of(code, number)
    if st is None:
        return [findings("J-04", "{} {}".format(code, number),
                         "目录无此号记录", "untestable", "目录未收录，不能据此定罪")]
    if st.get("rank") is None:
        return [findings("J-04", "{} {}".format(code, number),
                         "目录有该号但**无模式标注**（status 空）", "candidate",
                         "论文当模式用，目录无标注" + ("；原文：" + st.get("raw", "") if st.get("raw") else ""))]
    if name_claimed and st["name"] and st["name"].split()[-1].lower() != name_claimed.split()[-1].lower():
        return [findings("J-04", "{} {}".format(code, number),
                         "目录记的是 {} of *{}*，而论文声称是 *{}*".format(st["rank"], st["name"], name_claimed),
                         "candidate", "被占型：该号在目录里属于另一个物种")]
    return [findings("J-04", "{} {}".format(code, number),
                     "目录记 {} of *{}*".format(st["rank"], st["name"]), "ok", "与声称一致")]


def j05_reverse_by_name(name, claimed_numbers, source):
    """按学名反查目录，得到该种的模式号；与声称的号比对。

    用途：抓**漏登**与**异名混入**。
    规范示例：`CBS 125378` 是 *C. hymenocallidis* 的模式，不是 *C. siamense* 的。
    实测：本地 CBS 数据里 `CBS 125378` status = "holotype of Colletotrichum hymenocallidis"。
    """
    got = source.type_numbers_of(name)
    claimed = {n.replace(" ", "").upper() for n in claimed_numbers}
    gotset = {n.replace(" ", "").upper() for n in got if n}
    out = []
    if gotset and not (claimed & gotset):
        out.append(findings("J-05", name,
                            "目录按学名反查到的模式号 {} 与声称的 {} 无交集".format(sorted(gotset), sorted(claimed)),
                            "candidate", "漏登或异名混入"))
    elif gotset:
        out.append(findings("J-05", name, "一致（目录 {}）".format(sorted(gotset)), "ok", ""))
    return out


def j09_code_collision(code, source, paper_usage=""):
    """代码撞名：机构代码在登记表里的**类型**与文献里的**用法**不符。

    规范示例：原文写 "living culture LC…"，登记表记 `LC` = herbarium。
    ⚠️ 该输入错了会导致整批号走错分支（规范拾节把它列为已知陷阱）。
    """
    t = source.institution_type(code)
    if t is None:
        return [findings("J-09", code, "登记表无此代码", "untestable", "")]
    hit = "culture" in t and re.search(r"herbarium|specimen|标本", paper_usage or "", re.I)
    hit = hit or ("herbarium" in t and re.search(r"living culture|culture|培养", paper_usage or "", re.I))
    return [findings("J-09", code, "登记表类型 = {}；文献用法 = {}".format(t, (paper_usage or "")[:40]),
                     "candidate" if hit else "ok", "代码撞名" if hit else "")]


# ==================== 注册表 ====================
# ⚠️ 必须放在**所有被引用的函数定义之后**——RULES 在导入时就求值，
#    若某个 impl 函数尚未定义会 NameError（本项目踩过一次）。


# ---------------------------------------------------------------- J-14（单元版，推荐）

RE_EXTYPE_NUM = re.compile(r"(?:ex-type\s+(?:living\s+)?culture|ex-type)\s*[:=,\s]*"
                           r"([A-Z]{2,8}\s*[:\-]?\s*\d[\w.\-]*)", re.I)


def _norm_num(x):
    """号的正规化：去空白、去分隔符、**去掉尾部标点**。

    ⚠️ 尾部标点是实测 bug：同一个号在句末会带句点，`'GMBC0209'` 与 `'GMBC0209.'`
       于是成了两个不同的键，导致「一号两主」判不出来（E7 就卡在这）。
       号码里不会以标点结尾，故一律剥掉。
    """
    return re.sub(r"\s+", "", str(x)).replace(":", "").rstrip(".,;、，。")


def _nums(text):
    return {_norm_num(m.group(1)) for m in RE_EXTYPE_NUM.finditer(text or "")}


def j14_units(raw, genus, species_list=None):
    """J-14（单元版）：**同源内部矛盾** —— 图注里的 ex-type 号 vs 正文 Typus 段的号。

    为什么必须用单元版（实测教训）：
      · 压平文本 + 全局正则会抓到**别的图**的号（实测抓到 BRIP 15842a，纯误报）。
      · 图注的物种归属**不能靠文档序**——图是浮动的，caption 常出现在种标题之前。
        ⇒ 必须按**图注内文是否指名该种**来归属（caption 一般会写 "Colletotrichum X (...)"）。

    实测例（SiM 2022，*C. zhaoqingense*）：
      图注 `Colletotrichum zhaoqingense (… ex-type culture NN058985 …)`
      正文 Typus `ex-type culture CGMCC 3.20536 = LC13878 = NN071035`
      ⇒ 图注独有 NN058985 ⇒ 矛盾（正文优先）。
    """
    from .textutil import strip_tags, species_units, units
    us = units(raw)
    caps = [t for k, t in us if k == "caption"]
    su = species_units(raw, genus)
    body = {k: [t for t in v if t not in caps] for k, v in su.items()}
    out = []
    for ep in (species_list or sorted(su)):
        fig = set()
        for c in caps:
            if re.search(r"\b" + re.escape(genus) + r"\s+" + re.escape(ep) + r"\b", c, re.I):
                # 只取该种**自己的括号内**的号（避免抓到同图里的其他种）
                for m2 in re.finditer(re.escape(ep) + r"\s*\(([^)]*)\)", c, re.I):
                    fig |= _nums(m2.group(1))
                if not fig:
                    fig |= _nums(c)
        bd = set()
        for t in body.get(ep, []):
            bd |= _nums(t)
        only_fig = sorted(fig - bd)
        if fig and bd and only_fig:
            out.append(findings("J-14", "{}. {}".format(genus, ep),
                                "图注独有 {}（正文 Typus 为 {}）".format(only_fig, sorted(bd)[:3]),
                                "candidate", "图是转引层的转引；判定以正文 Typus 段为准"))
        elif fig and not bd:
            out.append(findings("J-14", "{}. {}".format(genus, ep),
                                "图注有 {} 而正文 Typus 段无号".format(sorted(fig)), "candidate", ""))
    return out


def j07_units(raw, genus):
    """J-07（单元版）：一号两主 —— 同一个号在不同物种名下被声称为 ex-type。

    为什么单元版更准：压平文本切块会跨表格，使同一个号落进多个物种块，
    从而产生大量**假阳性**（实测 12 条候选几乎全假）。按单元归属可避免。
    """
    from .textutil import species_units
    su = species_units(raw, genus)
    owners = {}
    for ep, texts in su.items():
        for t in texts:
            for x in _nums(t):
                owners.setdefault(x, set()).add(ep)
    two = {k: sorted(v) for k, v in owners.items() if len(v) > 1}
    return [findings("J-07", k, "一号两主：{}".format(v), "candidate",
                     "同侧内部同号多主（单元归属）") for k, v in sorted(two.items())]

def all_rules():
    """返回全部规则条目。

    ⚠️ **为什么是函数而不是模块级常量**：本注册表在**导入时**求值，
        若某个 `impl` 函数在其后定义就会 NameError。本项目为此连踩三次
        （rules 的 RULES、sources 的 KINDS，各两次）。改成函数后，
        求值推迟到调用时，与定义顺序彻底解耦。
    """
    return [
 {"id": "J-01", "name": "跨来源对账",              "auto": "offline",     "impl": j01_cross_source,           "note": "最高产；四态不可合并"},
 {"id": "J-02", "name": "前缀 vs 登记代码表",        "auto": "needs_table", "impl": None,                       "note": "需 8910 个机构代码表"},
 {"id": "J-03", "name": "号的属分布",               "auto": "needs_net",   "impl": None,                       "note": "层一实测 212→1（0.47%）；只作缩范围后的复核"},
 {"id": "J-04", "name": "目录模式标注",             "auto": "needs_table", "impl": j04_catalogue_status,       "note": "✅ 本地 CBS 数据 62,761 行即可跑；实测覆盖 E1/E4/E5/E6/incanum"},
 {"id": "J-05", "name": "按学名反查目录",            "auto": "needs_table", "impl": j05_reverse_by_name,        "note": "✅ 本地 CBS 数据即可跑；抓漏登与异名混入"},
 {"id": "J-06", "name": "原文模式指定段",            "auto": "offline",     "impl": "designations",             "note": "PMC 结构化 XML 可机读；见 designations 模块"},
 {"id": "J-07", "name": "一号两主",                "auto": "offline",     "impl": j07_units,                  "note": "单元版（推荐）：按 XML 单元归属物种，避免块跨表格造成的假阳性"},
 {"id": "J-08", "name": "批量规律异常",             "auto": "offline",     "impl": j08_batch_regularity,       "note": "与 J-13 互补；⚠️ 假阳性高——LC 等不定宽编号的位数差异不是异常，需按机构补『是否定宽』表"},
 {"id": "J-09", "name": "代码撞名（用法）",          "auto": "needs_table", "impl": j09_code_collision,         "note": "需 Biocollections 登记表；输入错了整批走错分支"},
 {"id": "J-10", "name": "代码撞名（采集地）",         "auto": "needs_net",   "impl": None,                       "note": "提议-待验证"},
 {"id": "J-11", "name": "strain/collection 不交叉引用","auto": "needs_net",  "impl": None,                       "note": "检索层可区分性"},
 {"id": "J-12", "name": "序列注释错",               "auto": "needs_net",   "impl": None,                       "note": "本项目未用"},
 {"id": "J-13", "name": "数字级近失配",             "auto": "offline",     "impl": j13_near_miss,              "note": "最高产的数字级判据"},
 {"id": "J-14", "name": "图 vs 正文",              "auto": "offline",     "impl": j14_units,                  "note": "单元版（推荐）：图注按内文指名归属；图可靠性低于正文"},
 {"id": "J-15", "name": "⚙️ 枚举范围（元规则）",     "auto": "offline",     "impl": j15_enumeration_scope,      "note": "报『检出 0』前必跑"},
 {"id": "J-16", "name": "记录声明 × 目录",          "auto": "needs_net",   "impl": None,                       "note": "同属错误，不依赖框架"},
 {"id": "J-17", "name": "⚙️ 传播计数带归属",        "auto": "needs_table", "impl": j17_citation_with_attribution, "note": "缺了会系统性高估传播"},
 {"id": "J-18", "name": "显眼度",                  "auto": "needs_table", "impl": j18_conspicuity,            "note": "仅跨属子类适用（10 条中 2 条）"},
 {"id": "J-19", "name": "号 ↔ Taxonomy 反查",       "auto": "needs_net",   "impl": None,                       "note": "补位/裁定，不做筛查（实测 30→0）"},
 {"id": "J-20", "name": "/note 的 Type status",    "auto": "offline",     "impl": j20_note_type_status,       "note": "必须用 GenBank 格式，FASTA 不带限定符"},
 {"id": "J-21", "name": "可检索性诊断",             "auto": "needs_net",   "impl": None,                       "note": "逐字段试 Title/Strain/All"},
    ]
