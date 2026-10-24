#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""属级类型材料流水线：给一个属名 → 自动跑一遍 → 出一张三列表

产出（主表，三列，单元格内嵌问题标注）
--------------------------------------
    种名 | 菌株号（含问题） | 模式菌（含问题）

  第 1 列 种名
      该属的物种（来自：① 框架论文；② NCBI Taxonomy 检索；③ 本地物种清单）

  第 2 列 菌株号（含问题）
      该种在【序列库记录】里声明的保藏号 —— 取 /culture_collection 与 /strain。
      单元格内以 「号 [标记]」 形式列出，标记说明该号的问题（见下"问题标记"）。

  第 3 列 模式菌（含问题）
      该种在【类型材料字段】里声明的号 —— 取 NCBI Taxonomy 的 type material
      与 GenBank 记录的 /type_material。同样内嵌问题标记。

  两列为什么要分开：**菌株号**是"序列库随手记的号"，**模式菌**是"被声明为类型的号"。
  错误恰恰长在两者的落差上 —— 实测本项目最值钱的一条错值（CBS 130417 漏登 / CBS 125378 混入）
  就是"接受名的类型材料槽位被异名的模式号占据"，只有把两列并排才看得见。

问题标记（grep 友好文字，勿改成 emoji——项目内 emoji 已被其它语义占用）
----------------------------------------------------------------------
    号不存在     该号在保藏目录/登记表中查无
    目录记为别种 目录把该号标为另一个种/属
    属异名       该号在 GenBank 的属分布与声明不符
    一株多名     同一保藏号在不同记录里被标了多个物种名
    位数孤例     同批位数与多数不一致
    仅一侧有     只有数据库/只有文献一侧声明
    未核         缺可核验的外部目录，只能标注

框架论文块（可选）
------------------
  --framework <文件>  投喂框架论文（PDF 抽出的 txt，或 PMC XML）
  投喂后额外产出 <genus>_framework.tsv，把论文侧声明与上面两列对照（J-01 四态）。

用法
----
  python -m typeaudit.pipeline --genus Colletotrichum --out out --limit 20
  python -m typeaudit.pipeline --genus Colletotrichum --framework SiM2022.txt --out out
  python -m typeaudit.pipeline --genus Diaporthe --corpus local/ --out out     # 离线语料
"""
import argparse, datetime, hashlib, io, json, os, re, sys, time

# 仓库根：本文件是 <repo>/src/typeaudit/pipeline.py，往上三级。
# ⚠️ 此前这里【没有】REPO 常量，某次改动引用了它却没定义 →
#    write_tables 每次都抛 NameError，而 _flush 的 try/except 把它吞了，
#    日志里只剩一行「增量落盘失败」，metrics.json 从此不再更新（静默）。
REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from collections import Counter, defaultdict

from . import sources as S
from . import rules as R

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

# ---------------------------------------------------------------- 问题标记
FLAG_ABSENT = "号不存在"
FLAG_MISLABEL = "记为别种"
FLAG_GENUS = "属异名"
# ⚠️ 这个标记【不由 check_one 产生】—— 它需要跨号、跨记录的证据
#    （同一保藏号在 GenBank 里被标成两个不同的种），由 tools/seq_axis.py
#    的「同号多名」一级扫描产出，再进「全部发现」的序列层。
#    在此之前它是个【死标记】：定义过、列进过词表、但从未被 append（实测 0 次调用）。
FLAG_MULTINAME = "一株多名"
FLAG_DIGIT = "位数孤例"
FLAG_ONLY_A = "仅数据库有"
FLAG_ONLY_B = "仅文献有"
# ---------- 目级分类查询（用于判定「号已改派」）----------
# 为什么需要：CBS 的低号段（101/102/103/142…）历史上被回收再分配，今天可能指向
# 完全不相干的真菌（如 Schizosaccharomyces、Naganishia）。此时目录那一行描述的是
# 另一株菌，与声明种「无从比对」——判成「目录记为别种」是冤枉保藏中心。
# 判据：两边分类取到【目】级；目不同 ⇒ 不可能是同一物种。
_MB_CLASS_CACHE = {}


def _mb_order(name):
    """学名 → 目名（MycoBank Classification 的第 7 级）；查不到返回空串。"""
    global _MB_CLASS_CACHE
    if not _MB_CLASS_CACHE:
        p = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(
            os.path.abspath(__file__)))), "data", "mb_class.tsv")
        if os.path.exists(p):
            with io.open(p, encoding="utf-8", errors="replace") as f:
                f.readline()
                for line in f:
                    k, _, v = line.rstrip("\n").partition("\t")
                    if k:
                        _MB_CLASS_CACHE[k.strip().lower()] = v
    parts = [x.strip() for x in _MB_CLASS_CACHE.get((name or "").strip().lower(), "").split(",") if x.strip()]
    return parts[6] if len(parts) > 6 else ""


FLAG_UNKNOWN = "未核"
FLAG_NODIR = "无公开数据库"      # 该机构无公开机读目录，无法核验（**不等于号不存在**）
FLAG_UNKNOWN_INST = "未知机构代码"   # 机构代码不在 Biocollections 登记表内 ⇒ **候选错值**（如 CGMG 少了 C）
FLAG_SPECIMEN = "标本号"    # 标本馆号（如 CBS H-21492），不该拿培养物目录核
FLAG_BASIONYM = "基名属"    # 目录记的是基名属（如 Vermicularia→Colletotrichum），非错
# 【机构旧代码别名】写论文时机构改过名/合并过的代码，不该当"疑似编造"。
# 依据（实测）：NCBI Taxonomy 对 Colletotrichum yunnanense 的类型材料同时列了
#   "AS 3.9617" 与 "CGMCC:3.9617" —— 同一管菌两个代码；
#   CGMCC 官网该记录的 history 写着 "←中国科学院微生物研究所"，
#   而 AS = Academia Sinica（中国科学院）旧代码。⇒ AS n 与 CGMCC n 同号。
# 只收有直接证据的映射；没证据的不猜，继续走"疑似非机构编号"。
INST_ALIAS = {
    "AS": "CGMCC",
}
FLAG_NOT_AN_INST = "疑似非机构编号"   # 不在登记表、且不像写错的代码 —— 多半是实验室/宿主标签被误当机构代码
FLAG_NOT_TYPE = "记为非模式菌株"      # 记录把该号当类型材料用，但机构目录标注为非模式菌株（实测 E2）
FLAG_SYNONYM = "异名（非错）"
FLAG_AUTHORITY = "权威分歧"
# 【模式地位不符】记录把某号当模式材料用，而目录给出的类型不是【主模式】。
# 主模式 = 命名法规里真正承载名字的那个（正模/新模式/选模式）；
# 补型、等模式、副模式、权威菌株等都是次要材料 —— 记录笼统写"模式材料"时并不等价。
# 归一化后再比：CBS 写 "Epitype"、ICMP 写 "ex Holotype"，大小写与连字符都不同。
# 主模式 = 真正承载名字的那份材料（含它的衍生培养物）。
def _lev(a, b):
    """两个短词的编辑距离（Levenshtein）。只用来判「像不像拼写错误」。"""
    a, b = (a or "").lower(), (b or "").lower()
    if a == b:
        return 0
    if not a or not b:
        return max(len(a), len(b))
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i]
        for j, cb in enumerate(b, 1):
            cur.append(min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (ca != cb)))
        prev = cur
    return prev[-1]


PRIMARY_TYPES = {"holotype", "exholotype", "neotype", "exneotype",
                 "lectotype", "exlectotype", "type", "extype"}


def _norm_rank(s):
    """把各家的模式类型写法归一：'ex Holotype' → 'exholotype'。"""
    return re.sub(r"[^a-z]", "", str(s or "").lower())


# 模式类型的白话说法 —— 标记里直接写出来，避免「类型层级待确认」这种抽象措辞。
_RANK_ZH = {
    "epitype": "补充模式", "exepitype": "补充模式的衍生培养物",
    "isotype": "等模", "exisotype": "等模的衍生培养物",
    "paratype": "副模", "exparatype": "副模的衍生培养物",
    "syntype": "全模", "exsyntype": "全模的衍生培养物",
    "authenticstrain": "权威菌株（作者鉴定过，非命名模式）",
}


def _rank_zh(raw):
    """给模式类型配一句白话；认不出就返回空串（标记里就只写原文）。"""
    k = _norm_rank(raw)
    if k in _RANK_ZH:
        return _RANK_ZH[k]
    if k.startswith("ex") and k[2:] in _RANK_ZH:
        return _RANK_ZH[k[2:]] + "的衍生培养物"
    return ""
FLAG_TYPE_RANK = "NCBI 记为模式材料"
# 目录记的种与声明种不同，而两者都是【现行接受名】⇒ 不是异名，须查原始发表。
# 与「异名（非错）」的区别：异名是"目录用的是已被并入的旧名"，这一档是两个都还活着。
FLAG_TWO_ACCEPTED = "并名冲突"   # 旧名「两皆接受名」—— 用户裁定改为此名，点明性质（名字该并给谁）
# NCBI 声明为模式号，但保藏机构公开目录查无 —— **不是错值**：
# 多半是「新种已保藏、机构目录尚未公开上线」（实测 14 条里 8 条是 2022 年发表的新种）。
FLAG_TAX_TYPE_ABSENT = "NCBI 声明为模式号，机构无此记录"      # 两家名录说法相反（NCBI 判独立 / MycoBank 判异名）      # 目录名与声明名【同属不同种】，且接受名就是声明种 —— 非错
# 只有这三家有【可核目录】。**错误判断只允许针对它们** ——
# 其余机构（CMM/MFLU/IMI/BRIP…）无从核验，代码写得再可疑也不能判「错」。
VERIFIABLE = ("CBS", "CGMCC", "ICMP")


def _cat_for(up, sources):
    """按机构代码取它的目录对象；没有可查目录的机构返回 None。"""
    if up == "CBS":
        return sources.get("cbs")
    if up == "CGMCC":
        return sources.get("cgmcc")
    if up == "ICMP":
        return sources.get("icmp")
    nc = sources.get("numcat")
    if nc is not None and hasattr(nc, "has_code"):
        try:
            if nc.has_code(up):
                return nc
        except Exception:
            return None
    return None


# 纠错重查的候选域：**严格限定这三家**。
# 放宽到「任何有目录的机构」会让判定太宽松 —— 候选域越宽，越容易冒出
# 「应为 XX」这种看似确定、实则无从验证的结论。（与使用者讨论后收紧。）
CORRECTABLE = ("CBS", "CGMCC", "ICMP")


def queryable_codes(sources):
    """
    纠错重查的候选域 —— 只在 CBS / CGMCC / ICMP 这三家里找。

    ⚠️ 自备目录（data/catalogue_<代码>.tsv，走 numcat）仍可用于**核验号是否存在**，
       但【不进入纠错的候选域】—— 否则每加一个机构，猜测空间就大一圈，判定越来越宽松。
    """
    out = set()
    if sources.get("cbs") is not None:
        out.add("CBS")
    if sources.get("cgmcc") is not None:
        out.add("CGMCC")
    if sources.get("icmp") is not None:
        out.add("ICMP")
    return out


def suggest_number_fix(code, num, sources, bio, max_tries=8):
    """
    纠错重查：号查不到时，看看是不是【机构代码前缀写错了】—— 数字原样放进去查。

    【为什么只做「代码级」纠错，不做「号码级」—— 这是一次实测的教训】
        号码级穷举（漏一位 / 换相邻两位 / 少个点）会给出**假纠正**：
            CBS 16995  → 「应为 CBS 6995」      ← 只是删一位后恰好存在
            ICMP 18578 → 「应为 ICMP 8578」     ← 18578 是真号，8578 只是恰好存在
        原因是**号码越短越可能存在**，「改一位后能查到」不构成证据。
        而代码级纠错有独立证据：该代码**不在**机构登记表里，而候选代码**在**。
        即便这样也只报「疑似」/「应为」——并且要求**用同一个号在候选代码下确实能查到**，
        否则不报（例：CGMG → 近似 CGMCC，但 CGMCC 3.14187 查无 ⇒ 不报纠正）。

    返回 (修正代码, 修正号, 说明) 或 None。**找不到就返回 None，绝不猜。**
    """
    up = re.sub(r"[^A-Z0-9]", "", (code or "").upper())
    qcodes = queryable_codes(sources)
    if up in qcodes:
        # 代码本身就能查 —— 查不到就是查不到，不做任何猜测（数字不会被改动）
        return None
    cands = []
    try:
        # only=qcodes：只在【有可查目录】的机构里找候选（先过滤再排序）
        near = (bio.nearest_codes(up, max_dist=2, min_len=4, limit=8, only=qcodes)
                if bio is not None else [])
    except Exception:
        near = []
    for c in near:
        if c in qcodes:
            # 数字原样带过去 —— 只修前缀，不改数字
            cands.append((c, num, "机构代码前缀应为 {}（该号在它下面确实查得到）".format(c)))
    tried = 0
    for c, v, why in cands:
        if tried >= max_tries:
            break
        cat = _cat_for(c, sources)
        if cat is None:
            continue
        tried += 1
        try:
            st = cat.status_of(c, v)
        except Exception:
            st = None
        if st is not None:
            return c, v, why
    return None


def _binomial_from_status(st):
    """从目录记录里抽出双名。优先用 name 字段，其次从 raw 文本里找 `of Genus species`。"""
    for src in ((st or {}).get("name") or "", (st or {}).get("raw") or ""):
        s = str(src).strip()
        if not s:
            continue
        m = re.search(r"\bof\s+([A-Z][a-z]+)\s+([a-z][a-z\-]+)", s)
        if m:
            return "{} {}".format(m.group(1), m.group(2))
        if re.match(r"^[A-Z][a-z]+\s+[a-z][a-z\-]+$", s):
            return s
    return None


def cbs_dot_rescue(cat, code, num, species, max_pos=6):
    """CBS「号.年」老号【漏点号】的识别 —— **只有能确认种名一致时才给「更正」**。

    背景（用户提出的关键问题）：「号找到了」≠「是同一株」。
    CBS 老号写作 `516.97`（516 号菌株的第 97 支），GenBank 里实测有记录写成 `51697`。
    插一个点号后那个号确实存在 —— 但它可能是**另一个种**的号。
    而 CBS 目录只有 24.1% 的条目带名称标注（status 非空），其余根本没法判。

    三种结论，严格按「我们实际能确认多少」给：
      · 目录那条【带种名且与声明种一致】 → `书写错误（应为 CBS 516.97）`
        这个形式会被 read_trusted 解析、把更正后的号收进受信表 ——
        所以它必须建立在「目录里那个号确实是这个种」之上，缺一不可。
      · 目录那条【带种名但是别的种】     → `号不存在（形近 CBS 516.97 但目录记的是 <名>）`
      · 目录那条【没有种名标注】         → `号不存在（形近 CBS 516.97 但目录无名无法确认）`

    后两种**绝不给更正**，也不进受信表 —— 因为无法排除「引到另一株菌」。
    返回标记字符串；返回 None 表示插点也救不回（仍是「号不存在」）。
    """
    bare = (num or "").strip()
    if not re.match(r"^\d{4,}$", bare):
        return None                      # 太短或含非数字，插点没有意义
    w = (species or "").split()
    sp_gen = w[0].lower() if w else ""
    sp_epi = w[-1].lower() if len(w) > 1 else ""
    for i in range(1, min(len(bare), max_pos + 1)):
        cand = bare[:i] + "." + bare[i:]
        try:
            st2 = cat.status_of(code, cand)
        except Exception:
            continue
        if st2 is None:
            continue
        nm = _binomial_from_status(st2)
        if not nm:
            return "号不存在（形近 {} {} 但该机构无名字可确认）".format(code, cand)
        ws = nm.split()
        g, e = ws[0].lower(), ws[-1].lower()
        if sp_epi and g == sp_gen and e == sp_epi:
            return "书写错误（应为 {} {}）".format(code, cand)
        return "号不存在（形近 {} {} 但该机构记的是 {}）".format(code, cand, nm)
    return None


def looks_like_typo(code, num, bio, max_dist=2, min_len=4):
    """
    判断一个「不在机构登记表里」的代码，更像「代码写错」还是「根本不是机构代码」。

    三条降级判据（任一成立 ⇒ 判为「疑似非机构编号」，降优先级）：
      ① 代码含小写字母 —— 真机构代码惯例全大写；`Paphi` 这种更像一个词（宿主属名）。
      ② 号码部分长度 ≤ 2 —— 真实保藏号通常更长；`2` / `I2` / `01` / `24` 像内部编号。
      ③ 登记表里找不到近似代码 —— 长度 ≥4、编辑距离 ≤2 的邻居一个都没有。

    否则判为「未知机构代码」（真错值候选）。

    ⚠️ 这是启发式，不是判据 —— 结果仍须人判。本机实测那 7 条：
       CGMG 3.14187 → 登记表里有 CGMCC（距离 2）、号码长、全大写 ⇒ 判为真错值候选，
                      与项目记录的 E3 一致；
       Paphi 2-1（含小写）、NK 24 / BPDI 2 / CSSN 2 / LLTA- 01（号 ≤2 位）⇒ 降级。
    """
    if re.search(r"[a-z]", code or ""):
        return False
    if len((num or "").strip()) <= 2:
        return False
    try:
        near = bio.nearest_codes(code, max_dist=max_dist, min_len=min_len, limit=1)
    except Exception:
        near = []
    return bool(near)


# 已知基名属（**快速路径，不是唯一依据**）。
# 真正的判定是按记录问 MycoBank：目录名的接受名落在声明属里即为基名属/异名。
# 本清单只用来省一次查询；它的召回缺口是实测过的 —— 只有 3 个属，
# 而 MycoBank 反推出 Colletotrichum 的同义属候选有 20 个，Gnomonia 那类会被清单漏掉。
BASIONYM_GENERA = {"vermicularia", "gloeosporium", "glomerella"}


def _src_ok(obj):
    """源是否**真有数据**（不是空载）。
    ⚠️ 这是本项目踩过两次的坑：`CBSCatalogue` 在文件缺失时【加载成功但为空】，
       `status_of()` 一律返回 None，调用方就会判「号不存在」——
       实测因此产出过 40 个纯假阳性（第一版），换成"不配目录"的入口又能复现一次。
       **规矩：n_rows()==0 的源一律不得用于判定，只能标「无目录」。**"""
    if obj is None:
        return False
    if hasattr(obj, "n_rows"):
        try:
            return obj.n_rows() > 0
        except Exception:
            return False
    return True          # api 类型（联网）视为可用


def norm_bare(code, number):
    r"""
    pipeline 用的归一键 —— **直接委托给 sources.norm_key，保证只有一份实现**。

    为什么要合并：原来这里是第二份实现，注释写「去分隔符」而代码只剥了国家标签，
    于是 `LLTA- 01` 与 `LLTA 01` 去重失效（用户实测同一个号占两行）。
    后来又发现它和 norm_key 对 `MFLU<THA>` 给出不同的键。
    两份实现 = 两套 bug，所以只留一份。

    键的规则（都在 sources.norm_key 里）：
        MFLU<THA>:090228  ==  MFLU 090228    ← 国家标签 + 冒号
        LLTA- 01          ==  LLTA 01        ← 空格
        CBS 100.06        !=  CBS 10006      ← **保留小数点**：这是两株不同的菌
    """
    return S.norm_key(code, number)


def _fmt(disp, flags):
    """单元格格式：`展示串 [问题1; 问题2]`（没问题时标 —）。
    ⚠️ 必须把 \n 压成空格——`/type_material` 的值常含换行，会撑破 TSV 的一行一记录。"""
    tag = "; ".join(flags) if flags else "—"
    disp = re.sub(r"\s+", " ", str(disp or "")).strip()
    return "{} [{}]".format(disp, tag).strip()


# ---------------------------------------------------------------- 阶段一：acquire
def acquire_species(genus, sources, corpus_dir=None, framework=None, limit=None,
                    email=None, species_list=None):
    """拿到该属的物种清单。四条路，按可靠性排：
       ① --species-list 显式清单（最可靠，推荐给你已经圈定范围的场景）
       ② 框架论文（若投喂）—— 它显式列出了该框架处理过的种
       ③ 本地语料目录（离线）
       ④ GenBank 检索（联网；从 esummary 的 Title 首二词提种名）

    ⚠️ 不要用 `strain_values_of()` 取种名——它返回的是 `/strain` 值（论文内部编号），
       里面全是 `C1315.2`、`LC0034` 这类东西，不是物种名。
    """
    raw_names = []
    if species_list:
        # ⚠️ 宽容三种写法（网页版曾经只能填「文件路径」，用户直接打字就被当成一个种名）：
        #    ① 文件路径       ② 多行文本（每行一个）       ③ 一行里多个种名（空格/逗号分隔）
        txt = ""
        if os.path.exists(species_list):
            txt = io.open(species_list, encoding="utf-8-sig").read()
            src_desc = "文件 " + species_list
        else:
            txt = species_list
            src_desc = "框内文本"
        # 用「属名 + 种加词」把一行里的多个种名切开
        for ln in txt.splitlines():
            ln = ln.strip()
            if not ln or ln.startswith("#"):
                continue
            hits = re.findall(r"\b[A-Z][a-z]+\s+[a-z][a-z\-]+", ln)
            if hits:
                raw_names.extend(hits)
            else:
                raw_names.append(ln)
        print("[acquire] 来自 {}".format(src_desc))
        if raw_names:
            print("[acquire] 解析出 {} 个种名（前 3：{}）".format(
                len(raw_names), ", ".join(raw_names[:3])))
    elif framework and os.path.exists(framework):
        t = io.open(framework, encoding="utf-8", errors="replace").read()
        raw_names = re.findall(re.escape(genus) + r"\s+([a-z][a-z\-]+)", t)
        print("[acquire] 来自框架论文：{} 个".format(len(set(raw_names))))
    elif corpus_dir and os.path.isdir(corpus_dir):
        for f in sorted(os.listdir(corpus_dir)):
            fp = os.path.join(corpus_dir, f)
            if os.path.isfile(fp):
                t = io.open(fp, encoding="utf-8", errors="replace").read()
                raw_names += re.findall(re.escape(genus) + r"\s+([a-z][a-z\-]+)", t)
        print("[acquire] 来自本地语料：{} 个".format(len(set(raw_names))))
    else:
        tx = sources.get("taxonomy")
        if tx is None:
            print("[acquire] ⚠️ 无 taxonomy 源；请用 --species-list / --framework / --corpus")
            return []
        # ⚠️ 分类学清单走 Taxonomy，**不要用 GenBank 标题抠**：为拿"属里有哪些种"
        #    去 efetch 数百条 GenBank 记录，成本高两个数量级（实测：旧路径分钟级、
        #    Taxonomy 路径秒级；旧路径还会把 lineage 名当种名）。
        try:
            raw_names = list(tx.species_in_genus(genus, retmax=3000, limit=None))
            print("[acquire] 来自 NCBI Taxonomy[Subtree]：{} 个种名".format(len(raw_names)))
        except Exception as e:
            print("[acquire] ⚠️ Taxonomy 取种名失败：{}".format(str(e)[:120]))

    # 归一：`Colletotrichum xxx` 形式，去重、排除种加词是 sp/aff/cf 的
    species = sorted({n.strip() for n in raw_names
                      if n and n.split()[-1] not in ("sp", "sp.", "aff", "cf", "x")})
    print("[acquire] 去重后 {} 个种名".format(len(species)))
    if limit:
        species = species[:limit]
        print("[acquire] --limit 生效，只处理前 {} 个".format(len(species)))
    return species


# ---------------------------------------------------------------- 阶段二：collect
def collect(species, genus, sources, note=None):
    """对该种取【菌株号】与【模式菌】两列，并逐号跑检查。

    ⚠️ 接口按 sources.py 的真实返回写（曾按推断写过 4 处，全错）：
      GenBankAPI.culture_collections_of() → list[str]      ← 保藏号，本列主键
      GenBankAPI.strain_values_of()       → list[str]      ← 论文内部编号，仅作旁证
      GenBankAPI.type_material_of()       → list[str]      ← /type_material 声明值
      TaxonomyAPI.type_material_of()      → (taxid, list)  ← 元组，不是 list
      CBSCatalogue.status_of()            → dict|None      ← dict，不是 str
    """
    gb = sources.get("genbank")
    tx = sources.get("taxonomy")
    cbs = sources.get("cbs")
    cbs_ok = _src_ok(cbs)          # 空载目录不得用于判定
    bio = sources.get("biocollections")
    bio_ok = _src_ok(bio)          # 机构登记表（判「未知机构代码」用）
    # ⚠️ species 已是全名（"Colletotrichum arxii"），不要再拼属名——
    #    拼了会得到 "Colletotrichum Colletotrichum arxii"，检索 0 命中。
    full = species if species.lower().startswith(genus.lower()) else "{} {}".format(genus, species)
    notes = []

    # ---- 证据链容器（每个种一份）----
    evid = {}          # 归一化号 → {source, field, query, window, cache, flags}
    fetches = []       # 每一次 GenBank 检索：(库, 检索式, 窗口, 服务端命中数)

    # NCBI Taxonomy 声明的「该种模式材料号」集合 —— 用来区分两种「查无」：
    # ① NCBI 说它是模式号、机构目录却没有 → 多半是新种未公开入库（**不确定**）
    # ② 谁都没说它是模式号、机构目录也没有 → 才是「号不存在」（**错值**）
    tax_type_set = set()

    def _mark(tag, fn, *a, **kw):
        """包住一次取数，把它触发的 esearch 记录进 fetches。

        这是「窗口够不够大」的实测依据 —— 例如某次检索命中 1,917 条而我们只取前 40，
        这里的 retmax 与 count 一对比就看得出来。
        """
        n0 = len(S.QUERY_LOG)
        try:
            return fn(*a, **kw)
        finally:
            for q in S.QUERY_LOG[n0:]:
                fetches.append(dict(q, retrieval=tag))


    # ---------- 第 2 列：菌株号 ----------
    cc_raw, st_raw = [], []
    if gb:
        # ⚠️ /culture_collection 是【稀有字段】：宽搜 200 条里普通分离株占绝大多数，
        #    实测 fructicola/gloeosporioides/siamense 宽搜 200 条含 cc 者为 0；
        #    加 `AND type` 限定后 40 条里有 4 条（10%）。故**两路合并**。
        try:
            cc_raw = list(_mark("菌株号·宽搜", gb.culture_collections_of, full, retmax=300))
        except Exception as e:
            notes.append("取保藏号(宽)失败:{}".format(str(e)[:30]))
        try:
            cc_raw += list(_mark("菌株号·类型材料搜", gb.culture_collections_of, full, retmax=120, extra="type"))
        except Exception as e:
            notes.append("取保藏号(type)失败:{}".format(str(e)[:30]))
        # 【覆盖声明】记录数超过扫描上限的种，明确写出来 —— 宁可说"没扫全"，
        # 也不能让读者以为"这已经是全部号"。
        try:
            _cov = getattr(gb, "last_coverage", None) if gb else None
            if _cov and _cov.get("total", 0) > _cov.get("scanned", 0):
                notes.append("覆盖不全：该种 GenBank 记录 {} 条，本次仅扫 {} 条".format(
                    _cov["total"], _cov["scanned"]))
        except Exception:
            pass
        try:
            st_raw = list(_mark("菌株号·/strain", gb.strain_values_of, full, retmax=200))
        except Exception as e:
            notes.append("取 strain 失败:{}".format(str(e)[:30]))
    else:
        notes.append("无 genbank 源")

    # ---------- 第 3 列：模式菌 ----------
    type_raw = []
    # 【两侧溯源】号是谁声称有的 —— 供「号不存在」时写明"声明方 / 查无方"。
    # 键用 norm_bare（与 tax_type_set 同一口径），值是该号的所有声明来源。
    claim_origin = {}

    def _note_claim(items, label):
        """记录一批号的声明来源（label 形如 'NCBI Taxonomy'）。"""
        for _v in (items or []):
            _m = re.match(r"^\s*([A-Za-z][A-Za-z\.\-]{0,12}(?:<[^>]*>)?)\s*[:#]?\s*([\w\.\-]+)\s*$",
                          str(_v or ""))
            if _m:
                claim_origin.setdefault(norm_bare(_m.group(1), _m.group(2)), set()).add(label)

    def _claim_text(key, fallback=None):
        """把声明来源拼成可读串。

        `fallback` 是「这一列天然来自哪」的兜底说明 —— 例如第 2 列的号来自
        GenBank 记录的 /strain /culture_collection 字段，它不在 claim_origin 里，
        但不该显示成"未标来源"（那会让人以为我们没有溯源）。
        """
        labs = sorted(claim_origin.get(key) or [])
        if labs:
            return "、".join(labs)
        return fallback or "未标来源"
    if tx:
        try:
            res = _mark("模式菌·Taxonomy", tx.type_material_of, full, retmax=3)
            taxid, names = res if isinstance(res, tuple) else (None, res or [])
            type_raw += list(names)
            _note_claim(names, "NCBI Taxonomy")
            # 同时记下「NCBI 说这些号是本种模式材料」，供「查无」时区分成因
            for _v in (names or []):
                _m = re.match(r"^\s*([A-Za-z][A-Za-z\.\-]{0,12}(?:<[^>]*>)?)\s*[:#]?\s*([\w\.\-]+)\s*$",
                              str(_v or ""))
                if _m:
                    tax_type_set.add(norm_bare(_m.group(1), _m.group(2)))
        except Exception as e:
            notes.append("取 Taxonomy 类型材料失败:{}".format(str(e)[:30]))
    else:
        notes.append("无 taxonomy 源")
    if gb:
        try:
            type_raw += list(_mark("模式菌·GenBank /type_material", gb.type_material_of, full, retmax=1500))
        except Exception as e:
            notes.append("取 /type_material 失败:{}".format(str(e)[:30]))
        try:
            # 上面那路返回的是学名字符串（不含号），对第 3 列的数字零贡献。
            # 这一路取同一批记录的 /culture_collection —— 才是第 3 列真正要的号，也是"谁声称有号"的来源。
            _n_before = len(type_raw)
            type_raw += list(_mark("模式菌·类型材料记录的 /culture_collection",
                                   gb.type_material_cc_of, full, retmax=1500))
            _note_claim(type_raw[_n_before:], "GenBank 类型材料记录 /culture_collection")
        except Exception as e:
            notes.append("取类型材料记录的 cc 失败:{}".format(str(e)[:30]))

    # ⚠️ 标记文字里【禁止】出现 ` ｜ `（全角竖线）—— 它是单元格内多个号之间的分隔符，
    #    写进标记会把标记劈成两半，语料匹配随之失效（实测 F30 因此「未出现在表里」）。
    #    需要分隔时用 ` ‖ `。
    # ---------- 逐号检查（**按机构分源**）----------
    # ⚠️ 最容易犯的错：拿 CBS 目录去核所有机构的号 ⇒ ICMP/MFLU/IMI 一律被误判"号不存在"。
    #    本项目实测第一版就因此产出 40 个纯假阳性。**没有目录的机构只能标「无目录」。**
    # ⚠️ 允许 <...>：`MFLU<THA>:090228` 这类带国家标签的号若不被识别，
    #    会走去重兜底路径，与剥标签后的 `MFLU090228` 归不到一起（实测去重失效）。
    RE_CODE = re.compile(r"^\s*([A-Za-z][A-Za-z\.\-]{0,12}(?:<[^>]*>)?)\s*[:#]?\s*([\w\.\-]+)\s*$")
    # 非号的长描述（/type_material 常写成「culture from holotype of X」），一律丢弃
    NOT_A_NUMBER = re.compile(r"(?i)^\s*(culture|from|of|type|ex-?type|isotype|holotype|"
                              r"neotype|epitype|paratype|specimen|voucher|strain|isolate)\b")

    # 【权威的管辖范围】—— 每个标记都要交代「我是问谁、问的哪个字段、窗口多大」。
    # 为什么：`号不存在` 其实是「该机构没有」∧「我查的字段覆盖了它」∧「窗口够大」∧ …
    # 本项目为此反复排查过 7 次（[Strain] 不含 /culture_collection、retmax 截断、
    # 缓存键少一维、归一化不一致、号在同种另一条记录里…）。把这几项写进证据链，
    # 下一次漏检才能一眼定位，而不是又花一轮。
    AUDIT_SOURCES = {
        "CBS": ("CBS 记录（本地 TSV）", "cbs_number → norm_key 精确匹配",
                "全表 62,741 行", "离线文件，不联网"),
        "CGMCC": ("CGMCC 官网 API", "cgmcc_strain_number 精确匹配",
                  "pageSize=10", "缓存 .cache/cgmcc/cgmcc_<号>.json"),
        "ICMP": ("ICMP 官方库（SCD）", "全文检索 query=<号>",
                 "pageNumber=0", "缓存 .cache/icmp/icmp_<号>.json"),
    }
    AUDIT_OTHER = ("机构登记表（Biocollections）", "code 精确匹配",
                   "9,672 条", "离线文件")

    def _code_map(up, num, claimed_type=False, _depth=0):
        """无目录时尝试【跨库映射】：别家机构的号 → CBS 号 → 再用 CBS 目录判一次。

        为什么需要：一管菌常同时存于多家机构（CBS 151555 = BRIP 15842a），
        GenBank 里作者只填手边那个号。不处理映射，「GenBank 无声明 CBS 号」
        会产出过半假阳性（本机实测 6 条里 3 条）。

        ⚠️ 两道闸，缺一不可（都是实测踩出来的）：
          ① 种名一致：映射到的 CBS 记录，其种名必须与声明种相同。
             否则会串物种 —— 实测 `C 824`（声明 alienum）的号串匹配到了
             CBS 115566，而那条记录是 Annulatascus fusiformis，另一个属。
          ② 号完全相等：CBSMap.code_map 内部按「整体相等」匹配，不用「包含」。
             实测 `BRIP 63641` 与 CBS 记录的 `BRIP 63641a` 只差一个尾缀，
             用包含判断会当成同一个号。

        返回 flag 列表（含一条「经映射」来源说明）；返回 None = 无法映射，
        调用方应保持「无目录」。**映射不是结论，是换一个权威来判定。**
        """
        if _depth > 0:
            return None
        sp_epi = species.split()[-1].lower() if species else ""
        for prov in (sources or {}).values():
            fn = getattr(prov, "code_map", None)
            if not callable(fn):
                continue
            try:
                r = fn(up, num)
            except Exception:
                continue
            if r is None:
                continue                     # 这个源没这个能力（不是「查无」）
            for m in (r.get("matches") or []):
                mm = re.search(r"\b([A-Z][a-z]+)\s+([a-z][a-z\-]+)\b", str(m.get("taxon") or ""))
                if not mm or not sp_epi or mm.group(2).lower() != sp_epi:
                    continue                 # 闸①：种名不符，放弃
                _d2, fl2 = check_one("{} {}".format(m["code"], m["number"]),
                                     claimed_type, _depth=1)
                if _d2 is None:
                    continue
                # 【写清两件事】映射到哪管菌 ≠ 那管菌是不是模式材料。
                # 取目标记录自己的记名与模式地位，写进标记：
                # 实测 CBS 112996 的 status 为空 —— 别名链能证明"是同一管菌"，
                # 但 CBS 并没有说它是模式材料。两者混在一起写会误导读者。
                _nm, _stt = "", ""
                try:
                    _cbs = (sources or {}).get("cbs")
                    if _cbs is not None and hasattr(_cbs, "status_of"):
                        _rec = _cbs.status_of(m.get("code"), str(m.get("number"))) or {}
                        # ⚠️ status_of() 返回的键是 rank / raw / unannotated，
                        #    **没有 'status' 这个键**。写成 _rec.get("status") 会永远得到空，
                        #    于是把「Epitype of ...」误报成「模式地位为空」（我自己踩过）。
                        _nm = (_rec.get("name") or "").strip()
                        if _rec.get("unannotated"):
                            _stt = "为空"
                        else:
                            _stt = "=" + ((_rec.get("raw") or _rec.get("rank") or "").strip()[:38])
                except Exception:
                    pass
                # 徽章只留最短的事实；详细叙述交给网页详情页（那里有空间）。
                # 长句塞进小徽章会竖排成一列，反而读不了（用户实测反馈）。
                return list(fl2) + ["经映射→{}".format(m["display"])]
        return None

    def check_one(raw, claimed_type=False, _depth=0, claim_src=None):
        """返回 (展示串, 问题标记列表)。"""
        raw_clean = re.sub(r"\s+", " ", str(raw or "")).strip()
        if not raw_clean or NOT_A_NUMBER.match(raw_clean):
            return None, None                     # 不是号 → 丢弃
        m = RE_CODE.match(raw_clean)
        if not m:
            return None, None
        code, num = m.group(1), m.group(2)
        up = code.upper().rstrip(".")
        fl = []
        # ---- 先解析：这个机构有没有可核验的目录？----
        # CBS 走 cbs 源；其余机构走「通用机构目录」（data/catalogue_<代码>.tsv）。
        # 【为什么要这一步】原来只有 CBS 一份目录，于是 ICMP / CMM / MFLU / IMI 的号
        # 一律标「无目录」= 无从核验。那是**适配器缺口**，不是数据缺陷，也不是「号不存在」。
        # 通用目录一就位，这些机构立刻可核，并且多出一档原来判不了的「号不存在」。
        cgmcc = sources.get("cgmcc") if sources else None
        cat = None
        if up == "CBS":
            cat = cbs if cbs_ok else None
        else:
            # 逐个问「这个机构你有没有目录」：通用机构目录（文件）→ ICMP（官方库）
            for prov in ((sources.get("numcat") if sources else None),
                         (sources.get("icmp") if sources else None)):
                if prov is None:
                    continue
                try:
                    fn = getattr(prov, "has_code", None) or getattr(prov, "claims", None)
                    if fn and fn(up):
                        cat = prov
                        break
                except Exception:
                    continue

        # 标本馆号（CBS H-xxxx 等）：标本号不走培养物目录
        if re.match(r"^H[-\d]", num) and up.startswith("CBS"):
            fl.append(FLAG_SPECIMEN)
        elif up == "CGMCC" and cgmcc is not None:
            # CGMCC 有官网目录可用，不该再标「无目录」。它给的不是名称而是
            # 【模式/非模式】标注 —— 恰好能查 E2 那类错：
            #     记录把该号当类型材料用，而官网标注「非模式菌株」。
            try:
                info = cgmcc.strain_of(num)
            except Exception:
                info = {"found": False, "error": "exception"}
            if not info.get("found"):
                if info.get("total") == 0:
                    # 是 NCBI 声明的模式号吗？是 →「不确定」（可能未公开入库），不是 →「号不存在」
                    _k = norm_bare(code, num)
                    if _k in tax_type_set:
                        fl.append(FLAG_TAX_TYPE_ABSENT)
                    else:
                        fl.append("{}（声明:{} ‖ 查无:CGMCC 官网返回 0 条）".format(
                            FLAG_ABSENT, _claim_text(_k, claim_src)))
                else:
                    fl.append(FLAG_NODIR)
            else:
                ist = (info.get("is_typestrain") or "")
                if claimed_type and ist and ("非模式" in ist):
                    fl.append("{} 记为非模式菌株".format(up))
                elif not ist:
                    fl.append(FLAG_UNKNOWN)
        elif cat is not None:
            # —— 有目录可核 ——
            try:
                st = cat.status_of(code.replace("<", "").replace(">", ""), num)
            except Exception:
                st = None
            if st is None:
                # 【插点号兜底】CBS 老号是「号.年」（如 516.97），实测有记录漏了点号。
                # 分三档、只有能确认种名一致才给「更正」——理由见 cbs_dot_rescue 的 docstring：
                # 「号找到了」≠「是同一株」。
                _dot = cbs_dot_rescue(cat, code, num, species) if up == "CBS" else None
                if _dot:
                    fl.append(_dot)          # 插点号救回的：按它自己的三档报
                elif norm_bare(code, num) in tax_type_set:
                    fl.append(FLAG_TAX_TYPE_ABSENT)
                else:
                    fl.append("{}（声明:{} ‖ 查无:{} 记录）".format(
                        FLAG_ABSENT, _claim_text(norm_bare(code, num), claim_src),
                        "CBS" if up == "CBS" else up))
            elif st.get("ambiguous"):
                # 键有歧义（目录里两个号压成同一个键）⇒ 无从核验，给不出记录就别给
                fl.append(FLAG_UNKNOWN)
            else:
                # 【模式地位不符】只在【记录声称该号是模式材料】时判；
                # 普通菌株本就不该是模式，报它没有意义。
                # ⚠️ 这一段【CBS 与 ICMP 共用】：CBS 的类型在 st["rank"]，
                #    ICMP 在 st["type_status"]（CGMCC 走另一分支的 is_typestrain）。
                #    只读 rank 会把 ICMP 的声明整个漏掉（用户指出过）。
                if claimed_type:
                    _rk_raw = (st.get("rank") or st.get("type_status") or "").strip()
                    if _rk_raw and _norm_rank(_rk_raw) not in PRIMARY_TYPES:
                        _zh = _rank_zh(_rk_raw)
                        # 【必须核对记录原文的等级】——只看"号出现在类型材料字段"不够。
                        # 全属实测 28 条这类里：10 条记录原文也是 epitype（本就一致）、
                        # 15 条记录根本没有 /type_material 字段（标记不成立）、
                        # 2 条记录写的是 holotype（真冲突）、1 条只笼统写 type material。
                        _tm_texts = []
                        # ⚠️ 这里【不能】用 sp_bare —— 它在下面几行才定义。
                        #    实测：用了就抛 NameError，被 except 吞掉 → 28 条全变「待核」，
                        #    看起来像"记录都没声明"，其实是自己代码错了。
                        _sp_bare_local = (species.split()[-1].lower() if species else "")
                        try:
                            _gb = sources.get("genbank") if sources else None
                            if _gb is not None and hasattr(_gb, "type_material_text"):
                                _recs = _gb.type_material_text(up, num) or []
                                _mine = [r["text"] for r in _recs
                                         if _sp_bare_local and
                                         _sp_bare_local in (r.get("organism") or "").lower()]
                                _tm_texts = _mine or [r["text"] for r in _recs]
                                if not _tm_texts and _recs:
                                    _tm_texts = [r["text"] for r in _recs]
                        except Exception as _e:
                            _tm_texts = []
                        _tm = " ".join(_tm_texts).lower()
                        if not _tm_texts:
                            fl.append("模式地位待核（{} 记为 {}，但 GenBank 记录未声明类型等级）".format(
                                up, _rk_raw))
                        elif "epitype" in _tm:
                            fl.append("记录原文亦为 epitype，与 {} 记为 {} 一致（非错）".format(
                                up, _rk_raw))
                        elif re.search(r"\bholotype\b|\bneotype\b|\blectotype\b", _tm):
                            fl.append("模式地位冲突（记录称正模，{} 记为 {}）".format(up, _rk_raw))
                        else:
                            fl.append("模式地位不精确（记录只写 type material，{} 记为 {}）".format(
                                up, _rk_raw))
                nm = (st.get("name") or "").lower()
                sp_bare = species.split()[-1].lower() if species else ""
                gen = (nm.split()[0] if nm.split() else "")
                if gen in BASIONYM_GENERA:
                    fl.append(FLAG_BASIONYM)       # 已知基名属，非错（快速路径）
                elif nm and sp_bare and sp_bare not in nm:
                    # 目录记的名与声明种不符。先别急着判「别种」——去问 MycoBank：
                    # 这个名字的【接受名链】有没有走到声明属？
                    #   走到了 ⇒ 它只是老属名（基名属/异名），非错。
                    #
                    # 【为什么把硬编码清单换成按记录问 MycoBank】
                    #     BASIONYM_GENERA 只有 3 个属，而 MycoBank 反推出 Colletotrichum
                    #     的同义属候选有 20 个（≥2 次）乃至 53 个。手工清单必然有召回缺口。
                    #     按记录问 MycoBank 是权威依据，不依赖任何手工清单。
                    # 【为什么必须走到不动点】
                    #     MycoBank 的 Current name 不一定是终点。实测链条：
                    #         Gnomonia cingulata → Glomerella cingulata var. cingulata
                    #                            → Colletotrichum gloeosporioides
                    #     只走一步会把 CBS 152.26 误判成「目录记为别种」（实测 9 条 vs 走完 10 条）。
                    acc = ""
                    myb_ok = False       # MycoBank 是否真的可用 —— 决定能不能下「别种」的结论
                    try:
                        myb = sources.get("mycobank") if sources else None
                        # ⚠️ 必须走 _src_ok：空载的 MycoBank 让 current_name 一律返回 None，
                        #    acc 因此为空，落到「目录记为别种」——把异名误判成错值。
                        #    与 CBS 空载目录曾产出 40 个假阳性是同一类坑（见 _src_ok 说明）。
                        if _src_ok(myb):
                            myb_ok = True
                            cur = " ".join(nm.split()[:2])
                            for _ in range(6):
                                nxt = (myb.current_name(cur) or "").strip()
                                nxt = " ".join(nxt.split()[:2]).lower()   # 去作者再查下一跳
                                if not nxt or nxt == cur.lower():
                                    break
                                cur = nxt
                            acc = cur.lower()
                    except Exception:
                        acc = ""
                    # ---- 两家名录对账（NCBI Taxonomy ↔ MycoBank）----
                    # 为什么先做这一步：异名判定过去只问 MycoBank 一家。两家说法相反时
                    # （实测 Colletotrichum analogum / celtidis：NCBI 各给独立 taxid，
                    #   MycoBank 判并入 camelliae / karsti），单问一家会把分歧误判成错值。
                    _auth = None
                    try:
                        if nm and species and tx is not None:
                            _auth = S.authority_verdict(nm, species, tx, myb)
                    except Exception:
                        _auth = None
                    if _auth and _auth.get("verdict") == "分歧":
                        # 「两边相反」是客观可验证的 ⇒ 判错值；
                        # 但「谁错」需第三方证据，所以标记里写明「方向待裁决」。
                        fl.append("{}（NCBI：{} ｜ MycoBank：{}；方向待裁决）".format(
                            FLAG_AUTHORITY, _auth.get("ncbi_label"), _auth.get("mycobank_label")))
                    elif _auth and _auth.get("verdict") == "无法对账" and myb_ok and not acc:
                        fl.append(FLAG_UNKNOWN)
                    elif acc and genus and acc.startswith(genus.lower() + " "):
                        # 接受名落在声明属里 ⇒ 非错。但必须分清是哪种非错：
                        #   目录名的属 == 声明属 → **异名**（同属不同种）
                        #        例：Colletotrichum hymenocallidis → C. siamense
                        #   目录名的属 != 声明属 → **基名属**（老属名）
                        #        例：Vermicularia capsici → Colletotrichum truncatum
                        # ⚠️ 这两者绝不能混：把同属异名标成「基名属」是**事实错误**
                        #    （「基名属」字面意思是目录记的是另一个属）。本机实测踩过：
                        #    CBS 125378 的目录名 Colletotrichum hymenocallidis 与声明名同属，
                        #    却一度被标成「基名属」。
                        # ⚠️ 还必须比【接受名的种加词】与声明种是否相同。
                        #    原来只判「接受名落在声明属里 ⇒ 异名」，于是
                        #    目录记 gloeosporioides、接受名也是 gloeosporioides、
                        #    而记录声明的是 alienum —— 三个都是现行接受名，
                        #    却被判成「异名（非错）」（实测 alienum/arboricola/boninense 共 6 例）。
                        #    现在与网页的 name_relation 用同一口径。
                        _acc_parts = acc.split()
                        _acc_epi = (_acc_parts[1] if len(_acc_parts) > 1 else "")
                        if _acc_epi and sp_bare and _acc_epi == sp_bare:
                            # 接受名就是声明种 ⇒ 目录名只是它的旧名/基名
                            fl.append(FLAG_SYNONYM if (gen and gen == (genus or "").lower())
                                      else FLAG_BASIONYM)
                        else:
                            _cat_nm = (st.get("name") or "")
                            _cp = _cat_nm.split()
                            # 【取"基础种加词"】跳过 f./var./subsp. 这类种下单元标记，
                            # 取第一个全小写的词。这一步解决三类误判（全属实测）：
                            #   · Colletotrichum acutatum f. pineum  → 基础种加词 acutatum
                            #   · Colletotrichum dematium var. minus → 基础种加词 dematium
                            #   · Colletotrichum Corda（只有属名+定名人）→ 取不到 → 记「只到属名」
                            _base_epi = ""
                            for _w in _cp[1:]:
                                if re.match(r"^[a-z][a-z\-]+$", _w):
                                    _base_epi = _w.lower()
                                    break
                            _is_infra = bool(re.search(r"\b(f\.|f\.sp\.|var\.|subsp\.|forma)\s",
                                                       _cat_nm + " ", re.I))
                            # 拼写变体：MycoBank 的 Name status 会明确写 Orthographic variant
                            _st_myb = ""
                            try:
                                _st_myb = (myb.status_of(_cat_nm) or "").strip() if myb is not None else ""
                            except Exception:
                                _st_myb = ""
                            _is_variant = _st_myb.lower().startswith("orthographic")

                            if not _base_epi:
                                # 目录只记到属名（第二词是定名人，不是种加词）⇒ 没法比种
                                fl.append("{} 只记到属名（{}），无法比对种".format(up, _cat_nm))
                            elif _base_epi == sp_bare:
                                # 基础种加词就是声明种 ⇒ 只是拼写/等级/异名的差别，非错
                                if _is_variant:
                                    fl.append("拼写变体（{} 记 {}，正名 {}，非错）".format(
                                        up, _cat_nm, acc or _base_epi))
                                elif _is_infra:
                                    fl.append("种下单元（{} 记为 {}，非错）".format(up, _cat_nm))
                                else:
                                    fl.append(FLAG_SYNONYM)
                            else:
                                # 真·两个不同的种。把"目录具体记的是什么"写清楚，
                                # 免得读者以为只是拼写差异（用户实测指出这一点）。
                                _extra = ""
                                if _is_infra:
                                    _extra = "（种下单元，基础种 {}）".format(_base_epi)
                                elif _is_variant:
                                    _extra = "（{} 是 {} 的拼写变体）".format(
                                        _cp[1] if len(_cp) > 1 else _base_epi, acc or "?")
                                # 若目录名在 MycoBank 里【根本查不到】，且与声明种只差 1–2 个字母
                                # ⇒ 很可能是拼写错误（实测 CBS 记 "protaea"，MycoBank 只有 "proteae"，
                                #    而 protaea 在 MycoBank 里查无此名）。
                                try:
                                    _in_myb = bool(myb and (myb.current_name(_cat_nm) or
                                                            (myb.status_of(_cat_nm) or "")))
                                except Exception:
                                    _in_myb = True
                                if (not _in_myb) and sp_bare:
                                    _d = _lev(_base_epi, sp_bare)
                                    if _d <= 2:
                                        _extra = "（MycoBank 查无此名，与声明种差 {} 个字母，疑拼写错误）".format(_d)
                                if _acc_epi == _base_epi:
                                    fl.append("{}（须查原始出处）‖{} 记 {}{}，MycoBank 亦接受该名".format(
                                        FLAG_TWO_ACCEPTED, up, _cat_nm, _extra))
                                else:
                                    fl.append("{}（须查原始出处）‖{} 记 {}{}，MycoBank 接受名 {}".format(
                                        FLAG_TWO_ACCEPTED, up, _cat_nm, _extra, acc))
                    elif myb_ok:
                        # 先判「号已改派」：目录所指生物与声明种不在同一个目
                        # ⇒ 两条记录根本没在说同一件事，不是记录冲突（非错）。
                        # 实测 CBS 101→Schizosaccharomyces japonicus、CBS 142→Naganishia albida。
                        _oc, _os = _mb_order(nm), _mb_order(species)
                        # 有性型/异名：目录名的现名（或声明种的现名）正好是对方
                        # 实测 Hypocrea crassa ↔ Trichoderma neocrassum（同一物种，两个属名）
                        _cur_sp = ""
                        try:
                            _cur_sp = " ".join((myb.current_name(species) or "").split()[:2]).lower()
                        except Exception:
                            _cur_sp = ""
                        _nm2 = " ".join(nm.split()[:2]).lower()
                        _pretty = " ".join([nm.split()[0].capitalize()] + nm.split()[1:]) if nm else nm
                        if _oc and _os and _oc != _os:
                            fl.append("{} 号已改派（目录现记为 {}，与声明种目级不同，非错）".format(
                                up, _pretty))
                        elif _cur_sp and _cur_sp == _nm2:
                            fl.append("{} 有性型/异名（目录记 {}，该种现名 {}，非错）".format(
                                up, _pretty, _cur_sp.split()[0].capitalize() + " " + _cur_sp.split()[1]
                                if len(_cur_sp.split()) > 1 else _cur_sp))
                        else:
                            fl.append("{} 记为别种（目录记 {}）".format(up, _pretty))
                    else:
                        # MycoBank 不可用 ⇒ 无法判断是异名还是别的种。
                        # 此时说「目录记为别种」是无据断言，只能标「未核」。
                        fl.append(FLAG_UNKNOWN)
                if st.get("unannotated"):
                    fl.append(FLAG_UNKNOWN)
                elif claimed_type and ("type_status" in st) and not (st.get("type_status") or "").strip():
                    # 记录把它当类型材料用，但该目录【没有】模式地位标注 ⇒ 无从核验。
                    # 注意：这和「非模式菌株」不是一回事 —— 目录没标 ≠ 不是模式。
                    fl.append(FLAG_UNKNOWN)
        elif up != "CBS":
            # 【第一步必须是查机构登记表】—— 机构真实存在（CMM/MFLU/IMI…）但没有机读目录
            # ⇒ 无从核验，标「无目录」。
            # ⚠️ 这一步必须在纠错之前！本机实测把顺序搞反后：
            #    CMM 是登记表里的机构，却因 _lev("CMM","ICMP")=2 被纠成「书写错误（应为 ICMP 3814）」
            #    —— 荒谬，而且更正号会被写进受信表，污染下游。
            #    「形近」只对【根本不在登记表里】的代码才有意义。
            typ = None
            if bio_ok:
                try:
                    typ = bio.institution_type(up)
                except Exception:
                    typ = None
            if typ:
                # 有登记、无目录 ⇒ 先试跨库映射；映射不到才标「无目录」
                _mapped = _code_map(up, num, claimed_type)
                if _mapped:
                    fl.extend(_mapped)
                else:
                    fl.append(FLAG_NODIR)
            else:
                # 代码确实不在登记表里。三条路：
                #   ① 形近 CBS/CGMCC/ICMP 之一，且【原号码】在更正后代码下确实查得到
                #      ⇒ 书写错误（已证实），给出更正号；更正号会被 read_trusted 收进受信表
                #   ② 形近三家之一，但那个号在它下面查无
                #      ⇒ 仍判书写错误（代码确非登记机构），只是给不出更正号
                #   ③ 既不形近三家、也不在登记表里
                #      ⇒ 降级为「疑似非机构编号」（不确定）。**不判错**——
                #        判了也无从验证，这正是候选域只留三家的原因。
                # ⚠️ 纠错**只修代码前缀，数字原样带过去**。号码级穷举会给出假纠正
                #    （实测 CBS 16995 → 6995、ICMP 18578 → 8578，后者本来就是真号）。
                qc = queryable_codes(sources)
                fix = suggest_number_fix(code, num, sources, bio)
                try:
                    near3 = (bio.nearest_codes(code, max_dist=2, min_len=4, limit=4, only=qc)
                             if (bio_ok and qc) else [])
                except Exception:
                    near3 = []
                if fix:
                    fl.append("书写错误（应为 {} {}）".format(fix[0], fix[1]))
                elif near3:
                    fl.append("书写错误（形近 {}，但该号在 {} 查无）".format(near3[0], near3[0]))
                else:
                    # 先看是不是机构旧代码（见 INST_ALIAS 的依据）
                    _ali = INST_ALIAS.get(up)
                    _done = False
                    if _ali:
                        try:
                            _d2, _fl2 = check_one("{} {}".format(_ali, num),
                                                  claimed_type, _depth=_depth + 1)
                            _txt2 = " ".join(_fl2 or [])
                            if _d2 is not None and not any(
                                    k in _txt2 for k in
                                    ("疑似非机构编号", "未知机构代码", FLAG_ABSENT)):
                                fl.append("机构旧代码（{} 为 {} 前身）".format(up, _ali))
                                fl.extend(_fl2 or [])
                                _done = True
                        except Exception:
                            pass
                    if not _done:
                        fl.append(FLAG_NOT_AN_INST)
        else:
            fl.append(FLAG_NODIR)
        # ---- 记录「这个结论依赖哪些合取项」----
        src, fld, win, cch = AUDIT_SOURCES.get(up, AUDIT_OTHER)
        evid[norm_bare(code, num)] = {
            "display": "{} {}".format(code, num),
            "source": src, "field": fld, "query": norm_bare(code, num),
            "window": win, "cache": cch,
            "catalog_loaded": (type(cat).__name__ if cat is not None else None),
            "claimed_type": bool(claimed_type),
            # 这条结论是不是靠跨库映射才得出的 —— 证据链的新合取项
            "mapped": any(str(x).startswith("经映射") for x in fl),
            "flags": list(fl),
        }
        return "{} {}".format(code, num), fl

    def _uniq_raw(vals):
        """按【归一键】把同一号的多种拼法收敛成一个代表值 —— **必须在判定之前做**。

        ⚠️ 为什么（本机实测 Colletotrichum adenophorae，影响面 4/482 格）：
          NCBI Taxonomy 把同一个号写成 4 种拼法：
              CGMCC 3.18952 / CGMCC-3.18952 / CGMCC3.18952 / CGMCC:3.18952
          旧顺序是「先逐个 check_one、再用 dedup 合并标记」，而 dedup 取的是**标记的并集**
          ⇒ 同一个号上挂了互相矛盾的两条结论：
              CGMCC 3.18952 [书写错误（形近 CGMCC，但该号在 CGMCC 查无）; 号不存在]
          一个号配两条结论，审稿时读不了；metrics 的 flag_counts 也被重复计数。
          正确顺序：先归一化收敛成 1 个代表值，再判一次。

        代表值的选法（确定性）：优先「纯字母代码 + 空格分隔」的规范形态。
        这不只是好看 —— `CGMCC-3.18952` 这种拼法会被解析成【别的代码】，
        判定结果因此不同；选规范形态同时保证了判定用的是最可信的那一次解析。
        """
        def score(v):
            s = re.sub(r"\s+", " ", str(v or "")).strip()
            m = RE_CODE.match(s)
            if not m:
                return (0, -len(s))
            letters = 1 if (m.group(1) or "").isalpha() else 0
            spaced = 1 if re.match(r"^[A-Za-z]+\s+\S", s) else 0
            return (2 + letters + spaced, -len(s))

        chosen = {}
        for v in vals:
            s = re.sub(r"\s+", " ", str(v or "")).strip()
            m = RE_CODE.match(s)
            k = norm_bare(m.group(1), m.group(2)) if m else re.sub(r"[\s:<>]", "", s).upper()
            if k not in chosen or score(v) > score(chosen[k]):
                chosen[k] = v
        return list(chosen.values())


    def dedup(pairs):
        """按【归一化键】去重（`CBS 132511` 与 `CBS:132511` 同一个号）。"""
        seen, out = {}, []
        for disp, fl in pairs:
            if disp is None:
                continue
            m = RE_CODE.match(disp.replace("(strain) ", ""))
            key = norm_bare(m.group(1), m.group(2)) if m else re.sub(r"[\s:<>]", "", disp).upper()
            if key in seen:
                old = seen[key]
                merged = sorted((set(old[1]) | set(fl)) - {FLAG_NODIR} or (set(old[1]) | set(fl)))
                old[1] = merged
                continue
            seen[key] = [disp, fl]
            out.append(seen[key])
        return [(d, f) for d, f in out]


    # ⚠️ 先 _uniq_raw（按归一键收敛拼法）再 check_one —— 顺序不能反，理由见 _uniq_raw。
    col_strain = dedup([check_one(x, claimed_type=False,
                        claim_src="GenBank 记录 /strain·/culture_collection")
                       for x in sorted(_uniq_raw(cc_raw))])
    # ⚠️ 原来只取前 3 个，实测会丢号：C. lini 的 /strain 有 15 个值，
    #    排序后 392-1、392-2 排在第 5、6 位 → 被截断，与旧表比对时成了"遗漏"。
    #    放到 6 个仍保留列的可读性，同时把"还有多少没显示"记进 notes，不静默截断。
    for x in sorted(set(st_raw))[:6]:
        col_strain.append(("(strain) " + x, []))
    col_type = dedup([check_one(x, claimed_type=True,
                       claim_src="类型材料声明")
                      for x in sorted(_uniq_raw(type_raw))])


    # ---------- 一株多名（同一保藏号被标多个种名）----------
    multiname = []
    return {
        "species": species,
        "col_strain": col_strain,
        "col_type": col_type,
        "notes": notes,
        "evidence": {"numbers": evid, "fetches": fetches},
    }


# ---------------------------------------------------------------- 阶段三：framework
def framework_columns(path, genus):
    """从框架论文抽出第三源：图标签 + 正文 Typus 段声明。"""
    if not path or not os.path.exists(path):
        return {}
    raw = io.open(path, encoding="utf-8", errors="replace").read()
    out = {}
    # 图标签：形如 <species> <CODE> <NUM>
    for m in re.finditer(re.escape(genus) + r"\s+([a-z][a-z\-]+)[^\n]{0,40}?([A-Z]{2,8})[\s:]?([\w\.\-]{3,})", raw):
        out.setdefault(m.group(1), []).append("{} {}".format(m.group(2), m.group(3)))
    # 正文 Typus 段
    for m in re.finditer(r"(?i)typus[^\n]{0,200}", raw):
        seg = m.group(0)
        sp = re.search(re.escape(genus) + r"\s+([a-z][a-z\-]+)", seg)
        if sp:
            out.setdefault(sp.group(1), []).extend(
                ["{} {}".format(a, b) for a, b in
                 re.findall(r"([A-Z]{2,8})[\s:]?([\w\.\-]{3,})", seg)])
    return {k: sorted(set(v)) for k, v in out.items()}


# ---------------------------------------------------------------- 阶段四：report
def write_tables(rows, out_dir, genus, fw=None, quiet=False):
    os.makedirs(out_dir, exist_ok=True)
    g = genus.lower()

    def fmt(cell):
        if not cell:
            return "—"
        # 单元格里的换行会拆行、把所有列错位（seq_axis 里实测踩过 1 例）。
        # 这里统一压成空格，属于"写档前必做"的清洗，不影响判定。
        s = " ｜ ".join(_fmt(disp, fl) for disp, fl in cell)
        return re.sub(r"[\s\u3000]+", " ", s).strip()

    p1 = os.path.join(out_dir, "{}_三列表.tsv".format(g))

    # 【护栏】部分运行不得覆盖完整结果表。
    # 实测事故：网页端「限制种数」默认 3，有人把输出目录指到 results/Colletotrichum
    # 后点「开始」→ 用 3 行覆盖了 408 行的正式结果表（文件从 136,825 字节掉到 678 字节）。
    # 这里改为：若新结果比已有表更少，写到 <name>.partial.tsv，原表一个字节都不动。
    _existing = 0
    if os.path.exists(p1):
        try:
            with io.open(p1, encoding="utf-8-sig") as _f:
                _existing = max(0, len(_f.read().splitlines()) - 1)
        except Exception:
            _existing = 0
    if _existing and len(rows) < _existing:
        p1 = os.path.join(out_dir, "{}_三列表.partial.tsv".format(g))
        if not quiet:
            print("[report] ⚠️ 本次只产出 {} 行，少于已存在的 {} 行 —— "
                  "为免覆盖完整结果表，改写 partial：{}".format(len(rows), _existing, p1))

    with io.open(p1, "w", encoding="utf-8-sig", newline="") as f:
        f.write("种名\t菌株号（含问题）\t模式菌（含问题）\n")
        for r in rows:
            f.write("{}\t{}\t{}\n".format(r["species"], fmt(r["col_strain"]), fmt(r["col_type"])))
    if not quiet:
        print("[report] 三列表 → {}  （{} 行）".format(p1, len(rows)))

    # ---- 覆盖声明：哪些种没扫全 ----
    # 为什么单独一张表而不是塞进三列表：三列表每格必须是「号＋[标记]」，
    # 塞非号内容会把网页卡片解析搅乱。覆盖情况是"这张表可信到什么程度"的元信息，
    # 单列一表最清楚，也不动主表结构。
    _cov_rows = []
    for r in rows:
        for _n in (r.get("notes") or []):
            if str(_n).startswith("覆盖不全"):
                _cov_rows.append((r["species"], str(_n)))
                break
    if _cov_rows:
        p3 = os.path.join(out_dir, "{}_覆盖声明.tsv".format(g))
        with io.open(p3, "w", encoding="utf-8-sig", newline="") as f:
            f.write("种名\t覆盖情况\n")
            for sp, n in _cov_rows:
                f.write("{}\t{}\n".format(sp, n))
        if not quiet:
            print("[report] 覆盖声明 → {} （{} 个种记录数超过扫描上限）".format(p3, len(_cov_rows)))

    if fw:
        p2 = os.path.join(out_dir, "{}_framework.tsv".format(g))
        with io.open(p2, "w", encoding="utf-8-sig", newline="") as f:
            f.write("种名\t框架侧声明\t数据库侧（菌株号）\t数据库侧（模式菌）\t对账\n")
            for r in rows:
                sp = r["species"].split()[-1]
                fwv = sorted(set(fw.get(sp, [])))
                db = sorted({c for c, _ in r["col_strain"] + r["col_type"]})
                if not fwv and not db:
                    v = "两侧均无"
                elif not fwv:
                    v = "仅数据库有"
                elif not db:
                    v = "仅框架有"
                elif set(fwv) & set(db):
                    v = "一致"
                else:
                    v = "不一致"
                f.write("{}\t{}\t{}\t{}\t{}\n".format(
                    r["species"], " ｜ ".join(fwv) or "—",
                    fmt(r["col_strain"]), fmt(r["col_type"]), v))
        if not quiet:
            print("[report] 框架对照表 → {}（J-01 四态：一致/不一致/仅一侧）".format(p2))

    p3 = os.path.join(out_dir, "{}_metrics.json".format(g))
    allc = [c for r in rows for c, _ in r["col_strain"] + r["col_type"]]
    allf = [f for r in rows for _, fl in r["col_strain"] + r["col_type"] for f in fl]
    # 【可复现性元数据】审稿人要复现时必须能答：哪天跑的、代码哪版、数据哪版。
    # 缺了这三样，跑出来的数字只是"某次运行的结果"，无法复现。
    def _fhash(p):
        try:
            h = hashlib.sha256()
            h.update(io.open(p, "rb").read())
            return h.hexdigest()[:12], datetime.datetime.fromtimestamp(
                os.path.getmtime(p)).strftime("%Y-%m-%d %H:%M")
        except Exception:
            return "", ""

    _src = os.path.join(os.path.dirname(os.path.abspath(__file__)), "pipeline.py")
    _vhash, _vtime = _fhash(_src)
    _snap_files = {}
    for _name in ("cbs_catalogue.tsv", "mycobank_full.tsv", "biocollections_codes.csv"):
        _q = os.path.join(REPO, "data", _name)
        if os.path.exists(_q):
            _hh, _tt = _fhash(_q)
            _snap_files[_name] = {"sha256_12": _hh, "mtime": _tt, "bytes": os.path.getsize(_q)}

    try:
        json.dump({
        "genus": genus, "n_species": len(rows), "n_numbers": len(allc),
        "generated": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "code_version": {"file": "src/typeaudit/pipeline.py", "sha256_12": _vhash, "mtime": _vtime},
        "data_snapshot": _snap_files,
        "denominator": {"primary": "n_numbers", "note": "错误率以号数（n_numbers）为分母；种数 n_species 作辅助"},
        "flag_counts": dict(Counter(allf)),
        "conventions": {
            "col2": "菌株号 = 序列库记录的 /culture_collection + /strain",
            "col3": "模式菌 = NCBI Taxonomy type material + GenBank /type_material",
            "flag_vocab": [FLAG_ABSENT, FLAG_MISLABEL, FLAG_UNKNOWN_INST, FLAG_GENUS,
                           FLAG_MULTINAME, FLAG_DIGIT, FLAG_ONLY_A, FLAG_ONLY_B,
                           FLAG_NODIR, FLAG_SPECIMEN, FLAG_BASIONYM, FLAG_SYNONYM,
                           FLAG_NOT_AN_INST, FLAG_NOT_TYPE, FLAG_UNKNOWN],
            "note": "机构词典不全时只漏检不误报 ⇒ 任何计数均为【下界】",
            "col2_coverage": ("第 2 列只覆盖 /culture_collection 字段存在的记录。"
                              "该字段在宽搜记录里罕见（实测 200 条里：fructicola/gloeosporioides/"
                              "siamense 均为 0，falcatum 2，arxii 18），故本列对多数种是【下界】，"
                              "不宜用『该列表为空』反推『该种无保藏号』。"),
            "col3_note": ("第 3 列含非保藏号写法（标本馆号 `CBS H-xxxx` 标『标本号』；"
                          "`Paphi 2-1` 这类个人编号标『无目录』）。"),
            "flag_meaning": {
                "号不存在": "该号在对应机构的公开记录里确实查无（仅在该机构有公开数据时判）",
                "目录记为别种": "该机构把该号标为另一个种/属（基名属已排除）",
                "基名属": "该机构记的是基名属（老属名，如 Vermicularia→Colletotrichum），非错",
                "异名（非错）": "机构记的名与声明名【同属不同种】，且 MycoBank 接受名就是声明种"
                                "（如 Colletotrichum hymenocallidis → C. siamense），非错",
                "标本号": "标本馆号，不该拿菌种保藏记录去核",
                "无目录": "该机构没有公开的机读数据库 ⇒ 无法核验（**不等于号不存在**）",
                "未核": "号存在，但该机构没标它的模式地位",
                "未知机构代码": "机构代码不在 Biocollections 登记表内 ⇒ **候选错值**（如 CGMG 少了 C）；但登记表未必收录齐全，须人判",
                "目录记为非模式": "记录把该号当类型材料用，但该机构（如 CGMCC 官网）标注为非模式菌株 ⇒ 候选错值",
                "疑似非机构编号": "不在登记表、且不像写错的代码（含小写字母 / 号码≤2位 / 无近似代码）"
                                  "⇒ 多半是实验室编号或宿主标签被误当机构代码，**优先看「未知机构代码」那一档**",
                "标本号": "标本馆号（CBS H-xxxx），不走菌种保藏记录",
                "基名属": "该机构记的是基名属，非错",
            },
        },
    }, io.open(p3, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    except Exception as _e:
        # 不许静默：元数据写不出来 = 这份产物【不可复现】，必须让人看见
        print("  ❌ metrics.json 写入失败（{}: {}）—— 本批产物缺少可复现性元数据".format(type(_e).__name__, _e))
    if not quiet:
        print("[report] 指标 → {}（含口径 conventions）".format(p3))

    # ---- 证据链（sidecar）----
    # 为什么不写进三列表：三列表是给人读的结论表，列一旦增加，所有下游脚本
    # 与已发基线全部要跟着改。合取项属于「审计过程」，单独落 JSONL，
    # 需要时按号 join 回去 —— 结论表一个字节都不动。
    p4 = os.path.join(out_dir, "{}_证据.jsonl".format(g))
    with io.open(p4, "w", encoding="utf-8", newline="\n") as f:
        for r in rows:
            ev = (r.get("evidence") or {})
            f.write(json.dumps({
                "species": r.get("species"),
                # 每个号：源 / 字段 / 检索式 / 窗口 / 缓存去向 / 结论标记
                "numbers": ev.get("numbers") or {},
                # 每次 GenBank 取数：库 / 检索式 / 窗口 / 服务端命中数
                "fetches": ev.get("fetches") or [],
                # 取数过程中被吞掉的异常 —— 从前只塞进返回字典、从不输出，
                # 于是「取数失败」与「真的查无」在结果里长得一模一样
                "notes": r.get("notes") or [],
            }, ensure_ascii=False) + "\n")
    if not quiet:
        print("[report] 证据链 → {}（每个标记依赖的源/字段/检索式/窗口/缓存）".format(p4))
    return p1, p2 if fw else None, p3, p4


def write_progress(out_dir, genus, done, total, finished, out_tsv=None):
    """写进度边车文件 `<属>_进度.json`。

    为什么必须有它（本机实测）：三列表本身「长得一样」—— 跑了 50/418 个种就停下的表，
    和「全属都查过了」的表，**格式与内容无从分辨**。不标注就会被当成完整结论。
    这跟前两轮的「暂无法判定」是同一类诚实边界问题：
    **部分结果必须自己声明是部分的，不能靠人记得。**
    """
    g = (genus or "").lower()
    os.makedirs(out_dir, exist_ok=True)
    p = os.path.join(out_dir, "{}_进度.json".format(g))
    d = {
        "genus": genus,
        "done": int(done), "total": int(total),
        "finished": bool(finished),
        "updated": time.strftime("%Y-%m-%d %H:%M:%S"),
        "out_tsv": out_tsv,
        "_说明": ("finished=false 表示这是【部分结果】：只跑了 done/total 个种。"
                  "不要把它当成全属结论；接着跑不会重复下载已完成的种。"),
    }
    try:
        json.dump(d, io.open(p, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    except Exception:
        pass
    return p


# ---------------------------------------------------------------- 主流程
def run(genus, out_dir="out", config=None, corpus_dir=None, framework=None,
        limit=None, email=None, live=True, species_list=None, ncbi_table=True):
    cfg = S.load(config or os.path.join("config", "sources.json"))
    sources = dict(cfg) if cfg else {}
    if live and email:
        os.environ.setdefault("NCBI_EMAIL", email)
    print("=" * 88)
    print("属级类型材料流水线：{}".format(genus))
    print("  可用源：{}".format(", ".join(sorted(sources)) or "（无外部源，仅联网检索）"))
    # 逐源体检：把"空载"显式打出来，别让用户以为核验做了
    for nm in sorted(sources):
        obj = sources[nm]
        if hasattr(obj, "n_rows"):
            n = 0
            try:
                n = obj.n_rows()
            except Exception:
                pass
            print("  [源] {:<14}{:<22}行数={:<8}{}".format(
                nm, type(obj).__name__, n,
                "⚠️ 空载 → 该源不参与判定，相关列将标「无公开数据库」" if n == 0 else "✅"))
        else:
            print("  [源] {:<14}{:<22}{}".format(nm, type(obj).__name__, "✅（联网）"))
    print("=" * 88)

    # ---- NCBI 结果表：把联网结果按属存表，重跑零下载 ----
    # 放在这里（acquire_species 之前），是为了让取种名那步的 Taxonomy 调用也进表。
    memo = None
    if ncbi_table:
        try:
            from . import ncbitable as _NT
            _root = os.path.dirname(os.path.dirname(os.path.dirname(
                os.path.abspath(__file__))))
            if ncbi_table is True:
                _tpath = os.path.join(_root, "data", "ncbi_{}.tsv".format(genus))
            else:
                _tpath = ncbi_table
            memo = _NT.MemoTable(_tpath, genus=genus)
            _n = S.attach_memo(sources, memo)
            _show = os.path.relpath(_tpath, _root) if os.path.isabs(_tpath) else _tpath
            print("  [ncbi表] {}：已有 {} 条 / 挂载 {} 个联网源".format(
                _show, len(memo.data), _n))
            if memo.load_errors:
                print("           ⚠️ 载入时有 {} 行读不通（已跳过，不影响其余）".format(
                    len(memo.load_errors)))
        except Exception as _e:
            # 表不可用绝不能拖垮主流程 —— 退回纯缓存模式即可
            print("  [ncbi表] ⚠️ 不可用（{}: {}）→ 退回纯缓存模式".format(
                type(_e).__name__, _e))
            memo = None

    species = acquire_species(genus, sources, corpus_dir, framework, limit, email, species_list)
    if not species:
        print("❌ 取不到物种清单；请用 --framework 投喂框架论文，或 --corpus 给本地语料")
        return 1

    fw = framework_columns(framework, genus) if framework else None
    rows = []
    total = len(species)
    _done = [0]
    _tbl = [None]

    def _flush(finished):
        """把【当前已完成的种】立刻落盘。

        为什么要每完成一个种就写一次：原来 write_tables 在循环【之后】才调用，
        于是点「停止」抛 KeyboardInterrupt 后，三列表**永远写不出来** ——
        抓到的东西没白费（缓存与 NCBI 表都已落盘），但拿不到半成品。
        实测确认过这个顺序问题。
        """
        try:
            _tbl[0] = write_tables(rows, out_dir, genus, fw, quiet=not finished)[0]
        except Exception as _e:
            print("  [report] ⚠️ 增量落盘失败（{}: {}）".format(type(_e).__name__, _e))
            return
        write_progress(out_dir, genus, _done[0], total, finished, _tbl[0])

    try:
        for i, sp in enumerate(species, 1):
            print("   [{:>3}/{}] {}".format(i, total, sp))
            rows.append(collect(sp, genus, sources))
            _done[0] = i
            _flush(False)          # ← 每完成一个种就落一次盘
    finally:
        # ⚠️ 必须用 finally：用户点「停止」时抛的是 KeyboardInterrupt。
        if memo is not None:
            memo.save(force=True)
            _st = memo.stats()
            print("  [ncbi表] 命中 {hits} / 联网新增 {misses} / 跳过 {skipped} → 共 {entries} 条".format(**_st))
            print("           {}".format(memo.path))
        _flush(_done[0] >= total and total > 0)      # 中断时把半成品标成「未完成」

    if _done[0] >= total and total > 0:
        print("\n完成。三列语义：种名 ｜ 菌株号（序列库侧，含问题）｜ 模式菌（类型材料侧，含问题）")
    else:
        print("\n⚠️ 这是【部分结果】：{}/{} 个种。已在表同目录标注 {}_进度.json".format(
            _done[0], total, genus.lower()))
    return 0


def annotate_links(src_dir, dst_dir=None, genus="Colletotrichum", sources=None,
                   icmp_budget=200, log=print):
    """【连接反哺】把连通组里已核实的结论借给「无目录」的号，另存一版三列表。

    为什么不直接改原表：反哺是**推断**，必须与原表并存、可对照。
    所以默认写到 <src_dir>_linked/，原表一字不动。

    返回统计 dict。
    """
    import shutil
    src_dir = src_dir if os.path.isabs(src_dir) else os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))), src_dir)
    dst_dir = dst_dir or (src_dir.rstrip("\\/") + "_linked")
    tsv = os.path.join(src_dir, "{}_三列表.tsv".format(genus.lower()))
    if not os.path.exists(tsv):
        return {"ok": False, "msg": "没有三列表：{}".format(tsv)}
    os.makedirs(dst_dir, exist_ok=True)
    if sources is None:
        sources = S.load(os.path.join(os.path.dirname(os.path.dirname(
            os.path.abspath(__file__))), "config", "sources.json"))

    # 1) 从表里抽号 + 现有判定
    lines = io.open(tsv, encoding="utf-8-sig").read().splitlines()
    nums, judge = [], {}
    RE_CELL = re.compile(r"^(.*?)\s*\[(.*)\]\s*$")
    for line in lines[1:]:
        f = line.split("\t")
        for ci in (1, 2):
            if len(f) <= ci or f[ci] in ("—", ""):
                continue
            for cell in f[ci].split(" ｜ "):
                m = RE_CELL.match(cell)
                if not m:
                    continue
                disp = m.group(1).strip()
                if disp.startswith("(strain)"):
                    continue
                fl = [x.strip() for x in m.group(2).split(";") if x.strip() and x.strip() != "—"]
                k = norm_bare(*re.match(r"^([A-Za-z][A-Za-z\.\-]{0,12})\s+([\w\.\-]+)$",
                                        disp).groups()) if re.match(
                    r"^([A-Za-z][A-Za-z\.\-]{0,12})\s+([\w\.\-]+)$", disp) else None
                if not k:
                    continue
                nums.append(disp)
                judge.setdefault(k, {"display": disp, "flags": fl})

    # 2) 分组（两座桥）
    log("[反哺] 号 {} 个，开始连通分组…".format(len(nums)))
    res = S.link_numbers(nums, sources, icmp_budget=icmp_budget)
    log("[反哺] 连通组 {} 组 ｜ 连接 {} 条".format(len(res["groups"]), len(res["edges"])))

    # 3) 推出可借用项
    inh = S.inherit_from_groups(res["groups"], judge)
    log("[反哺] 可借用的「无公开数据库」号：{} 个".format(len(inh)))

    # 4) 写出新表（逐格替换标记，其余不动）
    changed = 0
    out_lines = [lines[0]]
    for line in lines[1:]:
        f = line.split("\t")
        for ci in (1, 2):
            if len(f) <= ci or f[ci] in ("—", ""):
                continue
            cells = []
            for cell in f[ci].split(" ｜ "):
                m = RE_CELL.match(cell)
                if m:
                    disp = m.group(1).strip()
                    mm = re.match(r"^([A-Za-z][A-Za-z\.\-]{0,12})\s+([\w\.\-]+)$", disp)
                    k = norm_bare(*mm.groups()) if mm else None
                    if k and k in inh and "无目录" in m.group(2):
                        newf = [x for x in m.group(2).split(";") if x.strip() != "无目录"]
                        newf += inh[k]["flags"]
                        cells.append("{} [{}]".format(disp, "; ".join(newf)))
                        changed += 1
                        continue
                cells.append(cell)
            f[ci] = " ｜ ".join(cells)
        out_lines.append("\t".join(f))
    dst = os.path.join(dst_dir, "{}_三列表.tsv".format(genus.lower()))
    io.open(dst, "w", encoding="utf-8-sig", newline="\n").write("\n".join(out_lines) + "\n")
    # 组表也拷一份过去
    gt = os.path.join(src_dir, "号连接组.tsv")
    if os.path.exists(gt):
        shutil.copyfile(gt, os.path.join(dst_dir, "号连接组.tsv"))
    log("[反哺] 改了 {} 个单元格 → {}".format(changed, dst))
    return {"ok": True, "groups": len(res["groups"]), "edges": len(res["edges"]),
            "inherited": len(inh), "cells_changed": changed, "dst": dst,
            "stats": res.get("stats")}


def main(argv=None):
    ap = argparse.ArgumentParser(description="属级类型材料流水线：给属名 → 出三列表")
    ap.add_argument("--genus", required=True)
    ap.add_argument("--out", default="out")
    ap.add_argument("--config", default=None, help="源配置 json（默认 config/sources.json）")
    ap.add_argument("--framework", default=None, help="框架论文（PDF 抽出的 txt）→ 产出对照表")
    ap.add_argument("--corpus", default=None, help="本地语料目录（离线取种名）")
    ap.add_argument("--species-list", default=None, help="物种清单文件（每行一个种名，最可靠）")
    ap.add_argument("--limit", type=int, default=None, help="只处理前 N 个种（试跑用）")
    ap.add_argument("--email", default=os.environ.get("NCBI_EMAIL"))
    ap.add_argument("--ncbi-table", default="auto",
                    help="NCBI 结果表：auto=默认（data/ncbi_<属名>.tsv；命中零下载、"
                         "未命中联网并自动补写）；off=关闭；也可给一个具体路径")
    a = ap.parse_args(argv)
    _nt = a.ncbi_table
    if str(_nt).lower() in ("off", "no", "none", "0", "false"):
        _nt = False
    return run(a.genus, a.out, a.config, a.corpus, a.framework, a.limit, a.email,
               species_list=a.species_list, ncbi_table=_nt)


if __name__ == "__main__":
    sys.exit(main())
