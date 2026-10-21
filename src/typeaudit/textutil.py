# -*- coding: utf-8 -*-
"""文本工具：去标签、规范化、按「种名 + 作者」切物种处理块。

⚠️ 切块规则的血泪教训（实测 166 条目 / 78 属）：
    不要要求作者名以「大写字母 + 点」开头（`F. Liu` 那种）。
    Diaporthe 那篇写的是 `Diaporthe afzeliae Monkai and S. Lumyong, sp. nov.`
    ——作者是**完整姓氏、不带点**，严格规则会**一块都切不出来**（0 块）。
    默认用变体 B：种加词后面接任意大写字母**或左括号**即可。
    第二个坑（冒烟测试抓到）：重组合名写作 `Xylaria phyllocharis (Mont.) Fr.`，
    种加词后接的是 `(` 而不是大写字母——只要求大写字母会整类漏掉。
"""
import re

DEFAULT_GENUS = "Colletotrichum"

# 会被误当种加词的功能词（变体 D 用）
STOPWORDS = set("""species sp spp isolate isolates strain strains and the with from were was are can may
for this that these those all some new nov gen et al fig figs table tables note notes""".split())


def strip_tags(s):
    """XML/HTML -> 纯文本（去标签、压缩空白）。输入已是纯文本时等价于压缩空白。"""
    return re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", s))


def norm(s):
    """两侧统一规范化后再比：去空白/冒号/下划线/连字符/点，转大写。"""
    return re.sub(r"[\s:_\-\.\(\)]", "", s).upper()


def heading_rx(genus, variant="B"):
    """切块正则。variant: A 大写+点（不推荐）/ B 后接大写（默认）/ C 任意提及（易过度切分）。"""
    g = re.escape(genus)
    if variant == "A":
        return re.compile(r"\b" + g + r"\s+([a-z][a-z\-]{3,})\s+(?=[A-Z][a-z]*\.)")
    if variant == "C":
        return re.compile(r"\b" + g + r"\s+([a-z][a-z\-]{3,})\b")
    # B（默认）：种加词后接 大写字母 / 左括号 / `sp.` / `nom.`。
    # 三个分支各对应一个实测缺口：
    #   `(?=[A-Z(])` 的括号分支：重组合名 `Xylaria phyllocharis (Mont.) Fr.`
    #     ——括号里是基名作者，只要求大写字母会整类漏掉（冒烟测试抓到）。
    #   `sp\.` 分支：`Colletotrichum tongrenense sp. nov.` 这种**不写作者**的标题
    #     ——`sp.` 是小写，只要求大写字母会把整个块往后挪，造成块内丢掉模式指定段
    #     （E7 就是因此抽不到链）。
    #   `nom\.` 分支：`… nom. nov.` / `… nom. inval.` 同理。
    return re.compile(r"\b" + g + r"\s+([a-z][a-z\-]{3,})\s+(?=[A-Z(]|sp\.|nom\.)")


def blocks(text, genus=DEFAULT_GENUS, variant="B", stop=None):
    """把全文切成 {种加词: 该种的块文本}。同名保留最长块。

    ⚠️ 不要用固定字符窗口（如「种名后 N 字符」）——实测会越过种界，
       把相邻种的标本号串进来（v1 的 bug）。
    """
    rx = heading_rx(genus, variant)
    hits = [m for m in rx.finditer(text) if not (stop and m.group(1).lower() in stop)]
    out = {}
    for i, m in enumerate(hits):
        end = hits[i + 1].start() if i + 1 < len(hits) else len(text)
        ep, blk = m.group(1).lower(), text[m.start():end]
        if ep not in out or len(blk) > len(out[ep]):
            out[ep] = blk
    return out

def paras(raw):
    """按 XML 的 <p> 元素切段（保留顺序），返回 [纯文本段落]。

    ⚠️ **压平文本上的正则切块已到极限**——本项目实测：
      · 在 `.strip_tags()` 压平后的文本上用「种名 + 作者」正则切块，
        反复出现**块污染**（块跨表格/图注，同一号落进多个物种块）
        与**取块取错**（取到该种的首次提及，而 Typus 段不在那一段）。
      · 改用 <p> 段落结构后，`E9`（图注说 ex-type 是 NN058985、正文 Typus 说 NN071035）
        立刻可判：图注独有该号。而压平文本路线上怎么调正则都不稳。
    结论：**只要语料是结构化全文（PMC XML / JATS），就该按段落走，不要在压平文本上做正则切块。**
    边界：仍有材料不在 <p> 里（表格、列表、补充材料）——`E7`（一号两主）即因此未检出。
    """
    import re as _re
    return [strip_tags(p) for p in _re.findall(r"<p[^>]*>(.*?)</p>", raw, _re.S)]

# ==================== 文档序单元（结构化语料的首选切分路径） ====================
# ⚠️ 为什么不继续在压平文本上切块：实测反复出现两类病——
#    ① 块污染：块跨表格/图注，同一个号落进多个物种块；
#    ② 取块取错：取到该种的首次提及，而模式指定段不在那一段。
#    按 XML 的单元走可以从根上避免。文档序是关键：caption 必须晚于它所属的正文段。
# ⚠️ `<title>` 必须算单元：JATS 的**种级标题**就在它里面
#     （`<sec><title>Colletotrichum jishouense</title><p>…</p></sec>`）。
#     漏掉它 → 标题从未被看到 → 模式指定段被归到上一个种名下（E7 的归属错就出在这）。
#     另注：原文常把属名标成 `<italic>Colletotrichum</italic> jishouense`，
#     所以**在原始 XML 里直接搜字面串 "Colletotrichum jishouense" 是搜不到的**——
#     必须先按单元去标签再判断。
UNIT_PAT = re.compile(
    r"<title[^>]*>(.*?)</title>"
    r"|<p[^>]*>(.*?)</p>"
    r"|<caption[^>]*>(.*?)</caption>"
    r"|<td[^>]*>(.*?)</td>"
    r"|<th[^>]*>(.*?)</th>"                  # 表头单元格（原缺）
    r"|<list-item[^>]*>(.*?)</list-item>"
    r"|<def[^>]*>(.*?)</def>"                # 术语定义项
    r"|<fn[^>]*>(.*?)</fn>"                  # 脚注（模式材料有时写在脚注）
    r"|<boxed-text[^>]*>(.*?)</boxed-text>"
    r"|<disp-quote[^>]*>(.*?)</disp-quote>"
    r"|<supplementary-material[^>]*>(.*?)</supplementary-material>", re.S)
UNIT_KINDS = ["title", "p", "caption", "table_cell", "table_head", "list_item",
              "def", "footnote", "boxed_text", "disp_quote", "supplement"]

# 读全率自检用的关键词：这些内容若落在未覆盖区，必须报警
_CHECK_KEY = re.compile(r"(Holotype|holotype|Holotypus|ex-type|Typus|Isotype|Epitype|type material)", re.I)


# ==================== 精确的 span 版单元与读全率 ====================
# ⚠️ 教训：读全率必须用**精确 span** 判定。
#    我第一版用「关键词的上下文串是否整体出现在覆盖文本里」近似，
#    跨单元边界就误报，把 6 处虚报成 287 处。
def units_spans(raw):
    """返回 [(kind, text, start, end)]，start/end 为在 raw 中的字符偏移。"""
    out = []
    for m in UNIT_PAT.finditer(raw):
        for i, k in enumerate(UNIT_KINDS, start=1):
            g = m.group(i)
            if g is not None:
                t = strip_tags(g)
                if t.strip():
                    # 内层文本在 raw 中的实际起点：跳过开标签
                    inner = m.start(i)
                    out.append((k, t, inner, inner + len(g)))
                break
    return out


def coverage_report(raw):
    """**读全率自检**（精确 span 版）。

    报告：全文去标签字符数 / 被单元覆盖的字符数 / 覆盖率 /
          以及**落在未覆盖区的类型指定关键词**（这才是要紧的漏）。
    """
    all_text = strip_tags(raw)
    # ⚠️ 坐标系必须一致：单元文本长度 vs 全文去标签长度（**都去标签**）。
    #    我第一版拿「raw 里的 span 长度」比「去标签后的字符数」，算出 231% 的荒谬覆盖率。
    full = units_full(raw)
    covered_len = sum(len(t) for _, t, _, _ in full)
    # ⚠️ 两处口径必须一致：长度用 units_full，span 判定也必须用 units_full。
    #    我上一版长度用 units_full、判定用 units_spans，导致残差里的关键词永远显示"未覆盖"。
    sorted_spans = sorted((a, b) for _, _, a, b in full)
    import bisect
    starts = [a for a, _ in sorted_spans]
    def is_covered(pos):
        i = bisect.bisect_right(starts, pos) - 1
        return i >= 0 and sorted_spans[i][0] <= pos < sorted_spans[i][1]
    # 关键词必须用**去标签后**的文本定位；为此在 raw 内按相同顺序找
    # 简化而精确的做法：逐个 raw 命中判断
    unc = []
    for m in re.finditer(r"(Holotype|holotype|Holotypus|ex-type|Typus|Isotype|Epitype|type material)", raw):
        if not is_covered(m.start()):
            unc.append("[{}] {}".format(m.group(1), strip_tags(raw[max(0, m.start() - 70): m.start() + 50])[:72]))
    n_all = len(re.findall(r"(Holotype|holotype|Holotypus|ex-type|Typus|Isotype|Epitype|type material)", raw))
    return {"chars_all": len(all_text), "chars_covered": covered_len,
            "ratio": round(covered_len / len(all_text), 4) if all_text else 0,
            "kw_total": n_all, "kw_uncovered": len(unc), "uncovered_samples": unc[:6]}


def uncovered_mass_by_element(raw):
    """未覆盖的**字数**落在哪些元素里——用于判断漏得要不要紧。"""
    from collections import Counter
    TAG = re.compile(r"<(/?)([a-zA-Z][\w\-]*)([^>]*)>")
    tags = [(m.start(), m.group(1), m.group(2).lower(), m.group(3)) for m in TAG.finditer(raw)]
    spans = sorted((a, b) for _, _, a, b in units_spans(raw))
    import bisect
    starts = [a for a, _ in spans]
    def covered(pos):
        i = bisect.bisect_right(starts, pos) - 1
        return i >= 0 and spans[i][0] <= pos < spans[i][1]
    mass = Counter()
    stack = []
    for pos, slash, name, attrs in tags:
        if slash:
            if stack and stack[-1] == name:
                stack.pop()
        elif not attrs.rstrip().endswith("/"):
            stack.append(name)
        if not covered(pos):
            mass[stack[-1] if stack else "?"] += len(name) + 2
    return mass

# 排除区：元数据（参考文献表 / 文章头）不算正文
RE_EXCLUDE = re.compile(r"<ref-list[^>]*>.*?</ref-list>|<front[^>]*>.*?</front>", re.S)


def units_full(raw):
    """`units()` + **残差单元** —— 收全正文，包括裸写在 `<sec>` 下的文本。

    为什么必须补（实测）：JATS 里有裸文本直接写在 `<sec>` 下，例如
        `<sec><title>…</title> Holotype. MFLU 17-1201 …</sec>`
    这类内容不被 `UNIT_PAT` 命中，**静默漏掉 6/618 个类型指定关键词**
    （占 1%，但每一条都可能是一个模式凭证）。

    做法：把已识别单元占用的区间挖掉，剩下的空隙文本（排除参考文献表与文章头）
    作为 kind="residual" 的单元补回，并保持**文档序**。
    """
    spans = sorted((a, b, t) for _, t, a, b in units_spans(raw))
    excl = [(m.start(), m.end()) for m in RE_EXCLUDE.finditer(raw)]

    def in_excl(i, j):
        return any(a < j and i < b for a, b in excl)

    gaps = []
    cur = 0
    for a, b, _ in spans:
        if a > cur:
            gaps.append((cur, a))
        cur = max(cur, b)
    if cur < len(raw):
        gaps.append((cur, len(raw)))

    resid = []
    for a, b in gaps:
        if in_excl(a, b):
            continue
        t = strip_tags(raw[a:b])
        if len(t.strip()) >= 8:
            resid.append((a, b, "residual"))
    # 按位置排序输出
    merged = [(a, b, "unit") for a, b, _ in spans] + [(a, b, k) for a, b, k in resid]
    merged.sort()
    res = []
    for a, b, k in merged:
        t = strip_tags(raw[a:b])
        if t.strip():
            res.append((k if k == "residual" else "unit", t, a, b))
    return res


def units_all(raw):
    """正文全部文本（已识别单元 + 残差单元），按文档序。"""
    return units_full(raw)

def segment(raw, genus=DEFAULT_GENUS, is_xml=None):
    """**统一切分入口**：结构化 XML 走单元路径，纯文本走 blocks()。

    返回 {种加词: 文本}。两者语义一致，调用方无需关心底层。

    ⚠️ 本函数曾被误删一次（我删除旧版 coverage_report 时把它一起切掉了，
       gate.py 立刻 ImportError）。凡用「从 A 切到 B」的方式删代码，
       务必确认中间没有别的定义——或用 AST 级工具，别用字符串切片。
    """
    if is_xml is None:
        is_xml = bool(re.search(r"<p[^>]*>|<title[^>]*>", raw or ""))
    if is_xml:
        return {ep: " ".join(ts) for ep, ts in species_units(raw, genus).items()}
    return blocks(strip_tags(raw or ""), genus, "B")

def is_heading(text, genus=DEFAULT_GENUS, max_len=400):
    """判断一个单元是否为「种级标题」。返回种加词或 None。

    判据与 heading_rx 一致（种加词后接 大写/括号/sp./nom.），另加长度上限——
    正文里顺带提及不会短到 400 字符以内又恰好以属名开头。
    """
    if len(text) > max_len:
        return None
    m = re.match(r"\s*" + re.escape(genus) + r"\s+([a-z][a-z\-]{3,})\s*"
                 r"(?=[A-Z(]|sp\.|nom\.|,|$)", text)
    return m.group(1).lower() if m else None


def species_units(raw, genus=DEFAULT_GENUS):
    """按单元归属物种：{种加词: [该种名下的全部单元文本]}。

    归属规则：遇到种级标题即切换「当前物种」；其后单元依次归入，直到下一个种级标题。
    表格单元、列表项、残差文本也按同样规则归属——这正是 E7（一号两主，
    模式材料在表格/列表里）能取到的原因。

    ⚠️ `<title>` 必须算单元：JATS 的种级标题就在它里面
       （`<sec><title>Colletotrichum jishouense</title>…`）。漏掉它会让标题看不到，
       模式指定段被归到上一个种名下（E7 归属错就出在这）。
    """
    out, cur = {}, None
    for kind, txt in units(raw):
        ep = is_heading(txt, genus)
        if ep:
            cur = ep
            out.setdefault(cur, [])
        if cur is not None:
            out[cur].append(txt)
    return out

def units(raw):
    """按**文档顺序**产出 (kind, text)。kind ∈ UNIT_KINDS。

    嵌套已被自然处理：`<td><p>x</p></td>` 会作为 table_cell 整体命中，
    其内部的 <p> 不会被重复匹配（finditer 从匹配结束处继续）。
    """
    out = []
    for m in UNIT_PAT.finditer(raw or ""):
        for i, k in enumerate(UNIT_KINDS, start=1):
            g = m.group(i)
            if g is not None:
                t = strip_tags(g)
                if t.strip():
                    out.append((k, t))
                break
    return out
