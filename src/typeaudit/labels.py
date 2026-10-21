# -*- coding: utf-8 -*-
"""图标签解析：从**PDF 文本层**抽分类学修订图里的 `_T` 模式株标签。

为什么需要单独一个模块（关键实测）：
  分类学修订论文的模式株号，很多只出现在**系统发育树的图**上（标签带 `_T` 后缀）。
  而这类图在 PMC XML 里是**图片**——实测 SiM 2022：
      PMC XML（320,979 字符）：`_T` 出现 **0** 次
      PDF 文本层（178,230 字符）：`_T` 出现 **299** 次
  ⇒ 只用 PMC XML 会**完全拿不到**这一侧的号（副表3 的 `SiM2022_figure_T` 列全靠它）。

格式样例（PDF 文本层抽出的实际写法）：
    C. falcatum CBS 147945_T        C. axonopodi IMI 279189_T
    C. guizhouensis CGMCC 3.15112_T C. syzygicola MFLUCC 10-0624_T
"""
import re

# 图标签：属名缩写 + 种加词 + 号 + `_T`
RE_T_LABEL = re.compile(r"\b[A-Z]\.\s*([a-z][a-z\-]{2,})\s+([A-Z][A-Za-z]{1,11})\s*[:#]?\s*(\d[\w.\-]*)_T\b")

# 拼写变体表（规范 L-04）。实测：不建这张表，图标签与对账表会因拼写差异对不上——
# 副表3 里 4 个"看似未命中"的种全部是这一现象（axonopodis/axonopodi 等）。
SPELLING_VARIANTS = {
    "axonopodis": "axonopodi", "duyunense": "duyunensis", "laticiphilum": "laticiphila",
    "tabaci": "tabacum", "syzygiicola": "syzygicola", "guizhouense": "guizhouensis",
}


def canon(epithet):
    """把拼写变体归一到文献里实际使用的形（L-04）。"""
    e = (epithet or "").strip().lower()
    return SPELLING_VARIANTS.get(e, e)


def extract_T_labels(text):
    """从文本（PDF 层或任意纯文本）抽出 {种加词(归一后): set(号)}。"""
    out = {}
    for m in RE_T_LABEL.finditer(text or ""):
        ep = canon(m.group(1))
        out.setdefault(ep, set()).add((m.group(2) + m.group(3)).replace(" ", ""))
    return out


def labels_by_species(text):
    """返回 {种加词: '号'}（同种多号时用 ' | ' 连接，保持确定性顺序）。"""
    lab = extract_T_labels(text)
    return {ep: " | ".join(sorted(v)) for ep, v in sorted(lab.items())}
