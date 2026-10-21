# -*- coding: utf-8 -*-
"""培养物型：从分类学全文中抽「模式指定段」与等价保藏号链。

两条实测教训：
  A7 标签变体：SiM 2022 把类型指定写成拉丁 `Typus:`，不是 `Type:`。
      只搜 `Type:` 得 0 处 ⇒ 误判「该文无类型指定」= **假阴性**。
      ⇒ 必须扫标签变体表，并打印**每标签命中数**，让 0 命中可见。
  C10 等价链：`ex-type culture CFCC 54250 = ACCC 39709 = CNUCC 307307`
      三号同株。解析出这条链，才能把「论文图标签的号」与「GenBank 的号」连上。
      优先级：**正文有 Typus 段时，图标签让位于正文**。

适用边界（实测 166 条目 / 78 属）：本方法**只适用于培养物型分类学**。
  培养物型抽出率 63%；标本型（正文只写 `Holotype`）抽出率 **0%**——硬边界，改正则救不了。
  标本型请用 `specimen` 模块。
"""
import re

LABEL_VARIANTS = [
    "Typus", "Type", "Holotype", "holotype", "Epitype", "Isotype", "Neotype", "Paratype",
    "ex-type culture", "ex-type", "ex-holotype", "Type material", "Type status",
]

RE_CHAIN = re.compile(r"(?:ex-type culture|ex-type|ex-holotype)\s+([^);]{5,200})", re.I)
RE_TOKEN = re.compile(r"\b([A-Z][A-Za-z]{1,11})\s*[:：]?\s*(\d+[\w.\-]*)\b")
RE_HOLO = re.compile(r"\bholotype\s+([A-Z]{2,8}\s*\d+[\w.\-]*)", re.I)


def label_census(text, labels=LABEL_VARIANTS):
    """A7 自检：每个标签变体的命中数。0 命中的标签要显式提示，别直接下结论。"""
    return {lb: len(re.findall(re.escape(lb), text, re.I)) for lb in labels}


def parse_chains(fragment_text):
    """抽该文本内所有 ex-type 记录 [(原文片段, [号, ...])]。链长 1 = 只给了一个号。"""
    out = []
    for m in RE_CHAIN.finditer(fragment_text):
        toks = ["{}{}".format(a, b) for a, b in RE_TOKEN.findall(m.group(1))]
        if toks:
            out.append((m.group(1).strip(), toks))
    return out


def designation_for(blk):
    """在**该种自己的块内**找等价链与 holotype。返回 (链, holotype, 块)。"""
    if blk is None:
        return None, None, None
    ch = parse_chains(blk)
    ho = RE_HOLO.search(blk)
    return (ch[0][1] if ch else None), (ho.group(1).strip() if ho else None), blk
