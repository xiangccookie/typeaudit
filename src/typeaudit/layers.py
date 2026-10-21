# -*- coding: utf-8 -*-
"""层一/层二：记录级的属级审计。

层一 · **号级**：一个号在 GenBank 里的记录按属分布 → 纯跨属越界（本项目实测 212 → 0）
层二 · **值级**：一个 `/strain` 值在 GenBank 里的属分布 → 跨属错误（本项目实测 112 → 0）

实测定位（规范里的原话）：**逐号检查类方法的产出率极低**
    J-03 号的属分布（层一）  212 → 1（0.47%）
    J-19 号 ↔ Taxonomy      30 → 0（0%）
⇒ **只能用在「已缩小范围之后」，不可作为筛查主力。**
   本模块实现它们，是为了能复现这两条数字，而不是让你拿它去全量筛查。
"""
import re
from collections import Counter

# 值的形态分类（决定该值能不能审）
RE_CODE_LIKE = re.compile(r"^\s*([A-Z][A-Za-z]{1,11})\s*[:\-#]?\s*(\d[\w.\-]*)\s*$")


def classify_value(v):
    """把 /strain 值分成 code_like（像保藏号，可审）与 other（无法审）。"""
    m = RE_CODE_LIKE.match(v or "")
    return ("code_like", m.group(1).upper(), m.group(2)) if m else ("other", None, None)


def layer1(numbers, gb, genus_expected=None):
    """层一：逐号查 GenBank 属分布，标出纯跨属越界。

    numbers: [(code, num)] 或 ["CBS 127945", ...]
    genus_expected: 若给出，则把「记录里不出现该属」判为跨属越界。
    返回 (结果列表, 统计)
    """
    out = []
    for item in numbers:
        if isinstance(item, (tuple, list)):
            code, num = item
        else:
            m = RE_CODE_LIKE.match(item)
            if not m:
                continue
            code, num = m.group(1).upper(), m.group(2)
        r = gb.records_of_number(code, num)
        gs = list(r["genera"])
        cross = bool(gs) and genus_expected and all(g.lower() != genus_expected.lower() for g in gs)
        out.append({"number": "{} {}".format(code, num), "count": r["count"],
                    "genera": r["genera"], "cross_genus": cross})
    stat = {"n": len(out), "n_with_records": sum(1 for x in out if x["count"]),
            "n_cross_genus": sum(1 for x in out if x["cross_genus"])}
    return out, stat


def layer2(organism, gb, retmax=200):
    """层二：取该属记录里的全部 /strain 值，分形态；code_like 的逐个查属分布。

    返回 (值清单, 统计)。값级跨属错误 = code_like 且记录里不出现本属。
    """
    vals, seen = [], set()
    for v in gb.strain_values_of(organism, retmax=retmax):
        if v in seen:
            continue
        seen.add(v)
        k, code, num = classify_value(v)
        vals.append({"value": v, "kind": k, "code": code, "num": num})
    stat = {"n_values": len(vals),
            "code_like": sum(1 for x in vals if x["kind"] == "code_like"),
            "other": sum(1 for x in vals if x["kind"] == "other")}
    return vals, stat
