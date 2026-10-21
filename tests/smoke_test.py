# -*- coding: utf-8 -*-
"""冒烟测试：不联网、不需要任何数据文件，验证核心解析逻辑。

跑法：python tests/smoke_test.py     （在仓库根目录）
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))

from typeaudit import designations, gate, specimen, textutil   # noqa: E402

FAIL = []


def check(name, got, want):
    ok = got == want
    print("   [{}] {:<44}{}".format("OK  " if ok else "FAIL", name, "" if ok else "得到 {} 应为 {}".format(got, want)))
    if not ok:
        FAIL.append(name)


# ---- 合成语料：一个培养物型种（带等价链）+ 一个标本型种 ----
CULTURE_TEXT = (
    "Colletotrichum bambusicola C.L. Hou & Q.T. Wang, Mycologia 113: 452. 2021. "
    "Typus: China, Guangxi, on seeds of Phyllostachys edulis, Sep. 2016 (holotype CAF80001, "
    "ex-type culture CFCC 54250 = ACCC 39709 = CNUCC 307307). "
    "Colletotrichum afzeliae Monkai and S. Lumyong, sp. nov. "
    "Typus: Thailand, on leaves, 2018 (holotype MFLU 18-0001). "
)
# 前言段：属名出现次数须 ≥ gate.MIN_GENUS_MENTIONS，否则会被"样本真伪"轴先拦下。
# 手工数容易数错（本文件第一版就数错了），故**程序化满足**：按阈值重复足够多遍。
_PRE = " ".join(["Xylaria is a wood-inhabiting genus and Xylaria species were examined."]
                * (gate.MIN_GENUS_MENTIONS + 2))
SPECIMEN_TEXT = (
    _PRE + " "
    "Xylaria allima F. Liu & L. Cai, sp. nov. "
    "Holotype: China, Hainan, on wood, 2016, X. Liu 1234 (HAST 145965). "
    "Xylaria phyllocharis (Mont.) Fr. "
    "Holotype: Cuba, on wood (BPI 882429). "
)

print("=" * 92)
print("【冒烟测试：typeaudit 核心解析】")
print("=" * 92)

# 1) A7 标签变体：Typus 必须被算到
cen = designations.label_census(CULTURE_TEXT)
check("标签变体表能命中 Typus（=2）", cen["Typus"], 2)

# 2) 切块：变体 A 会漏（作者是完整姓氏 Monkai），变体 B 不漏
blkA = textutil.blocks(CULTURE_TEXT, "Colletotrichum", "A")
blkB = textutil.blocks(CULTURE_TEXT, "Colletotrichum", "B")
check("变体 A 切块数（作者带点者通过）", len(blkA), 1)
check("变体 B 切块数（完整姓氏也通过）", len(blkB), 2)

# 3) 等价链解析
ch, ho, _ = designations.designation_for(blkB.get("bambusicola"))
check("bambusicola 等价链", ch, ["CFCC54250", "ACCC39709", "CNUCC307307"])
check("bambusicola holotype", ho, "CAF80001")

# 4) 只引证、无 Typus 段的种不应抽出链（对照：afzeliae 有 Typus 但无 ex-type）
ch2, ho2, _ = designations.designation_for(blkB.get("afzeliae"))
check("afzeliae 无 ex-type 链", ch2, None)

# 5) 标本型：Ha 标本馆代码 + 编号
blkS = textutil.blocks(SPECIMEN_TEXT, "Xylaria", "B")
s_allima = specimen.type_sentence(blkS.get("allima"))
check("标本型类型句命中 Holotype", s_allima is not None and "Holotype" in s_allima, True)
cands = specimen.candidates(s_allima)
check("抽出标本馆代码 + 编号", (cands[0][0], cands[0][1]) if cands else None, ("HAST", "145965"))
s_phyl = specimen.type_sentence(blkS.get("phyllocharis"))
c2 = specimen.candidates(s_phyl)
check("BLI 型短代码也能抽", (c2[0][0], c2[0][1]) if c2 else None, ("BPI", "882429"))

# 6) 对账：规范化后应命中（模拟 GenBank voucher 写法差异）
hits = specimen.reconcile([("HAST", "145965", "T1")], ["HAST 145965", "MFLU 18-0001"])
check("规范化解掉写入格式差异（'HAST 145965'）", len(hits), 1)
hits2 = specimen.reconcile([("BPI", "882429", "T1")], ["JF-GUY-12-073", "528 (HAST, JF)"])
check("非同一采集则不命中（负结果复现）", len(hits2), 0)

# 7) C1 门：合成语料的基质判定与三轴输出
check("测试语料满足属名阈值（防测试自身退化）",
      gate.genus_mentions(SPECIMEN_TEXT, "Xylaria") >= gate.MIN_GENUS_MENTIONS, True)
fg_c = gate.full_gate(CULTURE_TEXT, "Colletotrichum")
check("培养物型语料判为 culture", fg_c["basis"], "culture")
check("培养物型语料不出『标本型不可审』", fg_c["verdict"].startswith("不可审 · 等价链法"), False)
fg_s = gate.full_gate(SPECIMEN_TEXT, "Xylaria")
check("标本型语料判为 specimen", fg_s["basis"], "specimen")
check("标本型语料被判不可审（硬边界）", "等价链法不适用" in fg_s["verdict"], True)
fg_fake = gate.full_gate("Ophelimus bipolaris sp. nov. " * 5, "Bipolaris")
check("假样本（属名仅出现在种加词）被判无效", fg_fake["verdict"].startswith("样本无效"), True)


# 8) J 规则库：用规范里的真实示例做单元检查
from typeaudit import rules as R   # noqa: E402
cv = R.coverage()
check("J 规则库覆盖 21 条", cv["total"], 21)
check("离线可跑规则 8 条", len(cv["offline_ready"]), 8)
check("J-13 命中近失配（E1 的 127945/147945）", len(R.j13_near_miss(["CBS 127945"], ["CBS 147945"])), 1)
check("J-13 不误报无关号", len(R.j13_near_miss(["CBS 127945"], ["CBS 999999"])), 0)
check("J-07 抓一号两主（E7 的 GMBC0209）",
      len(R.j07_one_number_two_species({"GMBC0209": {"jishouense", "tongrenense"}})), 1)
check("J-18 显眼度 E1 = 5/6", R.j18_conspicuity(6, 5)["conspicuity"], 0.8333)
check("J-18 显眼度 E6 = 0/9", R.j18_conspicuity(9, 0)["conspicuity"], 0.0)
_j01 = R.j01_cross_source({"a": "CBS 1"}, {"a": "CBS 2"})
check("J-01 分歧判为 candidate", _j01[0]["severity"], "candidate")
_j01b = R.j01_cross_source({"a": "CBS 1"}, {})
check("J-01 单侧判为 untestable（不混进可疑）", _j01b[0]["severity"], "untestable")
_j08 = R.j08_batch_regularity([{"group": "G", "num": "20514"}, {"group": "G", "num": "20513"},
                               {"group": "G", "num": "20515"}, {"group": "G", "num": "20516"},
                               {"group": "G", "num": "20517"}, {"group": "G", "num": "2051"}],
                              group_key="group")
check("J-08 位数孤例命中", len(_j08), 1)
_j08b = R.j08_batch_regularity([{"group": "L", "num": "951"}, {"group": "L", "num": "13850"}],
                               group_key="group")
check("J-08 小批不误报（门槛生效）", len(_j08b), 0)


# 9) 文档序单元：<title> 必须算单元（否则种级标题看不到）
from typeaudit.textutil import segment, species_units, units   # noqa: E402
_XML = ("<article><body><sec><title>Colletotrichum alpha</title>"
        "<p>Typus: (holotype HMAS 1, ex-type culture CGMCC 3.20501 = LC13801).</p></sec>"
        "<sec><title>Colletotrichum beta</title>"
        "<p>Typus: (holotype HMAS 2, ex-type culture CGMCC 3.20502 = LC13802).</p></sec>"
        "</body></article>")
_u = units(_XML)
check("<title> 被算作单元", any(k == "title" for k, _ in _u), True)
_su = species_units(_XML, "Colletotrichum")
check("种级标题切出两个种", sorted(_su), ["alpha", "beta"])
check("alpha 的 Typus 归到 alpha", "CGMCC 3.20501" in " ".join(_su["alpha"]), True)
check("beta 的 Typus 未混进 alpha", "CGMCC 3.20502" not in " ".join(_su["alpha"]), True)
_sg = segment(_XML, "Colletotrichum")
check("segment() 返回同结构", sorted(_sg), ["alpha", "beta"])
check("纯文本走 fallback 不报错", isinstance(segment("Colletotrichum alpha F. Liu, sp. nov.", "Colletotrichum"), dict), True)


# 10) 新模块：labels（PDF 图标签）/ reconcile（判定 ①–⑥）/ layers（层一、层二）
from typeaudit import labels as LB, reconcile as RC, layers as LY   # noqa: E402
_TXT = "C. falcatum CBS 147945_T C. axonopodi IMI 279189_T C. guizhouensis CGMCC 3.15112_T"
_lab = LB.extract_T_labels(_TXT)
check("PDF 图标签抽出 3 个种", len(_lab), 3)
check("标签号去空格", sorted(_lab["falcatum"]), ["CBS147945"])
check("拼写归一 axonopodis→axonopodi", LB.canon("axonopodis"), "axonopodi")
check("拼写归一后能对上（副表3 写 axonopodis，PDF 写 axonopodi）", "axonopodi" in _lab, True)

check("rule ② 共有号 → 一致", RC.verdict("CBS 1", "CBS 1 | CBS:1")[0], "✅ 一致")
check("rule ⑥ 仅一侧 → 不判为不一致", RC.verdict("", "CBS 1")[0].startswith("仅"), True)
check("rule ③ 标本↔培养物", RC.verdict("HMAS 1", "CBS 1")[0], "🔶 标本/活培养物配对")
check("rule ④ 两侧均培养物", RC.verdict("CBS 1", "CGMCC 2")[0], "🔶 疑似同株双保藏")
check("rule ① 人工裁定最高优先（覆盖机械规则）",
      RC.verdict("CBS 1", "CBS 1", species="X", manual_hits={"X"})[0], "🔴 有确认错值")
check("连字符号可解析（CGMCC-3.18952）", RC.split_nums("CGMCC-3.18952")[0][0], "CGMCC")

check("层二·code_like 形态", LY.classify_value("CBS 147945")[0], "code_like")
check("层二·实验室内码归 other", LY.classify_value("24-4")[0], "other")
check("层二·无空格也可解析", LY.classify_value("GMBC0209")[1], "GMBC")

print()
print("   失败项 = {}".format(len(FAIL)))
sys.exit(1 if FAIL else 0)
