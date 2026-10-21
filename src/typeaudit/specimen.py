# -*- coding: utf-8 -*-
"""标本型：抽「标本馆代码 + 编号」，并与 GenBank `/specimen_voucher` 对账。

⚠️ 这是一个**已实测的负结果**（78 属 / 166 条目）：
    候选种 12 个 → GenBank 有记录的仅 4 个 → 模式号被独立印证 **0 个**。
  三层根因（都不是"改正则"能解决的）：
    一、多数标本型新种在 GenBank 尚无序列；
    二、测序标本 ≠ 模式标本（Xylaria phyllocharis 模式 BPI 882429，
        GenBank 里的 voucher 是别的采集）；
    三、`/specimen_voucher` **无受控词表**——同一批里出现 `HAST 147515` /
        `528 (HAST, JF)` / `LIP:Bon 90132` / `UNIPAR:MR00378` 等写法。
  ⇒ 结论：标本型覆盖需要「标本馆号正规化层」，或改用 GBIF／馆藏目录作第三方源。
    本模块保留为可复现的负结果，**不建议**直接用于生产。

判据两轮校准（避免重犯）：
  宽松版：块级召回 98%，但 228 块抽出 2677 个候选（≈12/块），混入 TABLE:1 / ITS:19 / CYA:22
  硬白名单版：块级召回掉到 32%——白名单两头不讨好（真馆不在表里、作者姓氏混进来）
  分层版（当前）：T1 已知馆 / T2 括号+编号 / T3 弱，并剔除 GenBank 登录号前缀
"""
import io
import os
import re

from .textutil import norm

KNOWN_HERB = set("""HKAS HAST HMAS HMJAU K M NY BPI C G S W PC P UPS ZT BRIP IMI PDD TNS DAOM CUP
FH FLAS MICH MO MSC OSC PH PRM UC WIS LE LPS MEL TUR CM E F L LD O UME UPS CFMR DSM KRAM HAL
WA NIFS HLA GDGM HALL TENN VPI NYS CORT GENT UNIPAR MBT AMB LIP MFLU MFLUCC GJA GZU HGAS
MHHNU KUN PREM""".split())
ACC_PREFIX = set("""MW MT MN MK MH MG KY KJ KM JX JF GU EU AY DQ EF FJ HM HQ JN JQ KC KF KP KR
KT KX NR NG OM OP OR PP PQ PV PW PX PY PZ AB AJ AM AP AV""".split())
NOT_HERB = set("""TABLE TAB FIG FIGS ITS LSU SSU CYA MEA PDA OA SNA DNA RNA PCR NCBI BLAST HMM
OTU BOLD MLST RAPD AFLP ISSR SSR SNP CDS ORF ATP NADH COX RPB TEF TUB ACT GAPDH CHS SOD
SANGER ILLUMINA MISEQ HISEQ ACC NO COL STRAIN ISOLATE EX HOLO TYPE SPECIES NOV""".split())

RE_TYPE = re.compile(r"\b(Holotype|Holotypus|Isotype|Isotypus|Epitype|Neotype)\b\s*[:：]?\s*(.{0,400})", re.I)
RE_PAR = re.compile(r"\(([A-Za-z]{2,8})\s*[:,#]?\s*([0-9][0-9\-]{0,9})?\)")


def type_sentence(block):
    """取该块里的类型指定句（截断 420 字符）。"""
    m = RE_TYPE.search(block or "")
    return m.group(0)[:420] if m else None


def candidates(sentence):
    """分层候选 [(code, num, tier)]。T1 = 已知标本馆且编号 ≥2 位；其余丢弃。"""
    out, seen = [], set()
    for m in RE_PAR.finditer(sentence or ""):
        code, num = m.group(1).upper(), (m.group(2) or "")
        if code in ACC_PREFIX or code in NOT_HERB or (code, num) in seen:
            continue
        seen.add((code, num))
        if code in KNOWN_HERB and len(num) >= 2:
            out.append((code, num, "T1"))
    return out


def genbank_vouchers(path):
    """读一份 GenBank 平文件，返回 (记录数, voucher 列表)。"""
    if not path or not os.path.exists(path):
        return 0, []
    raw = io.open(path, encoding="utf-8", errors="replace").read()
    n = len(re.findall(r"^LOCUS", raw, re.M))
    vs = list(dict.fromkeys(re.findall(r'/specimen_voucher\s*=\s*"([^"]+)"', raw)))
    return n, vs


def reconcile(cands, vouchers):
    """把候选与 GenBank voucher 对账，返回命中列表 [(code, num, tier, voucher)]。"""
    hits = []
    for code, num, tier in cands:
        tgt = norm(code + num)
        for v in vouchers:
            if tgt and (tgt in norm(v) or norm(v) in tgt):
                hits.append((code, num, tier, v))
                break
    return hits
