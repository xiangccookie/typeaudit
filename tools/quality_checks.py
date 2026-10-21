#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
quality_checks.py —— 五条聚合检查（对应 data/known_errors.tsv 里的 F14–F18）。

【诚实边界（写死在这里，报告里必须原样带出去）】
    本工具核对的是【两列保藏号本身有没有问题】：
        第 2 列 菌株号 = 序列库随手记的号（GenBank /culture_collection + /strain）
        第 3 列 模式菌 = 被声明为类型的号（NCBI type material + GenBank /type_material）
    它能答的是这几问：
        · 这个号在机构目录里存在吗？        → 号不存在
        · 机构代码写错了吗？                → 未知机构代码 / 疑似非机构编号
        · 目录记的名与声明名对得上吗？      → 目录记为别种 / 基名属 / 异名（非错）
        · 这号在目录里有模式地位标注吗？    → 未核（没有标注 ≠ 不是模式菌株）
    它【不能】答的是这几问（缺的是另一条轴，不是实现 bug）：
        · 这个号的序列真的属于它声称的种吗？        → 需要序列轴
        · 这株菌到底是不是模式菌株？                → 需要原始文献/模式名录
        · 论文的声明对不对？                        → 需要文献轴
        · 号对但菌株错（标本装错）                  → 需要序列轴 + 实物核对

    所以：**本工具是高精度候选生成器，不是穷尽查错器。**
    空表 ≠ 没有错。

用法
    python tools/quality_checks.py                          # 用 out/ 里已拉到的 CBS 表
    python tools/quality_checks.py --cbs-tsv <路径> --genus Colletotrichum
产物
    out/quality_report.json   每条检查的命中数与样例（供 check_known_errors.py 消费）
"""
import sys
try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

import argparse, ast, io, json, os, re
from collections import Counter, defaultdict

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "src"))

SCOPE_INSIDE = [
    "号在该机构的公开记录里存在吗（号不存在）",
    "机构代码写错了吗（未知机构代码 / 疑似非机构编号）",
    "机构记的名与声明名对得上吗（记为别种 / 基名属 / 异名（非错））",
    "这个号有模式地位标注吗（未核 —— 没有标注不等于不是模式菌株）",
]
SCOPE_OUTSIDE = [
    "这个号的序列真的属于它声称的种吗（需序列轴）",
    "这株菌到底是不是模式菌株（需原始文献/模式名录）",
    "论文的声明对不对（需文献轴）",
    "号对但菌株错（需序列轴 + 实物核对）",
]


def dig(v):
    """拆嵌套值 [{'TemplateId':0,'Name':{'Value':'x'}}] → 'x'"""
    if v is None:
        return ""
    s = str(v).strip()
    if s[:2] in ("[{", "{'"):
        try:
            obj = ast.literal_eval(s)
            acc = []
            def walk(o):
                if isinstance(o, dict):
                    for k, vv in o.items():
                        if k == "Value" and not isinstance(vv, (dict, list)):
                            if str(vv).strip():
                                acc.append(str(vv))
                        else:
                            walk(vv)
                elif isinstance(o, (list, tuple)):
                    for vv in o:
                        walk(vv)
            walk(obj)
            return " | ".join(acc) if acc else s
        except Exception:
            return s
    return s


def load_records(tsv):
    rows = []
    with io.open(tsv, encoding="utf-8") as f:
        head = f.readline().rstrip("\n").split("\t")
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) == len(head):
                rows.append(dict(zip(head, p)))
    return rows


PLACEHOLDERS = {"[]", "[ ]", "?", "-", "—", "na", "n/a", "null", "none", "nan", "undefined"}


def chain_accepted(myb, name, steps=6):
    """沿 MycoBank 接受名链走到不动点，返回链条。MycoBank 的接受名未必是终点（实测过）。"""
    seq, cur = [name], name
    for _ in range(steps):
        nxt = (myb.current_name(cur) or "").strip()
        nxt = " ".join(nxt.split()[:2]).lower()
        if not nxt or nxt == cur.lower():
            break
        seq.append(nxt)
        cur = nxt
    return seq


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cbs-tsv", default=None, help="CBS 拉取结果（默认自动找 out/cbs_*.tsv）")
    ap.add_argument("--genus", default="Colletotrichum")
    ap.add_argument("--fields", default=os.path.join(REPO, "_cbs_fields.json"))
    ap.add_argument("--out", default=os.path.join(REPO, "out", "quality_report.json"))
    a = ap.parse_args()

    tsv = a.cbs_tsv
    if not tsv:
        # ⚠️ 必须优先挑【属相关】的那份。本机实测：目录里同时有
        #    cbs_Colletotrichum_*.tsv（按属拉）和 cbs_live_*.tsv（全量拉），
        #    而 sorted() 是大小写敏感的（大写 C < 小写 l），取 [-1] 恰好挑到全量那份，
        #    于是「属级召回缺口」永远算 0 —— 数据源选错，检查再对也没用。
        allf = [f for f in os.listdir(os.path.join(REPO, "out"))
                if re.match(r"^cbs_.*\.tsv$", f) and "diff" not in f]
        prefer = [f for f in allf if a.genus.lower() in f.lower()]
        pool = prefer or allf
        tsv = os.path.join(REPO, "out", sorted(pool)[-1]) if pool else None
    print("=" * 96)
    print("聚合质量检查（F14–F18）")
    print("=" * 96)
    print("  数据源：{}".format(tsv if tsv and os.path.exists(tsv) else "**没有 CBS 拉取结果**"))
    report = {}

    from typeaudit import sources as SN
    cfg = SN.load(os.path.join(REPO, "config", "sources.json"))
    myb = cfg.get("mycobank")
    cbs = cfg.get("cbs")

    records = load_records(tsv) if (tsv and os.path.exists(tsv)) else []
    TARGET = a.genus.lower()

    # ---------------- F14 属级召回缺口 ----------------
    # 【不能拿本表自身算】——本表就是「按属名单词筛」的结果，里面的正名属必然等于目标属，
    # 于是这个检查永远算 0。本机实测踩过这个坑。
    # 正确做法：**照号去测**——从语料取「应该属于本属、但按属名单词筛捞不到」的号，
    # 逐个查它 (a) 在不在按属名单词筛的结果里 (b) 名字经 MycoBank 是否落在目标属。
    #   (a)=不在 且 (b)=在  ⇒ 确认是召回缺口。
    in_pull = set()
    for r in records:
        in_pull.add(re.sub(r"[^A-Z0-9\.]", "", (r.get("cbs_number") or "").upper()))
    corpus = []
    cp = os.path.join(REPO, "data", "known_errors.tsv")
    if os.path.exists(cp):
        with io.open(cp, encoding="utf-8") as f:
            head = f.readline().rstrip("\n").split("\t")
            for line in f:
                q = line.rstrip("\n").split("\t")
                if len(q) >= len(head):
                    corpus.append(dict(zip(head, q)))
    miss = []
    probed = 0
    for row in corpus:
        if row.get("期望标记") != "召回缺口":
            continue
        cbs_name = (row.get("CBS名") or "").strip().lower()
        accepted = chain_accepted(myb, cbs_name)[-1] if (cbs_name and myb) else ""
        belongs = bool(accepted) and accepted.startswith(TARGET + " ")
        if not belongs:
            continue
        # 「号」列可以是逗号分隔的一组号
        for num in [x.strip() for x in (row.get("号") or "").split(",") if x.strip()]:
            probed += 1
            key = re.sub(r"[^A-Z0-9\.]", "", num.upper())
            if key not in in_pull:
                miss.append({"cbs_number": num, "cbs_name": cbs_name,
                             "accepted": accepted, "in_name_filtered_pull": False})
    report["F14"] = {"check": "属级召回缺口（照号测：按属名单词筛捞不到但在目标属内）",
                     "n": len(miss), "sample": miss[:10], "probed": probed,
                     "method": "语料里的号 → 查它在不在按属名单词筛的结果里 + MycoBank 判属"}

    # ---------------- F15 受控词表 ----------------
    # ⚠️ 这条【不能全自动】：没有参考词典就判不出拼写错。
    #    能做的是①导出完整词表供人过一遍 ②自动标出词表内部的近似对。
    vocab = {}
    near_pairs = []
    if os.path.exists(a.fields):
        try:
            flds = json.loads(io.open(a.fields, encoding="utf-8").read())
            if not isinstance(flds, list):
                flds = flds.get("Data") or []
            for f in flds:
                st = f.get("States") or []
                if st:
                    vocab[(f.get("UserName") or "").strip()] = list(st)
            # ⚠️ 相似对必须收窄到【同一字段内部】，并剔掉含数字/斜杠的取值。
            #    第一版跨全部字段两两比较，结果全是垃圾对：
            #        '1 months' ~ '2 months'、'0.007' ~ '0.014'、'0 months' ~ '12 months'
            #    —— 这类「相似」没有任何信息量，只会把真正的可疑项淹掉。
            for fname, vals in vocab.items():
                toks = sorted({str(t).strip() for t in vals if str(t).strip()})
                for i in range(len(toks)):
                    for j in range(i + 1, len(toks)):
                        x, y = toks[i].lower(), toks[j].lower()
                        if re.search(r"[0-9/\\]", x + y):
                            continue
                        if abs(len(x) - len(y)) > 2 or min(len(x), len(y)) < 5:
                            continue
                        if _lev(x, y, 2) <= 2:
                            near_pairs.append((fname, toks[i], toks[j]))
        except Exception as e:
            print("  [F15] 读字段表失败：{}".format(str(e)[:60]))
    report["F15"] = {"check": "受控词表（半自动：自动出清单，人工判定拼写）",
                     "auto_can_detect": False,
                     "n": len(vocab), "n_vocab_fields": len(vocab),
                     "n_tokens": sum(len(v) for v in vocab.values()),
                     "near_pairs": near_pairs[:10],
                     "vocab": vocab,
                     "boundary": "拼写错需要参考词典才能自动判；本检查只导出完整词表并标相似对，须人工过一遍"}

    # ---------------- F16 权威分歧 ----------------
    agree = []
    for r in records:
        nm = (r.get("taxon_name") or "").strip().lower()
        if not nm or not myb:
            continue
        final = chain_accepted(myb, nm)[-1]
        if final and final != nm:
            g1, g2 = nm.split()[0], final.split()[0]
            if g1 != g2:
                agree.append({"cbs_number": r.get("cbs_number", ""), "cbs_name": nm,
                              "accepted": final})
    by_genus = Counter(x["accepted"].split()[0] for x in agree)
    report["F16"] = {"check": "权威分歧（CBS 名的属 ≠ MycoBank 接受名的属）",
                     "n": len(agree), "by_accepted_genus": dict(by_genus.most_common(10)),
                     "sample": agree[:10]}

    # ---------------- F17 目录数据质量 ----------------
    cov = {}
    if cbs is not None:
        n_all = len(getattr(cbs, "rows", []) or [])
        n_with = sum(1 for r in getattr(cbs, "rows", [])
                     if (r.get("has_status") or "").strip() == "1")
        cov = {"rows": n_all, "with_status": n_with,
               "without_status": n_all - n_with,
               "coverage_pct": round(100.0 * n_with / n_all, 1) if n_all else None,
               "parsed_as_number": getattr(cbs, "n_parsed", None),
               "skipped_non_number": getattr(cbs, "n_skipped", None),
               "collision_keys": len(getattr(cbs, "collisions", {}) or {})}
    report["F17"] = {"check": "CBS 记录自身质量（覆盖率 / 非号码行 / 撞车键）", "n": cov.get("without_status", 0),
                     "detail": cov}

    # ---------------- F18 占位符被当数据 ----------------
    ph = []
    for r in records:
        for col in ("taxon_name", "status", "taxon_full", "organism", "other_coll"):
            v = dig(r.get(col)).strip()
            if v.lower() in PLACEHOLDERS:
                ph.append({"cbs_number": r.get("cbs_number", ""), "column": col, "value": v})
    report["F18"] = {"check": "占位符被当数据（[]/?/-/NA/null…）", "n": len(ph), "sample": ph[:10]}

    # ---------------- F28 共享外部保藏号（同株重录 / 并名可疑） ----------------
    # 两条记录指向【同一个外部保藏号】：
    #     种名相同 ⇒ 同株重录候选（同一个菌被登记了两次）
    #     种名不同 ⇒ 并名可疑（同一个外部号被归到两个种）
    # ⚠️ 分词必须按「代码+号码」整体匹配，且【号码部分要允许 `/` 和 `.`】。
    #    本机实测：号码字符类漏了 `/`，于是 `PD 88/858`、`PD 88/673`… 全被截成 `PD88`，
    #    7 条不同记录挤成同一个键，凭空造出 9 条假发现。
    #    判别特征：假发现【成簇】出现，且共享键总是同一段前缀（PD88/PD92…）；
    #              真共享的键是完整号（FMR6273、HKUCC10928），不会挤在一起。
    TOK = re.compile(r"([A-Za-z][A-Za-z\.\-]{0,12})\s*[:#]?\s*([0-9][\w\.\-/]{0,20})")
    _nn = lambda x: re.sub(r"[^A-Z0-9]", "", (x or "").upper())
    share = defaultdict(set)
    for r in records:
        for code, num in TOK.findall(r.get("other_coll") or ""):
            c = code.upper().rstrip(".")
            if c in ("AND", "OR", "THE", "SEE"):
                continue
            share[_nn(c + num)].add(_nn(r.get("cbs_number")))
    bycbs = {_nn(r.get("cbs_number")): r for r in records}
    sh_same, sh_diff = [], []
    for k, v in sorted(share.items()):
        if len(v) < 2:
            continue
        taxa = {_nn(bycbs[c].get("taxon_name")) for c in v if c in bycbs}
        rec = {"shared_number": k, "n": len(taxa),
               "records": [{"cbs": bycbs[c].get("cbs_number"),
                            "taxon": bycbs[c].get("taxon_name"),
                            "other_coll": bycbs[c].get("other_coll")}
                           for c in sorted(v) if c in bycbs]}
        (sh_same if len(taxa) == 1 else sh_diff).append(rec)
    report["F28"] = {"check": "共享外部保藏号（同种名=同株重录 / 异种名=并名可疑）",
                     "n": len(sh_same) + len(sh_diff),
                     "n_same_taxon": len(sh_same), "n_diff_taxon": len(sh_diff),
                     "sample": (sh_diff + sh_same)[:8],
                     "multi_other_coll": sum(1 for r in records if (r.get("other_coll") or "").strip())}

    # ---------------- F29 并名关系自洽性（MycoBank 自查） ----------------
    # A→B 且 B→A 逻辑上不可能 —— 两个不能都是「当前名」。成环同理。
    # 这是**新错误类型**：纯 MycoBank 自查，不依赖任何序列，几秒跑完。
    # ⚠️ 变量名不能叫 `a` —— main() 里的 `a` 是 argparse 命名空间，
    #    在这里覆盖它会让后面 `a.out` 崩掉（本机实测踩过）。
    # ⚠️ 归并必须用【集合】而不是「每名留一个代表」：
    #    MycoBank 里一个两词名可能对应多个全名、各有不同的接受名，
    #    「留一个代表」的做法取决于字典遍历顺序 ⇒ 计数不稳定（实测 13 vs 49）。
    # 归并口径由适配器统一提供（species_cur）—— 检查项自己截断会导致计数不稳定
    cur = dict(getattr(myb, "species_cur", {}) or {}) if myb is not None else {}
    # 互相并名：X 的接受名里含 Y，且 Y 的接受名里含 X —— 与遍历顺序无关
    recip = sorted({tuple(sorted((x, y)))
                    for x, ys in cur.items() for y in ys
                    if x in cur.get(y, ()) and x < y})
    same_genus = [x for x in recip if x[0].split()[0] == x[1].split()[0]]
    diff_genus = [x for x in recip if x[0].split()[0] != x[1].split()[0]]
    report["F29"] = {"check": "并名关系自洽性（互相并名 = 逻辑不可能）",
                     "n": len(recip), "n_same_genus": len(same_genus),
                     "n_diff_genus": len(diff_genus),
                     "diff_genus": [{"a": a, "b": b} for a, b in diff_genus],
                     "same_genus_sample": [{"a": a, "b": b} for a, b in same_genus[:12]],
                     "mapping": len(cur),
                     "boundary": "本检查只看 MycoBank 自身的关系是否自洽，不判断哪个名字对"}

    # ---------------- F19 归一化撞车 ----------------
    # 这条盯的是【工具自己】的归一键：删点曾把 CBS 100.06 与 CBS 10006 压成同一个键，
    # 目录查询时会静默张冠李戴。现在键保留小数点，撞车只剩真正的同号重复。
    n_coll = (cov.get("collision_keys") if cov else None)
    report["F19"] = {"check": "归一化撞车（保留小数点后剩下的真重复键）",
                     "n": n_coll if isinstance(n_coll, int) else 0,
                     "detail": {"collision_keys": n_coll,
                                "note": "修前 2,709 组（删点造成）；现在只剩真重复"}}

    # ---------------- 输出 ----------------
    def line(fid, label):
        R = report[fid]
        print("  {:<5} {:<40} 命中 {}{}".format(
            fid, label, R["n"],
            "   ← 尚未实现/半自动" if fid == "F15" else ""))

    print()
    line("F14", "属级召回缺口")
    line("F15", "受控词表（半自动）")
    line("F16", "权威分歧")
    line("F17", "CBS 记录自身质量")
    line("F18", "占位符被当数据")
    line("F19", "归一化撞车")
    line("F28", "共享外部保藏号")
    line("F29", "并名关系自洽性")
    print()
    print("  F14 样例：")
    for x in report["F14"]["sample"][:6]:
        print("     {:<14} {:<30} → 接受名 {}".format(
            x["cbs_number"], x["cbs_name"][:30], x["accepted"][:30]))
    print("  F16 按接受名属分布：{}".format(report["F16"]["by_accepted_genus"]))
    print("  F17 CBS 记录：{}".format(report["F17"]["detail"]))
    print("  F18 样例：{}".format(report["F18"]["sample"][:4]))
    print("  F28 共享 {} 个外部号（同种名 {} / 异种名 {}）：".format(
        report["F28"]["n"], report["F28"]["n_same_taxon"], report["F28"]["n_diff_taxon"]))
    for x in report["F28"]["sample"][:4]:
        print("     {:<14} {}".format(x["shared_number"], " ｜ ".join(
            "{}:{}".format(r["cbs"], (r["taxon"] or "")[:26]) for r in x["records"])))
    print("  F29 互相并名 {} 对（异属 {} ⇒ 逻辑不可能；同属 {} ⇒ 多为拼写变体对）".format(
        report["F29"]["n"], report["F29"]["n_diff_genus"], report["F29"]["n_same_genus"]))
    for x in report["F29"]["diff_genus"][:4]:
        print("     {} ⇄ {}".format(x["a"], x["b"]))
    _f15 = report["F15"]
    print("  F15 词表字段数 {} / 取值数 {}；相似对 {} 组".format(
        _f15["n_vocab_fields"], _f15["n_tokens"], len(_f15["near_pairs"])))

    print()
    print("=" * 96)
    print("诚实边界（报告须原样带出）")
    print("=" * 96)
    print("  能答：")
    for x in SCOPE_INSIDE:
        print("     + " + x)
    print("  不能答：")
    for x in SCOPE_OUTSIDE:
        print("     - " + x)
    report["_scope"] = {"inside": SCOPE_INSIDE, "outside": SCOPE_OUTSIDE,
                        "statement": "本工具是高精度候选生成器，不是穷尽查错器。空表 ≠ 没有错。"}

    d = os.path.dirname(a.out)
    if d and not os.path.isdir(d):
        os.makedirs(d)
    with io.open(a.out, "w", encoding="utf-8") as f:
        json.dump(report, f, ensure_ascii=False, indent=1)
    print("\n  已写 {}".format(a.out))
    return 0


def _lev(a, b, cap):
    if a == b:
        return 0
    if abs(len(a) - len(b)) > cap:
        return cap + 1
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i]
        for j, cb in enumerate(b, 1):
            cur.append(min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (ca != cb)))
        prev = cur
    return prev[-1]


if __name__ == "__main__":
    sys.exit(main())
