#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
check_known_errors.py —— 拿【已发现的错误】当回归语料，逐条验证工具能否命中。

【为什么必须有这个】
    没有它，「工具能不能查出所有错误」只能靠嘴说。有了它，每个发现都变成一条
    可重跑的断言：命中 / 未命中 / 尚未实现。**未命中的条目就是待补的检查项** ——
    这是把「确查所有」从口号变成可验证陈述的唯一办法。

【它怎么测】
    不复制判定逻辑（复制的测试测的是复制品），而是**真跑一遍流水线**，
    再从产出的三列表里读标记。语料里每个 case 自带「宣告种」，所以能一次性
    把需要的种拼成 --species-list 跑完。

用法
    python tools/check_known_errors.py                  # 跑流水线 + 核验（需联网）
    python tools/check_known_errors.py --no-run         # 只核验已有的三列表
    python tools/check_known_errors.py --out out_known  # 换输出目录
"""
import sys
try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

import argparse, io, json, os, re, subprocess

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CORPUS = os.path.join(REPO, "data", "known_errors.tsv")

# 尚未实现的聚合检查：列在这里，让缺口一眼看得见，而不是悄悄算通过
AGG_MISSING = {"AGGREGATE_MISSING"}
AGG_OK = {"AGGREGATE_OK"}


def load_corpus(path):
    rows = []
    with io.open(path, encoding="utf-8") as f:
        head = f.readline().rstrip("\n").split("\t")
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= len(head):
                rows.append(dict(zip(head, p)))
    return rows


def norm(x):
    return re.sub(r"[^A-Za-z0-9\.]", "", (x or "").upper())


def parse_table(path):
    """
    三列表 → {(种名规范化, 号规范化): 标记集合}

    ⚠️ 标记必须按【集合】收集，不能存单个字符串：
       同一个号会同时出现在「菌株号」列和「模式菌」列，而两列的标记可以不同。
        实测 CBS 169.59：菌株号列是 [未核]、模式菌列是 [—]。
        用字典直接赋值的话，后写的会覆盖先写的，于是核验结果随列序变化 ——
        那是验证器的 bug，不是被验证对象的 bug。
    """
    got = {}
    if not os.path.exists(path):
        return got
    with io.open(path, encoding="utf-8-sig") as f:
        f.readline()
        for line in f:
            fld = line.rstrip("\n").split("\t")
            if len(fld) < 3:
                continue
            sp = norm(fld[0])
            for col in (1, 2):
                for cell in fld[col].split("｜"):
                    m = re.match(r"^(.*?)\s*\[(.*?)\]\s*$", cell.strip())
                    if not m:
                        continue
                    num, flags = m.group(1).strip(), m.group(2).strip()
                    num = re.sub(r"^\(strain\)\s*", "", num)     # 去掉 (strain) 前缀
                    got.setdefault((sp, norm(num)), set()).update(
                        x.strip() for x in flags.split(";") if x.strip())
    return got


def observed_flags(table):
    """
    扫三列表，返回实际出现过的标记集合 —— 用来和语料覆盖做差集。

    【为什么需要它 —— 这是「会不会有新东西」的唯一客观判据】
        语料是【验证器】不是【发现器】：它只能告诉你已知检查有没有失效。
        要判断「跑完有没有语料覆盖不到的东西」，只能拿【本次实际出现的标记集合】
        去减【语料覆盖的标记集合】。差集非空 ⇒ 有东西需要新增语料行 / 新增检查。

        本机实测的差集样例：「目录记为别种」出现在旧代码跑出的三列表里，
        而语料里只有新代码的「异名（非错）」—— 同一个号、两次跑出不同标记，
        差集就是这么来的。它不是 bug，是入口。
    """
    from collections import Counter
    got = Counter()
    if not os.path.exists(table):
        return got
    with io.open(table, encoding="utf-8-sig") as f:
        f.readline()
        for line in f:
            fld = line.rstrip("\n").split("\t")
            if len(fld) < 3:
                continue
            for col in (1, 2):
                for cell in fld[col].split("｜"):
                    m = re.match(r"^(.*?)\s*\[(.*?)\]\s*$", cell.strip())
                    if not m:
                        continue
                    for flag in [x.strip() for x in m.group(2).split(";") if x.strip()]:
                        if flag != "—":            # 「—」是「无问题」，不是错误类
                            got[flag] += 1
    return got


def run_pipeline(species_list, out_dir, email):
    spf = os.path.join(out_dir, "_species.txt")
    io.open(spf, "w", encoding="utf-8").write("\n".join(species_list) + "\n")
    env = dict(os.environ)
    env["PYTHONPATH"] = os.path.join(REPO, "src")
    cmd = [sys.executable, "-u", "-m", "typeaudit.pipeline",
           "--genus", "Colletotrichum", "--species-list", spf,
           "--out", out_dir, "--config", os.path.join(REPO, "config", "sources.json"),
           "--email", email]
    print("$ " + " ".join(cmd[2:]))
    r = subprocess.run(cmd, cwd=REPO, env=env, timeout=3600,
                       text=True, encoding="utf-8", errors="replace", capture_output=True)
    for ln in (r.stdout or "").split("\n"):
        if "[源]" in ln or "可用源" in ln or "[acquire]" in ln or "[report]" in ln:
            print("   " + ln.strip())
    if r.returncode != 0:
        print((r.stderr or "")[-1200:])
    return r.returncode


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(REPO, "out_known"))
    ap.add_argument("--no-run", action="store_true")
    ap.add_argument("--email", default=os.environ.get("NCBI_EMAIL") or "audit@example.org")
    a = ap.parse_args()

    rows = load_corpus(CORPUS)
    cases = [r for r in rows if r.get("检验方式") == "case"]
    # 检验方式=findings ⇒ 到「全部发现」表里核对（名字层/记录层/目录层的结论不在三列表里）
    finds = [r for r in rows if r.get("检验方式") == "findings"]
    agg = [r for r in rows if r.get("检验方式") not in ("case", "findings")]

    if not a.no_run:
        sps = sorted(set(r["宣告种"] for r in (cases + finds) if r.get("宣告种")))
        print("=" * 96)
        print("用 {} 个种真跑一遍流水线，再逐条核验".format(len(sps)))
        print("=" * 96)
        if not os.path.isdir(a.out):
            os.makedirs(a.out)
        run_pipeline(sps, a.out, a.email)
        # 生成「全部发现」表：findings 类语料要在那张表里核
        try:
            import importlib.util as _ilu
            _sp = _ilu.spec_from_file_location("_cf", os.path.join(REPO, "tools", "collect_findings.py"))
            _m = _ilu.module_from_spec(_sp)
            _sp.loader.exec_module(_m)
            _old = sys.argv
            sys.argv = ["collect_findings.py", "--out", a.out, "--genus", "colletotrichum"]
            try:
                _m.main()
            finally:
                sys.argv = _old
        except Exception as _e:
            print("⚠️ 生成全部发现失败（findings 语料将无法核对）：{}".format(_e))

    table = os.path.join(a.out, "colletotrichum_三列表.tsv")
    got = parse_table(table)
    print()
    print("=" * 96)
    print("逐条核验（语料：{}）".format(os.path.relpath(CORPUS, REPO)))
    print("=" * 96)
    n_hit = n_miss = 0
    for r in cases:
        key = (norm(r["宣告种"]), norm(r["号"]))
        aset = got.get(key)
        exp = r["期望标记"]
        mark = "; ".join(sorted(aset)) if aset else "—"
        # 期望值以 * 结尾 ⇒ 前缀匹配。用于「标记里带动态值」的情况
        # （如 `书写错误（形近 CGMCC，但该号在 CGMCC 查无）` 的括号内容会随数据变）。
        # 期望值以 ! 开头 ⇒ 【负样本】：这条标记**必须不出现**。
        #   为什么需要：现有语料全是"要报出来"，没有一条防"不该报的别报"，
        #   而实测两个 bug 恰恰是把对的报成错的（84 个一致被误报成异名）。
        negative = exp.startswith("!")
        target = exp[1:] if negative else exp
        if negative:
            if aset is None:
                mark, verdict = "（号不在表里，负样本无法验证）", "未命中"
                n_miss += 1
                hit = False
            else:
                bad = any((a.startswith(target[:-1]) if target.endswith("*") else a == target)
                          for a in aset)
                hit = not bad
                mark = "; ".join(sorted(aset))
        elif target.endswith("*"):
            hit = any(a.startswith(target[:-1]) for a in (aset or ()))
        else:
            hit = (aset is not None) and (target in aset)
        if aset is None and not negative:
            mark, verdict = "（未出现在表里）", "未命中"
            n_miss += 1
        elif hit:
            verdict = "命中"
            n_hit += 1
        elif not (negative and aset is None):
            verdict = "未命中"
            n_miss += 1
        print("  {:<5} {:<10} {:<16} 期望={:<14} 实际={:<16} {}".format(
            r["id"], r["宣告种"].replace("Colletotrichum ", "C. "), r["号"],
            exp, mark, verdict))

    # ---------------------------------------------------------------- findings 类核验
    fp = os.path.join(a.out, "colletotrichum_发现汇总.tsv")
    frows = []
    if os.path.exists(fp):
        _L = io.open(fp, encoding="utf-8").read().splitlines()
        if _L:
            _h = _L[0].split("\t")
            for _l in _L[1:]:
                _d = dict(zip(_h, _l.split("\t")))
                frows.append((_d.get("对象", ""), _d.get("类别", "")))
    if finds:
        print()
        print("=" * 96)
        print("名字层 / 记录层 / 目录层核验（读 {}）".format(os.path.basename(fp)))
        print("=" * 96)
        for r in finds:
            num = (r.get("号") or "").upper().replace(" ", "")
            epi = (r.get("宣告种") or "").split()[-1].lower()
            cats = [c for o, c in frows
                    if num and num in o.upper().replace(" ", "") and epi in o.lower()]
            exp = r.get("期望标记") or ""
            negative = exp.startswith("!")
            target = (exp[1:] if negative else exp).rstrip("*")
            found = any(c == target or c.startswith(target) for c in cats)
            hit = (not found) if negative else found
            if hit:
                n_hit += 1
            else:
                n_miss += 1
            print("  {:<5} {:<10} {:<16} 期望={:<22} 实际={:<24} {}".format(
                r.get("id", ""), (r.get("宣告种") or "").replace("Colletotrichum ", "C. "),
                r.get("号", ""), exp, ("、".join(cats) if cats else "（无对应发现）"),
                "命中" if hit else "未命中"))

    # 聚合检查：读 quality_checks.py 的产物来判命中，不靠人工标注
    qr_path = os.path.join(REPO, "out", "quality_report.json")
    qr = {}
    if os.path.exists(qr_path):
        try:
            qr = json.loads(io.open(qr_path, encoding="utf-8").read())
        except Exception as e:
            print("  ! 读 {} 失败：{} —— 聚合类将全部报「尚未实现」".format(qr_path, e))
            qr = {}
    else:
        print("  ! 没有 {}，先跑 tools/quality_checks.py".format(qr_path))
    print()
    print("=" * 96)
    print("聚合型发现（读 out/quality_report.json 判定）")
    print("=" * 96)
    n_agg_ok = n_agg_miss = n_agg_half = 0
    for r in agg:
        fid = r["id"]
        R = qr.get(fid) or {}
        lo = (r.get("期望下限") or "").strip()
        if lo == "半自动":
            st, half = "半自动（自动出清单，人工判拼写）", True
            n_agg_half += 1
        elif R:
            n = R.get("n") or 0
            ok = (n >= int(lo)) if lo.isdigit() else True
            st = "命中 n={}{}".format(n, "" if ok else "（低于下限 {}）".format(lo))
            if ok: n_agg_ok += 1
            else: n_agg_miss += 1
        else:
            st, _ = "**尚未实现**", None
            n_agg_miss += 1
            half = False
        print("  {:<5} {:<12} {}".format(fid, st, r["发现"][:64]))

    # ---------------- 语料覆盖落差：本次出现但语料未覆盖的标记 ----------------
    got = observed_flags(table)
    # 负样本（!开头）不是"期望出现的标记"，不能算进覆盖集
    corpus_set = set(r["期望标记"].strip() for r in cases
                     if r.get("期望标记", "").strip()
                     and not r["期望标记"].strip().startswith("!"))
    corpus_set |= set(r.get("期望标记", "").strip() for r in agg if r.get("期望标记", "").strip())

    # ⚠️ 落差比对必须认 `*` 前缀约定，否则会两头报假：
    #    语料写 `书写错误*`（前缀），实际标记是 `书写错误（形近 CGMCC…）`（字面）——
    #    精确比对会把实际标记报成「语料未覆盖」、同时把 `书写错误*` 报成「有但没出现」。
    #    本机实测踩过这一对自相矛盾的假报。
    def _covered(flag):
        for e in corpus_set:
            if flag == e:
                return True
            if e.endswith("*") and flag.startswith(e[:-1]):
                return True
        return False

    def _observed(entry):
        if entry.endswith("*"):
            return any(k.startswith(entry[:-1]) for k in got)
        return entry in got

    only_obs = {k: v for k, v in got.items() if not _covered(k)}
    # 「—」是「无问题」不是标记，别把它算进「语料有但没出现」
    only_corpus = sorted(x for x in corpus_set
                         if not _observed(x)
                         and x not in ("—", "召回缺口", "权威分歧", "占位符"))
    print()
    print("=" * 96)
    print("语料覆盖落差（跑全属时，这一段才是「有没有新东西」的答案）")
    print("=" * 96)
    print("  本次实际出现的标记（{} 种）：".format(len(got)))
    for k, v in got.most_common():
        print("     {:<18} {:>5} 次{}".format(k, v, "" if _covered(k) else "   ← 语料未覆盖"))
    print()
    print("  ★ 语料未覆盖的标记：{}".format(
        "、".join("{}（{} 次）".format(k, v) for k, v in sorted(only_obs.items(), key=lambda x: -x[1]))
        if only_obs else "（无 —— 本次没有语料覆盖不到的东西）"))
    if only_obs:
        print("     含义：这些标记出现了，但语料里没有对应的期望行。")
        print("     下一步：给 data/known_errors.tsv 补行 —— 但**必须同时写「依据」**，")
        print("             依据要来自工具之外（登记表/目录/MycoBank/文献/人工判读）。")
        print("             若把依据写成「工具就是这么输出的」，那条测试永远通过、什么都没验证。")
        print("     如果某个标记属于【新类型】（不是已知类型的实例），那就是加轴的信号。")
    # ---------------- 依据审计：这是「回归测试会不会退化成同义反复」的闸门 ----------------
    # 如果某条语料的「期望值」只是抄了工具自己的输出，那条测试永远通过、且什么都没验证。
    # 所以要求每条必须写「依据」，且依据必须来自工具之外（登记表/目录/MycoBank/文献/人工判读）。
    n_ev = sum(1 for r in cases + finds + agg if (r.get("依据") or "").strip())
    n_no = sum(1 for r in cases + finds + agg if not (r.get("依据") or "").strip())
    print()
    print("  依据审计：{} 条有外部依据 / {} 条缺依据".format(n_ev, n_no))
    if n_no:
        print("  ! 缺依据的条目：{}".format(
            "、".join(r["id"] for r in cases + agg if not (r.get("依据") or "").strip())))
        print("    缺依据 = 该条只记录了工具当前的行为，**不构成验证**。")
        print("    补依据时才允许把期望值定下来；依据必须来自工具之外，")
        print("    **绝不能写成「因为工具这么输出」** —— 那是同义反复，会让测试永远通过。")
    print()
    print("  语料有、本次没出现：{}".format("、".join(only_corpus) if only_corpus else "（无）"))
    if only_corpus:
        print("     含义：要么本批种里没有这类问题，要么检查失效了 —— 值得抽查确认。")

    n_agg_missing = n_agg_miss
    print()
    print("=" * 96)
    print("结论")
    print("=" * 96)
    # ⚠️ 分母必须含 findings 那批：它们与个案类共用 n_hit / n_miss 计数，
    #    只写 len(cases) 会算出 115% 这种不可能的数（实测踩过）。
    print("  逐条类：命中 {} / 未命中 {} / 共 {} 条（三列表 {} + 发现表 {}）".format(
        n_hit, n_miss, len(cases) + len(finds), len(cases), len(finds)))
    print("  聚合类：命中 {} / 半自动 {} / 未命中 {} / 共 {} 条".format(
        n_agg_ok, n_agg_half, n_agg_missing, len(agg)))
    tot = len(cases) + len(finds) + len(agg)
    ok = n_hit + (len(agg) - n_agg_missing)
    print("  总覆盖率：{} / {} = {:.0f}%".format(ok, tot, 100.0 * ok / max(1, tot)))
    print()
    print("  ★ 「尚未实现」那几条就是待补的检查项 —— 它们不是 bug，是缺口清单。")
    return 0 if (n_miss == 0 and n_agg_missing == 0) else 2


if __name__ == "__main__":
    sys.exit(main())
