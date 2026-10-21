#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""一键补全：把三张"旁证表"自动跑出来，再汇总成一张发现表。

为什么分两种节奏（这是设计要点，不是偷懒）
--------------------------------------
  · 四方核验、反向审计 —— **按种**做，可以在采集完该种后立刻做，开销小（每物种 1–3 次查询）
  · 号连接            —— **必须等号表成型**：跨号、跨种求连通分量，
                         而且要逐个号查 CBS 的其他机构号 / ICMP 的 Duplicates，
                         开销按【号数】增长（1.5–4 秒/号），所以单独开关 `--link`

用法
----
  python tools/enrich.py --out out_audit --genus colletotrichum            # 快档（不含号连接）
  python tools/enrich.py --out out_audit --genus colletotrichum --link     # 全档（含号连接）
  python tools/enrich.py --out out_audit --skip-reverse                    # 跳过反向审计
"""
import argparse, io, json, os, subprocess, sys, time

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(REPO, "src"))

def run_py_script(script, args, log, label=None):
    """在【进程内】执行另一个脚本的 main()。

    ⚠️ 为什么不能 subprocess + sys.executable：
       程序被 PyInstaller 打包后，`sys.executable` 指向 **exe 本身**，不是 Python。
       于是 Popen([sys.executable, "xxx.py", ...]) 实际上是「又启动一次 exe」——
       它检测到服务已在运行、开个浏览器、退出码 0。
       结果：功能显示成功、实际什么都没做（实测「一键补全」就是这样假成功的）。
    """
    import contextlib, importlib.util, sys as _sys
    p = os.path.join(HERE, script)
    label = label or script
    log("{} …".format(label))
    t0 = time.time()

    class _W:
        def write(self, s):
            if s and s.strip():
                log("   " + s.rstrip()[:400])
            return len(s)
        def flush(self):
            pass

    old_argv = _sys.argv
    try:
        spec = importlib.util.spec_from_file_location(
            "_inproc_" + os.path.splitext(os.path.basename(script))[0], p)
        m = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(m)
        _sys.argv = [p] + [str(x) for x in (args or [])]
        # ⚠️ stderr 也要捕获：argparse 的错误、traceback 都走 stderr，
        #    只重定向 stdout 会让「因参数错而退出」看起来像「跑完但没输出」（实测 rc=2 却零日志）。
        with contextlib.redirect_stdout(_W()), contextlib.redirect_stderr(_W()):
            rc = m.main()
        log("  → 完成（{}，{:.0f}s）".format(rc, time.time() - t0))
        return rc
    except SystemExit as e:
        log("  → 退出（{}，{:.0f}s）".format(e.code, time.time() - t0))
        return e.code
    except Exception as e:
        log("  ❌ {}: {}".format(type(e).__name__, str(e)[:200]))
        return 1
    finally:
        _sys.argv = old_argv


SPEC_COLS = ["种", "MycoBank接受名", "NCBI_taxid", "NCBI类型材料号", "CBS类型号",
             "CGMCC该名号", "CGMCC模式号", "矛盾数", "待核数", "提示数",
             "矛盾明细", "待核明细", "名称分歧"]


def read_species(out_dir, genus):
    p = os.path.join(out_dir, "{}_三列表.tsv".format(genus.lower()))
    if not os.path.exists(p):
        return []
    return [l.split("\t")[0].strip()
            for l in io.open(p, encoding="utf-8-sig").read().splitlines()[1:] if l.strip()]


def read_numbers(out_dir, genus):
    """从三列表两列里抽全部号（含标记的去掉标记）。"""
    p = os.path.join(out_dir, "{}_三列表.tsv".format(genus.lower()))
    seen, out = set(), []
    if not os.path.exists(p):
        return out
    import re as _re
    for line in io.open(p, encoding="utf-8-sig").read().splitlines()[1:]:
        f = line.split("\t")
        for ci in (1, 2):
            if len(f) <= ci or f[ci] in ("—", ""):
                continue
            for cell in f[ci].split(" ｜ "):
                m = _re.match(r"^(.*?)\s*\[", cell)
                disp = (m.group(1) if m else cell).strip()
                if not disp or disp.startswith("(strain)"):
                    continue
                k = _re.sub(r"\s+", "", disp).upper()
                if k in seen or not _re.search(r"\d", k):
                    continue
                seen.add(k)
                out.append(disp)
    return out


def step_species(out_dir, genus, srcs, log):
    """四方核验：按种问 MycoBank / NCBI / CBS / CGMCC 并对账。"""
    from typeaudit import sources as S
    names = read_species(out_dir, genus)
    log("四方核验：{} 个种".format(len(names)))
    rows = []
    t0 = time.time()
    for i, nm in enumerate(names, 1):
        try:
            st = S.species_status(nm, srcs)
        except Exception as e:
            log("  [{}] {} 失败 {}".format(i, nm, type(e).__name__))
            continue
        nc = st.get("number_check") or {}
        rows.append({
            "种": nm,
            "MycoBank接受名": ((st.get("myco") or {}).get("current") or "")[:90],
            "NCBI_taxid": str((st.get("ncbi") or {}).get("taxid") or ""),
            "NCBI类型材料号": "; ".join((nc.get("ncbi_type_numbers") or [])[:5]),
            "CBS类型号": "; ".join((nc.get("cbs_type_numbers") or [])[:4]),
            "CGMCC该名号": "; ".join((nc.get("cgmcc_all") or [])[:4]),
            "CGMCC模式号": "; ".join((nc.get("cgmcc_type") or [])[:3]),
            "矛盾数": len(st.get("conflicts") or []),
            "待核数": len(st.get("pending") or []),
            "提示数": len(st.get("notes") or []),
            "矛盾明细": " ｜ ".join((st.get("conflicts") or [])[:2])[:300],
            "待核明细": " ｜ ".join((st.get("pending") or [])[:2])[:300],
            "名称分歧": json.dumps(st.get("name_conflict"), ensure_ascii=False) if st.get("name_conflict") else "",
        })
        if i % 20 == 0:
            log("  [{}/{}] {:.0f}s".format(i, len(names), time.time()-t0))
    p = os.path.join(out_dir, "四方核验.tsv")
    with io.open(p, "w", encoding="utf-8", newline="") as f:
        f.write("\t".join(SPEC_COLS) + "\n")
        for r in rows:
            f.write("\t".join(str(r.get(c, "")).replace("\t", " ") for c in SPEC_COLS) + "\n")
    n_c = sum(1 for r in rows if r["矛盾数"])
    n_p = sum(1 for r in rows if r["待核数"])
    log("  → {}（{} 种，矛盾 {}，待核 {}，{:.0f}s）".format(
        os.path.basename(p), len(rows), n_c, n_p, time.time()-t0))
    return len(rows)


def step_link(out_dir, genus, srcs, log, icmp_budget=200):
    from typeaudit import sources as S
    nums = read_numbers(out_dir, genus)
    log("号连接：{} 个号（ICMP 预算 {}）".format(len(nums), icmp_budget))
    t0 = time.time()
    res = S.link_numbers(nums, srcs, icmp_budget=icmp_budget,
                         progress=lambda n: log("   ICMP 已查 {} 个".format(n)))
    import re as _re
    p = os.path.join(out_dir, "号连接组.tsv")
    with io.open(p, "w", encoding="utf-8", newline="") as f:
        f.write("组号\t成员数\t成员\t机构代码\t连接来源\n")
        for i, g in enumerate(res.get("groups") or [], 1):
            codes = sorted({_re.match(r"^([A-Za-z]+)", x).group(1).upper()
                            for x in g if _re.match(r"^([A-Za-z]+)", x)})
            vias = sorted({(e.get("via") or "") for e in (res.get("edges") or [])
                           if e.get("a") in g and e.get("b") in g})
            f.write("{}\t{}\t{}\t{}\t{}\n".format(i, len(g), "; ".join(g),
                                                  "; ".join(codes), "; ".join(vias)))
    log("  → {}（{} 组，孤立 {}，{:.0f}s）".format(
        os.path.basename(p), len(res.get("groups") or []),
        len(res.get("singletons") or []), time.time()-t0))
    return len(res.get("groups") or [])



def step_tree(out_dir, genus, log, max_pairs=3, phylo_bin=None):
    """④ 自动建树：对「两皆接受名」的候选逐对跑三级裁决，直接给判断。

    触发条件很窄：只有【目录记名与声明种都是现行接受名】那一档才需要建树
    —— 其余类别（异名、无目录、次要模式…）建树也答不了，白烧时间。

    成本：每对 1–5 分钟（NCBI 取数 + MAFFT + IQ-TREE），所以有上限 max_pairs，
    被跳过的会如实报出来 —— 让读者知道"没跑"不等于"跑了没结论"。
    """
    import json as _json, os as _os, re as _re
    p = _os.path.join(out_dir, "{}_发现汇总.tsv".format(genus.lower()))
    if not _os.path.exists(p):
        log("建树：没有发现汇总表，跳过")
        return 0
    lines = io.open(p, encoding="utf-8").read().splitlines()
    if not lines:
        return 0
    hdr = lines[0].split("\t")
    cand = []
    for l in lines[1:]:
        d = dict(zip(hdr, l.split("\t")))
        if "两皆接受名" not in (d.get("类别") or ""):
            continue
        obj = d.get("对象") or ""
        fnd = d.get("发现") or ""
        if " / " not in obj:
            continue
        declared = obj.split(" / ")[0].strip().split()[-1]
        # 两种写法都要认：名字层写「目录记名 X」，号层标记写「目录记 X」。
        # 只写「目录记名」会把号层那批全漏掉（实测候选 0 对）。
        # 认三种写法：老版「目录记 X」、新版「CBS 记 X / ICMP 记 X」、以及名字层「目录记名 X」。
        # ⚠️ 只写「目录记」会在改名后一条都匹配不上（实测候选从 6 对掉到 0 对）。
        m = _re.search(r"(?:CBS|ICMP|目录|机构)\s*记\s*名?\s*([A-Z][a-z]+\s+[a-z][a-z\-]+)", fnd)
        if not m:
            continue
        cat_epi = m.group(1).split()[-1].lower()
        if cat_epi and cat_epi != declared.lower():
            cand.append((declared, cat_epi, d.get("对象")))
    log("建树：候选 {} 对（两皆接受名）".format(len(cand)))
    # 【环境确认】先过 auto_verdict 的 probe；不通就整段跳过，绝不硬跑。
    try:
        import importlib.util as _ilu2
        _sp2 = _ilu2.spec_from_file_location(
            "_av", os.path.join(HERE, "auto_verdict.py"))
        _av = _ilu2.module_from_spec(_sp2); _sp2.loader.exec_module(_av)
        _ok, _why = _av.phylo_ready(phylo_bin or
                                    "/home/yelu/miniforge3/envs/phylo/bin")
    except Exception as _e2:
        _ok, _why = False, "环境确认脚本本身出错：{}".format(_e2)
    if not _ok:
        log("  ❌ 环境确认未通过，跳过建树：{}".format(_why))
        log("     （不硬跑：环境不全时建树的失败会被误读成'数据不够'）")
        return 0
    log("  ✔ " + _why)
    if not cand:
        return 0
    use = cand[:max_pairs]
    if len(cand) > len(use):
        log("  只跑前 {} 对，其余 {} 对未跑（max_pairs={}）".format(
            len(use), len(cand) - len(use), max_pairs))
    pairs_file = _os.path.join(out_dir, "_待建树对.tsv")
    with io.open(pairs_file, "w", encoding="utf-8", newline="") as f:
        for a, b, _o in use:
            f.write("{}\t{}\n".format(a, b))
    for a, b, obj in use:
        log("  · {}  vs  {}".format(a, b))
    run_py_script("auto_verdict.py",
                  ["--pairs-file", pairs_file,
                   "--out", _os.path.join(REPO, "tree_auto"),
                   "--min-cols", "120", "--per-name", "150"],
                  log, "自动三级裁决（{} 对）".format(len(use)))
    # 回读裁决结果
    sp = _os.path.join(REPO, "tree_auto", "自动裁决汇总.tsv")
    if _os.path.exists(sp):
        log("  裁决结果：")
        for l in io.open(sp, encoding="utf-8").read().splitlines()[1:]:
            f = l.split("\t")
            log("    {:<28} 层级 {} ｜ {} 位点 ｜ {}".format(
                f[0][:28], f[1], f[2], f[3]))
    return len(use)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="out_audit")
    ap.add_argument("--genus", default="colletotrichum")
    ap.add_argument("--link", action="store_true", help="含号连接（慢，按号数增长）")
    ap.add_argument("--skip-reverse", action="store_true")
    ap.add_argument("--skip-summary", action="store_true")
    ap.add_argument("--icmp-budget", type=int, default=200)
    ap.add_argument("--email", default=os.environ.get("NCBI_EMAIL") or "audit@example.org",
                    help="NCBI 要求的邮箱（序列轴、反向审计都要）")
    ap.add_argument("--skip-seq", action="store_true", help="跳过序列轴（同号多名扫描）")
    ap.add_argument("--seq-max", type=int, default=80,
                    help="序列轴最多检查多少个种（约 6 号/种，0.5 秒/号）")
    ap.add_argument("--phylo-bin",
                    default="/home/yelu/miniforge3/envs/phylo/bin",
                    help="WSL 里 mafft/iqtree2/trimal 所在目录")
    ap.add_argument("--tree", type=int, default=3,
                    help="自动建树的候选对数上限（0=不建树）；每对 1–5 分钟")
    a = ap.parse_args()
    if getattr(a, "email", ""):
        os.environ["NCBI_EMAIL"] = a.email

    out_dir = a.out if os.path.isabs(a.out) else os.path.join(REPO, a.out)
    from typeaudit import sources as S
    srcs = S.load(os.path.join(REPO, "config", "sources.json"))

    def log(s):
        print(s, flush=True)

    t0 = time.time()
    # ⚠️ 必须写进【真实环境变量】，不能只写子进程 env：
    #    改成进程内调用后，"传给子进程的 env" 这条路径已经不存在了，
    #    而 reverse_audit 靠 os.environ["NCBI_EMAIL"] 拿邮箱，缺了就静默 return 2。
    if not os.environ.get("NCBI_EMAIL"):
        os.environ["NCBI_EMAIL"] = "researcher@example.org"
        log("（NCBI_EMAIL 未设，暂用 researcher@example.org）")
    log("=== 一键补全 {} / {} ===".format(os.path.basename(out_dir), a.genus))
    step_species(out_dir, a.genus, srcs, log)
    if a.link:
        step_link(out_dir, a.genus, srcs, log, a.icmp_budget)
    else:
        log("号连接：跳过（加 --link 才跑；按号计费，1.5–4 秒/号）")
    if not a.skip_reverse:
        run_py_script("reverse_audit.py", ["--genus", "Colletotrichum"], log, "反向审计")
    if not a.skip_summary:
        run_py_script("collect_findings.py", ["--out", a.out, "--genus", a.genus],
                  log, "汇总发现")
    # 序列轴：同一保藏号在 GenBank 里被标成几个种（便宜、证据直接）
    # 放在汇总之后、建树之前：它本身就是"并名"的最直接证据，建树是它的补充。
    if not getattr(a, "skip_seq", False):
        log("序列轴（同号多名）…")
        run_py_script("seq_axis.py",
                      ["--out", a.out, "--genus", a.genus, "--mode", "same_number",
                       "--max-species", str(getattr(a, "seq_max", 80)),
                       "--email", a.email or "audit@example.org"],
                      log, "序列轴·同号多名")
        run_py_script("seq_axis.py",
                      ["--out", a.out, "--genus", a.genus, "--mode", "type_collision",
                       "--max-species", str(getattr(a, "seq_max", 80)),
                       "--email", a.email or "audit@example.org"],
                      log, "分类库·一号多模式")
        run_py_script("collect_findings.py", ["--out", a.out, "--genus", a.genus],
                      log, "重新汇总（并入序列层）")
    else:
        log("序列轴：跳过（--skip-seq）")

    if a.tree and a.tree > 0:
        step_tree(out_dir, a.genus, log, max_pairs=a.tree, phylo_bin=a.phylo_bin)
    else:
        log("建树：跳过（--tree 0）")
    log("=== 全部完成：{:.0f}s ===".format(time.time()-t0))
    return 0


if __name__ == "__main__":
    sys.exit(main())
