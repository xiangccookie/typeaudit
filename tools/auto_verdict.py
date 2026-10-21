#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# --- console utf8 bootstrap（Windows 控制台默认 GBK，输出 ↔ ⚠️ 🔴 会崩；详见 src/typeaudit/_console.py）
import sys as _sys
for _s in (_sys.stdout, _sys.stderr):
    try:
        _s.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass
"""自动并名裁决：争议名 → 自己下序列 → 三级降级 → 带强度标签的结论。

为什么要「自动」+「降级」
------------------------
手工给 type_synonym_verdict.py 喂一个连号段（--acc-range）只适合「一批提交」的情形；
而 MycoBank 与 NCBI 的分歧名散落在几十年、不同论文里。所以这里改成：
按【种名】自动抓记录 → 按 /culture_collection 归组 → 凑位点。

三级降级（用户指定）
--------------------
  Tier A  8 位点  ITS/ACT/CAL/CHS/GAPDH/GS/TUB2/SOD2  → 结论强度：较强
  Tier B  5 位点  ITS/TUB2/ACT/CHS/GAPDH              → 结论强度：中等
  Tier C  1 位点  仅 ITS                              → 结论强度：仅提示（不得作为裁决）

⚠️ 强度必须随降级一起降。理由：本项目实测同一面板下 8 个基因会给出 8 个不同的最近邻，
   单位点相似度不足以定种（Kennedy et al. 2022 亦指出单基因相似度常致误判）。
   所以 Tier C 的结论【只能写成「提示」，不能写「确认」】。

用法
----
  python tools/auto_verdict.py --pair celtidis karsti --workdir tree_auto/celtidis_karsti
  python tools/auto_verdict.py --pairs-file data/争议对.tsv --out tree_auto
  # 只跑某一级：--tier 5
"""
import argparse, io, json, os, re, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(REPO, "src"))

# 每级的【基因集合】是显式写死的：Tier C 就是「只用 ITS」，不是「任一位点」。
TIER_GENES = {
    "A": ["ITS", "ACT", "CAL", "CHS", "GAPDH", "GS", "TUB2", "SOD2"],
    "B": ["ITS", "TUB2", "ACT", "CHS", "GAPDH"],
    "C": ["ITS"],
}
TIERS = [(8, "A", "较强"), (5, "B", "中等"), (1, "C", "仅提示")]


def load_verdict_module():
    """按路径加载 type_synonym_verdict.py（它不在包内，只能这样 import）。"""
    import importlib.util
    p = os.path.join(HERE, "type_synonym_verdict.py")
    spec = importlib.util.spec_from_file_location("tsv_mod", p)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def fetch_raw_for_names(names, email, per_name=60, sleep=1.2, log=print, strains=None):
    """按种名抓 GenBank 记录 → {version: {L,defn,cc,strain,host,tm}}（与 parse_gb 同格式）。"""
    from Bio import Entrez
    Entrez.email = email
    mod = load_verdict_module()
    raw = {}
    for name in names:
        term = '"{}"[Organism]'.format(name)
        h = Entrez.esearch(db="nuccore", term=term, retmax=per_name)
        r = Entrez.read(h); h.close()
        ids = r.get("IdList") or []
        log("  [取数] {:<34} 检索 {} 条 → 取 {}".format(
            name.replace("Colletotrichum ", ""), r.get("Count"), len(ids)))
        for k in range(0, len(ids), 100):
            chunk = ids[k:k + 100]
            for attempt in range(3):
                try:
                    h2 = Entrez.efetch(db="nuccore", id=",".join(chunk),
                                       rettype="gb", retmode="text")
                    txt = h2.read(); h2.close()
                    break
                except Exception as e:
                    log("    批 {} 第 {} 次失败: {}".format(k, attempt + 1, str(e)[:60]))
                    time.sleep(3)
            else:
                continue
            raw.update(mod.parse_gb(txt))
            time.sleep(sleep)
    # 定向补抓：按「模式株号 + 种名」把关键那株的记录补齐
    # 为什么必须补：按种名抓到的记录会摊薄在很多菌株上，模式株自己往往凑不齐位点。
    for nm, st in zip(names, (strains or [None, None])):
        if not st:
            continue
        num = re.sub(r"[^A-Za-z0-9]", "", st)
        term = '"{}"[All Fields] AND "{}"[Organism]'.format(
            st.replace(":", " "), nm)
        try:
            from Bio import Entrez
            h = Entrez.esearch(db="nuccore", term=term, retmax=300)
            r = Entrez.read(h); h.close()
            ids = r.get("IdList") or []
            log("  [补抓] {:<10} 检索 {} 条 → 取 {}".format(st, r.get("Count"), len(ids)))
            before = len(raw)
            for k in range(0, len(ids), 100):
                chunk = ids[k:k + 100]
                h2 = Entrez.efetch(db="nuccore", id=",".join(chunk), rettype="gb", retmode="text")
                txt2 = h2.read(); h2.close()
                raw.update(mod.parse_gb(txt2))
                time.sleep(sleep)
            log("        新增 {} 条记录".format(len(raw) - before))
        except Exception as e:
            log("        补抓失败：{}".format(type(e).__name__))
    return raw


def type_strain_targets(names, email, sleep=1.2, log=print):
    """自动取每个种的【模式株保藏号】（形如 "CBS:132134"），供 verdict 定位两片叶。

    来源：该名的 type_material 记录里的 /culture_collection 字段（不是 /type_material，
    后者只有"culture from holotype of X"这样的文字，不含号码）。
    """
    from Bio import Entrez
    Entrez.email = email
    mod = load_verdict_module()
    out = []
    for n in names:
        got = None
        try:
            h = Entrez.esearch(db="nuccore",
                               term='"{}"[Organism] AND type_material[All Fields]'.format(n),
                               retmax=20)
            r = Entrez.read(h); h.close()
            ids = r.get("IdList") or []
            if ids:
                time.sleep(0.4)
                h2 = Entrez.efetch(db="nuccore", id=",".join(ids[:20]),
                                   rettype="gb", retmode="text")
                recs = mod.parse_gb(h2.read()); h2.close()
                for ver, rec in recs.items():
                    for c in (rec.get("cc") or []):
                        c = c.strip()
                        if re.match(r"^[A-Za-z]{2,8}[\s:]*\d", c):
                            got = re.sub(r"[\s]+", ":", c, count=1)
                            break
                    if got:
                        break
                if not got:
                    # 退化：很多记录只填 /strain 不填 /culture_collection
                    for ver, rec in recs.items():
                        for c in (rec.get("strain") or []):
                            c = c.strip()
                            if re.match(r"^[A-Za-z]{2,8}[\s:]*\d", c):
                                got = re.sub(r"[\s]+", ":", c, count=1)
                                break
                        if got:
                            break
        except Exception as e:
            log("  取模式号失败 {}: {}".format(n, type(e).__name__))
        log("  模式株 {} → {}".format(n.replace("Colletotrichum ", ""), got or "（未取到）"))
        out.append(got)
        time.sleep(sleep)
    return out


def loci_of_strain(raw, code, mod):
    """某株（按保藏号数字段匹配）在 raw 里覆盖了几个位点。返回 (位点数, 命中的记录数)。"""
    import re as _re
    n = mod.norm_code(code)
    num = _re.sub(r"^0+", "", _re.sub(r"^[A-Z]+", "", n))
    if not num:
        return 0, 0
    genes, nrec = set(), 0
    for ver, r in raw.items():
        keys = [mod.norm_code(c) for c in (r.get("cc") or [])] + \
               [mod.norm_code(s) for s in (r.get("strain") or [])]
        if not any(num in _re.sub(r"^0+", "", k) for k in keys if k):
            continue
        g = mod.gene_of(r.get("defn"))
        if g:
            genes.add(g); nrec += 1
    return len(genes), nrec


def run_tier(mod, workdir, min_genes, targets, email, phylo_bin, min_cols, genes=None):
    """跑 prep→align→tree→verdict，返回 verdict 是否拿到两片叶。"""
    a = argparse.Namespace(workdir=workdir, min_genes=min_genes, min_cols=min_cols,
                           pair=targets, email=email, sleep=1.2, force=False,
                           species_map=None, phylo_bin=phylo_bin, out=None,
                           genes=genes)
    for stage in ("prep", "align", "tree"):
        try:
            rv = getattr(mod, "stage_" + stage)(a)
        except Exception as e:
            return {"ok": False, "stage": stage, "err": "{}: {}".format(type(e).__name__, str(e)[:120])}
        if stage == "tree" and rv is False:
            return {"ok": False, "stage": "tree", "err": "未生成树文件"}
    try:
        v = mod.stage_verdict(a)
    except Exception as e:
        return {"ok": False, "stage": "verdict", "err": "{}: {}".format(type(e).__name__, str(e)[:120])}
    # ⚠️ 必须看【返回值】：verdict 可能正常结束但没得出任何结论（目标没出现在树里）。
    #    原来只判断「有没有抛异常」，于是把「没裁决」报成了「✅ 跑通」——假成功。
    if v is False:
        return {"ok": False, "stage": "verdict", "err": "未能定位目标菌株或未得出结论"}
    return {"ok": True}



def phylo_ready(bin_dir, log=print):
    """跑树前的【环境确认】。返回 (ok, 说明)。

    ⚠️ 为什么必须先确认，而不是"跑起来再看"：
        WSL 或 phylo 环境缺失时，失败不会以"环境缺了"的样子出现，而是以
        空比对、「至少 3 条序列」、找不到 tree.treefile、序列字符数不足……
        等间接形式出现 —— 每一条都像是"数据不够"，排查极费时，而且容易被
        误读成生物学结论。**宁可提前拒绝，也不要硬跑出一个错结论。**

    检查三件事：
      ① wsl 命令可用，且能进 bash
      ② bin_dir 下 mafft / iqtree2 / trimal 三个都可执行
      ③ 每个都能 --version（确认不是坏链接/架构不符）
    """
    import subprocess as _sp
    def _run(cmd, timeout=40):
        try:
            r = _sp.run(["wsl", "-e", "bash", "-lc", cmd],
                        capture_output=True, text=True, timeout=timeout,
                        encoding="utf-8", errors="replace")
            return r.returncode, (r.stdout or "") + (r.stderr or "")
        except FileNotFoundError:
            return -1, "本机没有 wsl 命令"
        except Exception as e:
            return -2, "{}: {}".format(type(e).__name__, e)
    ld = str(bin_dir or "").rstrip("/")
    # ① WSL 可用
    rc, out = _run("echo WSL_OK")
    if rc != 0 or "WSL_OK" not in out:
        return False, "WSL 不可用（{}）".format((out or "无输出").strip()[:120])
    # ② 三个二进制存在且可执行
    rc, out = _run("for b in mafft iqtree2 trimal; do "
                   "test -x '{}/'$b && echo OK:$b || echo MISSING:$b; done".format(ld))
    miss = [ln.split(":", 1)[1] for ln in out.splitlines() if ln.startswith("MISSING:")]
    if miss:
        return False, "{} 缺可执行文件：{}（目录 {}）".format(
            "phylo 环境", "、".join(miss), ld)
    # ③ 能跑出版本号
    rc, out = _run("'{}/mafft' --version 2>&1 | head -1; "
                   "'{}/iqtree2' --version 2>&1 | head -1; "
                   "'{}/trimal' --version 2>&1 | head -1".format(ld, ld, ld), timeout=90)
    if rc != 0:
        return False, "三个工具都能找到，但运行失败（{}）".format((out or "").strip()[:120])
    vers = [l.strip() for l in out.splitlines() if l.strip()][:3]
    return True, "环境就绪：{}".format(" ｜ ".join(vers))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pair", nargs=2, default=None, help="两个种加词，如 celtidis karsti")
    ap.add_argument("--pairs-file", default=None, help="每行 a<TAB>b 的文件")
    ap.add_argument("--workdir", default=None)
    ap.add_argument("--out", default=os.path.join(REPO, "tree_auto"))
    ap.add_argument("--tier", type=int, default=None, help="只跑指定位点数（8/5/1）")
    ap.add_argument("--per-name", type=int, default=60)
    ap.add_argument("--min-cols", type=int, default=150)
    ap.add_argument("--phylo-bin", default="/home/yelu/miniforge3/envs/phylo/bin")
    ap.add_argument("--email", default=os.environ.get("NCBI_EMAIL", "researcher@example.org"))
    ap.add_argument("--probe", action="store_true",
                    help="只做环境确认（WSL + mafft/iqtree2/trimal），不建树")
    a = ap.parse_args()

    # 【环境确认】跑树前必须先过这一关 —— 不通就不跑，别硬跑出错结论。
    ok, why = phylo_ready(a.phylo_bin)
    if a.probe:
        print(("OK  " if ok else "FAIL") + " " + why)
        return 0 if ok else 3
    if not ok:
        print("❌ 环境确认未通过，拒绝建树：")
        print("   " + why)
        print("   处置：修好 WSL/phylo 环境，或用 --phylo-bin 指向正确的 bin 目录。")
        print("   提示：本工具不会在环境不全时强行建树（那会产出一堆看似'数据不够'的假结论）。")
        return 3
    print("✔ " + why)

    pairs = []
    if a.pair:
        pairs = [tuple(a.pair)]
    elif a.pairs_file:
        for line in io.open(a.pairs_file, encoding="utf-8"):
            f = line.strip().split("\t")
            if len(f) >= 2 and f[0] and not f[0].startswith("#"):
                pairs.append((f[0].strip(), f[1].strip()))
    if not pairs:
        print("没有输入对（--pair 或 --pairs-file）"); return 2

    mod = load_verdict_module()
    tiers = [(t, lab, s) for t, lab, s in TIERS if (a.tier is None or t == a.tier)]

    rows = []
    for (na, nb) in pairs:
        names = ["Colletotrichum " + na.split()[-1], "Colletotrichum " + nb.split()[-1]]
        wd = a.workdir or os.path.join(a.out, "{}__{}".format(na.split()[-1], nb.split()[-1]))
        os.makedirs(wd, exist_ok=True)
        print("\n" + "="*74)
        print("▌ {}  vs  {}".format(names[0], names[1]))
        targets = type_strain_targets(names, a.email)
        rawp = os.path.join(wd, "raw.json")
        if os.path.exists(rawp):
            raw = json.load(io.open(rawp, encoding="utf-8"))
            print("  [取数] 命中缓存：{} 条记录".format(len(raw)))
        else:
            raw = fetch_raw_for_names(names, a.email, per_name=a.per_name, strains=targets)
            json.dump(raw, io.open(rawp, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
            print("  [取数] 写入 {} 条 → {}".format(len(raw), rawp))

        got = None
        import shutil as _sh
        # —— 目标株补齐：若模式株在位点表里缺席，再按「株号 + 类型材料」定向取一次 ——
        #    实测 karsti 模式株 CBS 132134 的 ITS 存在于 NR_144790.1（from TYPE material），
        #    但按种名抓的前 150 条里没有它，导致树里找不到这片叶。
        for nm, tg in zip(names, targets):
            if not tg:
                continue
            ng0, _ = loci_of_strain(raw, tg, mod)
            if ng0 >= 1:
                continue
            try:
                from Bio import Entrez
                Entrez.email = a.email
                term = '"{}"[All Fields] AND "{}"[Organism] AND type_material[All Fields]'.format(
                    tg.replace(":", " "), nm)
                h = Entrez.esearch(db="nuccore", term=term, retmax=100)
                r = Entrez.read(h); h.close()
                ids = r.get("IdList") or []
                print("  [补齐] {} 命中 {} 条，取 {} 条".format(tg, r.get("Count"), len(ids)))
                for k in range(0, len(ids), 100):
                    h2 = Entrez.efetch(db="nuccore", id=",".join(ids[k:k+100]),
                                       rettype="gb", retmode="text")
                    raw.update(mod.parse_gb(h2.read())); h2.close()
                    time.sleep(1.2)
                ng1, nr1 = loci_of_strain(raw, tg, mod)
                print("         补齐后覆盖 {} 个位点（{} 条记录）".format(ng1, nr1))
                json.dump(raw, io.open(rawp, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
            except Exception as e:
                print("  [补齐] 失败：{}".format(type(e).__name__))

        # 先算两株各自覆盖的位点数，取较小值 → 这才是【可用的最高层级】
        # 用户规则：凑够 8 用 8，不够用 5，再不行用 ITS 单基因。
        avail = 8
        for tgt in targets:
            if not tgt:
                avail = min(avail, 1); continue
            ng, nr = loci_of_strain(raw, tgt, mod)
            print("  [层级判定] {:<12} 覆盖 {} 个位点（{} 条记录）".format(tgt, ng, nr))
            avail = min(avail, ng)
        print("  [层级判定] 可用最高层级 = ≥{} 位点".format(max(avail, 1)))
        tiers_use = [(g, lab, st) for (g, lab, st) in tiers if g <= max(avail, 1)] or tiers[-1:]
        for (min_g, lab, strength) in tiers_use:
            # ⚠️ 每级必须用【独立子目录】：
            #    Tier A 若一株都没凑齐，会写出一个空的 seqs.json；
            #    下一级看到「缓存已存在」就复用这个空文件，于是永远 0 条序列。
            #    实测就是这样三级全废的。
            twd = os.path.join(wd, "tier{}".format(lab))
            os.makedirs(twd, exist_ok=True)
            rp_src = os.path.join(wd, "raw.json")
            rp_dst = os.path.join(twd, "raw.json")
            if os.path.exists(rp_src) and not os.path.exists(rp_dst):
                _sh.copyfile(rp_src, rp_dst)
            print("\n  ── Tier {}（≥{} 位点，强度：{}）｜ 工作目录 {}".format(
                lab, min_g, strength, os.path.basename(twd)))
            res = run_tier(mod, twd, min_g, targets, a.email, a.phylo_bin, a.min_cols,
                           genes=TIER_GENES.get(lab))
            if res.get("ok"):
                got = (lab, min_g, strength)
                print("  ✅ Tier {} 跑通".format(lab))
                break
            print("  ❌ Tier {} 失败（{}）：{}".format(lab, res.get("stage"), res.get("err")))
        if not got:
            rows.append(dict(pair="{} vs {}".format(na, nb), tier="—", n_loci=0,
                             strength="无法裁决", note="三级全部失败"))
        else:
            rows.append(dict(pair="{} vs {}".format(na, nb), tier=got[0], n_loci=got[1],
                             strength=got[2], note=""))
        print("\n  结果：{}".format(rows[-1]))

    outp = os.path.join(a.out, "自动裁决汇总.tsv")
    os.makedirs(a.out, exist_ok=True)
    with io.open(outp, "w", encoding="utf-8", newline="") as f:
        f.write("争议对\t命中层级\t位点数\t结论强度\t备注\n")
        for r in rows:
            f.write("{}\t{}\t{}\t{}\t{}\n".format(r["pair"], r["tier"], r["n_loci"],
                                                  r["strength"], r["note"]))
    print("\n汇总 → {}".format(outp))
    return 0


if __name__ == "__main__":
    sys.exit(main())
