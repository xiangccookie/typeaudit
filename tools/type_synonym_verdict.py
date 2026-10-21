#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# --- console utf8 bootstrap（Windows 控制台默认 GBK，输出 ↔ ⚠️ 🔴 会崩；详见 src/typeaudit/_console.py）
import sys as _sys
for _s in (_sys.stdout, _sys.stderr):
    try:
        _s.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass
"""并名裁决：多基因建树 → 一致性判定 → 三档标记

为什么需要它
------------
「A 是不是 B 的异名」不能用单基因的『谁离谁最近』来判——本项目实测：
同一面板下 8 个基因给出 8 个不同的最近邻，单基因拓扑不可引用。
稳的用法是问「各基因是否独立地一致支持同一判断」，即**一致性原则**。

本脚本的五阶段（可断点续跑，每阶段产物落盘）
--------------------------------------------
  fetch    号段检索 + 批量取回 GenBank 记录        → raw.json
  prep     按 /culture_collection 归组 → 菌株×基因矩阵 → 取齐全集 → 写 fasta
  align    MAFFT 比对 + trimal 修剪（含列数闸门）    → panel_<gene>.trim.fa
  tree     串联 + IQ-TREE 建树                      → concat.fa / tree.treefile
  verdict  距离分布 + 拓扑检验 + 逐基因一致性 → 三档标记

用法
----
  python type_synonym_verdict.py fetch   --acc-range JX009400:JX010500 --workdir W
  python type_synonym_verdict.py prep    --workdir W --min-genes 8
  python type_synonym_verdict.py align   --workdir W
  python type_synonym_verdict.py tree    --workdir W
  python type_synonym_verdict.py verdict --workdir W --pair ICMP:18578 ICMP:18642
  # 或一次跑完：把上面五条换成 all 并给齐参数

环境依赖（关键）
----------------
* Windows 侧：Python + Biopython
* WSL 侧：conda 环境 `phylo` 内的 iqtree2 / mafft / trimal
  （本项目实测版本：iqtree2 2.4.0、mafft v7.526）
  用 --phylo-bin 指定 bin 目录，默认 /home/<user>/miniforge3/envs/phylo/bin
* NCBI 限流：每请求 ≥1.0s，且一律走本地缓存（重跑不重复请求）

分组键铁律
----------
**一律用 /culture_collection（保藏号），绝不用物种名** —— 物种名会随分类学修订漂移。
保藏号写法必须归一化（去空格/冒号/点号/连字符/国家标签 <THA>），否则
`CBS 124.22` 与 `CBS124.22`、`MFLUCC:090551` 与 `MFLUCC<THA>090551` 会被当成两株。
"""
import argparse, csv, io, itertools, json, os, re, statistics, subprocess, sys, time

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

GENES = ["ITS", "ACT", "CAL", "CHS", "GAPDH", "GS", "TUB2", "SOD2"]
GENE_PAT = [
    (r"\bactin\b|\(ACT\)", "ACT"), (r"beta-tubulin 2|\(TUB2\)|tub2", "TUB2"),
    (r"glyceraldehyde|GAPDH", "GAPDH"), (r"calmodulin|\(CAL\)", "CAL"),
    (r"chitin synthase|CHS", "CHS"), (r"superoxide dismutase|SOD2", "SOD2"),
    (r"glutamine synthetase|\(GS\)", "GS"), (r"internal transcribed|ITS", "ITS"),
]

# ---------------------------------------------------------------- 工具函数
def norm_code(s):
    """保藏号/菌株号归一：剥国家标签、去分隔符、统一大写。只用于【比较】，不用于展示。"""
    s = re.sub(r"<[^>]*>", "", s or "")
    return re.sub(r"[\s:：.\-]", "", s).upper()

def gene_of(defn):
    for pat, name in GENE_PAT:
        if re.search(pat, defn or "", re.I):
            return name
    return None

def strain_key(rec):
    """分组键：优先 ICMP，其次其它常见前缀，最后 /strain。"""
    pref = ["ICMP", "CBS", "IMI", "DAR", "BRIP", "CGMCC", "MFLUCC", "MAFF", "NBRC"]
    for p in pref:
        for c in rec.get("cc", []):
            if norm_code(c).startswith(p):
                return norm_code(c)
    st = rec.get("strain") or []
    return norm_code(st[0]) if st else "?"

def rd(h):
    d = h.read()
    return d.decode("utf-8", "replace") if isinstance(d, bytes) else d

def load_species_map(path):
    """读「物种名纠正表」：两列 TSV —— 保藏号 \t 纠正后的种名（可用 # 注释）。
    为什么需要它：GenBank 的 DEFINITION 物种名会随分类学修订漂移，
    同一株的不同基因记录可能还留着旧名（实测 ICMP:10643 有 camelliae /
    Glomerella cingulata f.sp. camelliae / C. gloeosporioides f.sp. camelliae 三种），
    按 GenBank 名分组会把不同种并到一起、把同种拆开。
    **纠正必须按【保藏号】而不是按种名** —— 按种名映射会误伤同名下的其它株。"""
    if not path or not os.path.exists(path):
        return {}
    out = {}
    for line in io.open(path, encoding="utf-8-sig"):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) >= 2 and f[0].strip():
            out[norm_code(f[0])] = f[1].strip()
    return out

def parse_gb(text):
    """GenBank flat file → {version: {defn, cc, strain, host, tm}}
    注意：^LOCUS\\s+ 必须用 \\s+（LOCUS 后是多空格对齐，\\s 只吃一个）。"""
    out = {}
    for b in re.split(r"(?m)^LOCUS\s+", text)[1:]:
        m = re.match(r"(\S+)\s+(\d+)\s+bp", b)
        if not m:
            continue
        base, L = m.group(1), int(m.group(2))
        vm = re.search(r"(?m)^VERSION\s+(\S+)", b)
        ver = vm.group(1) if vm else base
        defn = (re.search(r"DEFINITION\s+(.+?)(?=\n\s*ACCESSION)", b, re.S) or [None, ""])[1]
        get = lambda k: re.findall(r'/' + k + r'="((?:[^"]|"[^"]*")*?)"', b)
        out[ver] = dict(L=L, defn=re.sub(r"\s+", " ", defn).strip(),
                        cc=get("culture_collection"), strain=get("strain"),
                        host=get("host"), tm=re.sub(r"\s+", " ", ";".join(get("type_material"))))
    return out

def wsl(cmd, timeout=1800):
    """在 WSL 的 phylo 环境里执行 shell 命令。"""
    full = "export PATH={}:$PATH; {}".format(PHYLO_BIN, cmd)
    r = subprocess.run(["wsl", "-e", "bash", "-lc", full],
                       capture_output=True, text=True, encoding="utf-8",
                       errors="replace", timeout=timeout)
    return (r.stdout or "") + (r.stderr or "")

def to_wsl(p):
    """Windows 路径 → WSL 路径（/mnt/c/...）。"""
    p = os.path.abspath(p)
    m = re.match(r"^([A-Za-z]):[\\/](.*)$", p)
    if not m:
        return p.replace("\\", "/")
    return "/mnt/{}/{}".format(m.group(1).lower(), m.group(2).replace("\\", "/"))

# ---------------------------------------------------------------- 阶段一：fetch
def stage_fetch(a):
    from Bio import Entrez
    Entrez.email = a.email
    cache = os.path.join(a.workdir, "raw.json")
    out_meta = os.path.join(a.workdir, "fetch_meta.json")
    if os.path.exists(cache) and not a.force:
        print("[fetch] 已有缓存，跳过（--force 可重取）：{}".format(cache)); return
    os.makedirs(a.workdir, exist_ok=True)
    term = "{}[ACCN] AND Colletotrichum[Organism]".format(a.acc_range)
    h = Entrez.esearch(db="nuccore", term=term, retmax=5000)
    rec = Entrez.read(h); h.close()
    ids = rec.get("IdList", [])
    print("[fetch] 检索式 {}  → count={} 取到 {}".format(term, rec.get("Count"), len(ids)))
    raw = {}
    B = 150                       # 实测：单次 150 条 gb 稳定；22 条以上不能用 gbc（返回 XML）
    for k in range(0, len(ids), B):
        part = ids[k:k + B]
        for attempt in range(3):
            try:
                txt = rd(Entrez.efetch(db="nuccore", id=",".join(part),
                                       rettype="gb", retmode="text"))
                break
            except Exception as e:
                print("   批 {} 第 {} 次失败：{}".format(k, attempt + 1, str(e)[:80]))
                time.sleep(3)
        else:
            continue
        raw.update(parse_gb(txt))
        print("   批 {:>4}-{:<4} 累计 {}".format(k, k + len(part), len(raw)))
        time.sleep(a.sleep)
    json.dump(raw, io.open(cache, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    json.dump({"term": term, "count": rec.get("Count"), "n_records": len(raw)},
              io.open(out_meta, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    print("[fetch] 完成：{} 条 → {}".format(len(raw), cache))

# ---------------------------------------------------------------- 阶段二：prep
def stage_prep(a):
    raw = json.load(io.open(os.path.join(a.workdir, "raw.json"), encoding="utf-8"))
    cell, sp_of = {}, {}
    for ver, r in raw.items():
        g = gene_of(r.get("defn"))
        if not g:
            continue
        k = strain_key(r)
        cell.setdefault(k, {})
        if g not in cell[k] or r.get("L", 0) > cell[k][g][1]:
            cell[k][g] = (ver, r.get("L", 0))
        m = re.search(r"Colletotrichum\s+([a-z\-]+)", r.get("defn", ""))
        sp_of.setdefault(k, m.group(1) if m else "sp")
    # —— 物种名纠正（可选，按保藏号）——
    smap = load_species_map(getattr(a, "species_map", None))
    if smap:
        n_fix = 0
        for k in list(sp_of):
            if k in smap:
                print("   改标 {}：{} → {}".format(k, sp_of[k], smap[k]))
                sp_of[k] = smap[k]
                n_fix += 1
        print("[prep] 物种名纠正表生效：{} 株被改标".format(n_fix))
    full = [k for k in cell if len(cell[k]) >= a.min_genes]
    print("[prep] 菌株 {} 株；≥{} 基因齐全 {} 株".format(len(cell), a.min_genes, len(full)))
    # —— 取序列（缓存；这一步容易漏：矩阵里存的是【登记号】，fasta 里必须是【序列】）——
    accs = sorted({cell[k][g][0] for k in full for g in cell[k]})
    scache = os.path.join(a.workdir, "seqs.json")
    if os.path.exists(scache) and not a.force:
        seqs = json.load(io.open(scache, encoding="utf-8"))
        print("[prep] 序列缓存命中：{} 条".format(len(seqs)))
    else:
        from Bio import Entrez
        Entrez.email = a.email
        seqs = {}
        for i in range(0, len(accs), 160):
            chunk = accs[i:i + 160]
            for attempt in range(3):
                try:
                    txt = rd(Entrez.efetch(db="nuccore", id=",".join(chunk),
                                           rettype="fasta", retmode="text"))
                    break
                except Exception as e:
                    print("   取序列批 {} 第 {} 次失败：{}".format(i, attempt + 1, str(e)[:70]))
                    time.sleep(3)
            else:
                continue
            cur = None
            for ln in txt.splitlines():
                if ln.startswith(">"):
                    cur = ln[1:].split()[0]
                    seqs[cur] = ""
                elif cur:
                    seqs[cur] += ln.strip().upper()
            print("   取序列 {:>4}-{:<4} 累计 {}".format(i, i + len(chunk), len(seqs)))
            time.sleep(a.sleep)
        json.dump(seqs, io.open(scache, "w", encoding="utf-8"), ensure_ascii=False)
        print("[prep] 序列已缓存：{}".format(scache))

    # 写 fasta（标签 = C<种前4字母>_<保藏号>；内容是【序列】）
    os.makedirs(a.workdir, exist_ok=True)
    for g in GENES:
        recs = []
        for k in sorted(full, key=lambda x: (sp_of[x], x)):
            if g in cell[k]:
                acc = cell[k][g][0]
                sq = seqs.get(acc) or seqs.get(acc.split(".")[0])
                if not sq:
                    print("      ⚠️ {} 的 {} 取不到序列（{}）".format(k, g, acc))
                    continue
                recs.append(("C{}_{}".format(sp_of[k][:4], k), sq))
        if not recs:
            continue
        p = os.path.join(a.workdir, "panel_{}.fa".format(g))
        io.open(p, "w", encoding="ascii", newline="\n").write(
            "".join(">{}\n{}\n".format(n, v) for n, v in recs))
        print("   {:<7}{:>4} 条 → {}".format(g, len(recs), os.path.basename(p)))
    json.dump({"cell": {k: {g: v[0] for g, v in vv.items()} for k, vv in cell.items()},
               "species": sp_of, "full": full},
              io.open(os.path.join(a.workdir, "panel_index.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    print("[prep] 完成")

# ---------------------------------------------------------------- 阶段三：align
def stage_align(a):
    """WSL 只负责跑 mafft / trimal；**列数一律用 Python 的 AlignIO 测**。
    教训：用 `seqkit stats | awk '$8'` 经两层引号传递会取错字段，
    而 `awk NR==2` 数行在序列折行时更会得到假值——测量工具本身必须先被校验。"""
    from Bio import AlignIO
    need = [g for g in GENES if os.path.exists(os.path.join(a.workdir, "panel_{}.fa".format(g)))]
    w = to_wsl(a.workdir)
    bad = []
    for g in need:
        cmd = ("cd '{w}' && mafft --auto --quiet panel_{g}.fa > panel_{g}.aln.fa 2>/dev/null && "
               "trimal -in panel_{g}.aln.fa -out panel_{g}.trim.fa -gt 0.8 2>/dev/null && echo DONE"
               ).format(w=w, g=g)
        out = wsl(cmd, timeout=900)
        ca = ct = 0
        try:
            ca = AlignIO.read(os.path.join(a.workdir, "panel_{}.aln.fa".format(g)),
                              "fasta").get_alignment_length()
        except Exception as e:
            print("   {:<7}比对结果读取失败：{}".format(g, str(e)[:60]))
        try:
            ct = AlignIO.read(os.path.join(a.workdir, "panel_{}.trim.fa".format(g)),
                              "fasta").get_alignment_length()
        except Exception as e:
            print("   {:<7}修剪结果读取失败：{}".format(g, str(e)[:60]))
        ok = ct >= a.min_cols
        if not ok:
            bad.append((g, ct))
        print("   {:<7}比对={:<5} 修剪={:<5} {}".format(g, ca, ct, "✅" if ok else "❌ 列数不足"))
    if bad:
        print("\n   ⚠️ 以下基因列数低于闸门 {}：{}".format(a.min_cols, bad))
        print("      若列数塌缩，通常是有序列过短导致比对错位——检查各序列长度后再决定是否剔除。")
    print("[align] 完成（闸门 {} 列）".format(a.min_cols))

# ---------------------------------------------------------------- 阶段四：tree
def stage_tree(a):
    from Bio import SeqIO
    w = to_wsl(a.workdir)
    # 串联 —— 必须取【各基因的序列交集】
    # ⚠️ 原来这里是 `order = order or ids`，即拿第一个基因的名单当全体：
    #    · 若第一个基因只有 1 条，串联就只有 1 条 → IQ-TREE 报「至少 3 条序列」；
    #    · 若各基因成员不同，串联长度参差 → IQ-TREE 报「序列字符数不足」。
    #    实测这两条都在 celtidis vs karsti 上撞到了，所以改成求交集。
    per_gene = {}
    # 层级可限定基因集合：Tier C「只用 ITS」时只取 ITS 面板，
    # 否则会把 ITS/ACT/CAL/... 一起求交集 —— 实测交集会塌成 1 条，树建不出来。
    use_genes = getattr(a, "genes", None) or GENES
    for g in use_genes:
        p = os.path.join(a.workdir, "panel_{}.trim.fa".format(g))
        if not os.path.exists(p):
            continue
        recs = list(SeqIO.parse(p, "fasta"))
        if len(recs) >= 2:                 # 单条序列的基因没有比对信息，不参与串联
            per_gene[g] = {r.id: str(r.seq).upper() for r in recs}
    if not per_gene:
        print("   ❌ 没有任何 panel_*.trim.fa（或都只有 1 条）—— prep 阶段没有菌株凑齐位点，无法建树")
        print("      多半是位点太散。可降低 --min-genes，或按模式株号补抓该菌株的记录。")
        return
    ids_sets = [set(v) for v in per_gene.values()]
    common = set.intersection(*ids_sets)
    print("[tree] 参与串联的基因 {} 个：{}".format(len(per_gene), list(per_gene)))
    print("[tree] 各基因共有菌株 {} 条（交集）；各基因条数：{}".format(
        len(common), {g: len(v) for g, v in per_gene.items()}))
    if len(common) < 3:
        print("   ❌ 交集只有 {} 条，低于建树下限 3 条 → 拒绝建树".format(len(common)))
        print("      原因：各基因的菌株集合几乎不重叠（每株只测了一两个位点）。")
        print("      处置：降级到更少的位点（如只用 ITS），或按菌株补抓齐各位点。")
        return
    order = sorted(common)
    concat = {n: "".join(per_gene[g][n] for g in per_gene) for n in order}
    cp = os.path.join(a.workdir, "concat.fa")
    io.open(cp, "w", encoding="ascii", newline="\n").write(
        "".join(">{}\n{}\n".format(n, s) for n, s in concat.items()))
    total = len(next(iter(concat.values())))
    print("[tree] 串联 {} 条 × {} bp".format(len(concat), total))
    if total < a.min_cols * 2:
        print("   ❌ 闸门：串联仅 {} bp，远低于预期（各基因列数需先过 align 阶段闸门），拒绝建树".format(total))
        return
    out = wsl("cd '{}' && iqtree2 -s concat.fa -m MFP -B 1000 -T AUTO -pre tree -redo "
              "> iqtree.stdout 2>&1; echo EXIT=$?; "
              "grep -E 'Number of parsimony informative sites|Best-fit model according' tree.iqtree"
              .format(w), timeout=3600)
    print(out.strip())
    tf = os.path.join(a.workdir, "tree.treefile")
    if not os.path.exists(tf):
        print("   ❌ 没有生成 tree.treefile —— IQ-TREE 未成功建树（看上面的 iqtree.stdout）")
        return False
    print("[tree] 完成 → {}".format(tf))
    return True

# ---------------------------------------------------------------- 阶段五：verdict
def stage_verdict(a):
    from Bio import Phylo
    core = [c for c in (a.pair or []) if c]
    tree = Phylo.read(os.path.join(a.workdir, "tree.treefile"), "newick")
    leaves = [t.name for t in tree.get_terminals()]
    sp2 = {}
    for l in leaves:
        sp2.setdefault(l.split("_")[0], []).append(l)

    def find(code):
        """按保藏号的【数字部分】匹配叶标签。
        叶标签形如 Csiam_18578（归一后 CXIAM18578），而目标写 ICMP:18578（归一后 ICMP18578）
        —— 前缀对不上，所以必须比数字。"""
        n = norm_code(code)
        num = re.sub(r"^[A-Z]+", "", n)          # ICMP18578 → 18578
        num = re.sub(r"^0+", "", num)            # 去前导零：CBS 0123 → 123
        if not num:
            return None
        # 第一轮：数字段包含（叶标签形如 Ckars_CBS132134）
        hits = [l for l in leaves if num and num in re.sub(r"^0+", "", norm_code(l))]
        # 第二轮：整串包含（叶标签可能直接用了菌株名）
        if not hits:
            hits = [l for l in leaves if n in norm_code(l)]
        if len(hits) > 1:
            print("   ⚠️ {} 匹配到多片叶：{}（取第一个）".format(code, hits))
        return hits[0] if hits else None
    A, B = (find(core[0]), find(core[1])) if len(core) >= 2 else (None, None)
    print("[verdict] 目标：{} → {} ｜ {} → {}".format(core[0] if core else "?", A,
                                                   core[1] if len(core) > 1 else "?", B))
    if not A or not B:
        print("   ❌ 目标菌株未出现在树中，无法裁决"); return False

    # 距离分布
    intra = {}
    for sp, ls in sp2.items():
        if len(ls) >= 2:
            intra[sp] = [tree.distance(x, y) for x, y in itertools.combinations(ls, 2)]
    inter = [tree.distance(x, y) for s1, s2 in itertools.combinations(sorted(sp2), 2)
             for x in sp2[s1] for y in sp2[s2]]
    dAB = tree.distance(A, B)
    spA = A.split("_")[0]
    si = intra.get(spA, [])

    # 拓扑：A、B 的最小共同祖先是否只含同种
    cl = tree.common_ancestor(A, B)
    desc = [t.name for t in cl.get_terminals()]
    alien = [d for d in desc if not d.startswith(spA)]

    # 判定
    tag, why = "并名-4 未裁·缺数据", ""
    if len(si) < 5:
        why = "{} 种内仅 {} 对，不足以定义种内分布".format(spA, len(si))
    elif dAB <= max(si) and not alien:
        tag = "并名-2 异株·遗传支持"
        why = ("争议对 {:.5f} ≤ {} 种内 max {:.5f}；最小共同祖先含 {} 叶、无外种混入，bootstrap={}"
               .format(dAB, spA, max(si), len(desc), cl.confidence))
    elif dAB > max(si):
        tag = "并名-3 异株·遗传不支持"
        why = "争议对 {:.5f} > {} 种内 max {:.5f}".format(dAB, spA, max(si))
    else:
        why = "最小共同祖先混入外种：{}".format(",".join(alien[:4]))

    print("\n" + "=" * 92)
    print("【裁决】{}".format(tag))
    print("   {}".format(why))
    print("=" * 92)
    print("   争议对 dist = {:.5f}".format(dAB))
    print("   {} 种内 {} 对：min={:.5f} 中位={:.5f} max={:.5f}".format(
        spA, len(si), min(si), statistics.median(si), max(si)))
    if inter:
        print("   全体种间 {} 对：min={:.5f} 中位={:.5f} max={:.5f}".format(
            len(inter), min(inter), statistics.median(inter), max(inter)))
    print("   多株物种单系性：")
    mono = 0
    for sp, ls in sorted(sp2.items()):
        if len(ls) < 2:
            continue
        c = tree.common_ancestor(ls)
        ds = [t.name for t in c.get_terminals()]
        ex = [d for d in ds if not d.startswith(sp)]
        mono += (not ex)
        print("      {:<8}n={:<3}{}  bootstrap={}".format(
            sp, len(ls), "✅单系" if not ex else "❌非单系(混入{})".format(len(ex)), c.confidence))
    print("   ⇒ 阴性对照：{}/{} 个多株物种单系（方法可分辨）".format(
        mono, sum(1 for ls in sp2.values() if len(ls) >= 2)))

    out = os.path.join(a.workdir, "verdict_{}.tsv".format(
        re.sub(r"\W+", "", core[0] + "_" + core[1]) if len(core) >= 2 else "result"))
    io.open(out, "w", encoding="utf-8-sig", newline="").write(
        "项\t值\n目标A\t{}\n目标B\t{}\n争议对距离\t{:.5f}\n{}种内max\t{:.5f}\n判定\t{}\n依据\t{}\n"
        .format(core[0], core[1], dAB, spA, max(si) if si else float("nan"), tag, why))
    print("\n   裁决表 → {}".format(out))

# ---------------------------------------------------------------- main
# ---------------------------------------------------------------- 阶段六：plot
    return True


def stage_plot(a):
    """出图（**一律用英文标签**——本机 CJK 字体不保证可用，中文会变方框）。
       tree_verdict.{png,pdf}  ：标注树图（目标种单系群高亮 + 目标两株标记）
       distance_dist.{png,pdf} ：种内 vs 种间距离分布 + 争议对位置"""
    from Bio import Phylo
    import itertools as it, statistics as stt
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    figdir = os.path.join(a.workdir, "figures")
    os.makedirs(figdir, exist_ok=True)
    tree = Phylo.read(os.path.join(a.workdir, "tree.treefile"), "newick")
    tree.ladderize()
    leaves = [t.name for t in tree.get_terminals()]
    sp2 = {}
    for l in leaves:
        sp2.setdefault(l.split("_")[0], []).append(l)

    core = [c for c in (a.pair or []) if c]
    def find(code):
        n = norm_code(code); num = re.sub(r"^[A-Z]+", "", n)
        hits = [l for l in leaves if num and num in norm_code(l)]
        return hits[0] if hits else None
    A, B = (find(core[0]), find(core[1])) if len(core) >= 2 else (None, None)
    spA = A.split("_")[0] if A else None
    group = set(sp2.get(spA, []))
    dAB = tree.distance(A, B) if (A and B) else None

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 6.2})

    # ---- 图一：标注树图 ----
    fig, ax = plt.subplots(figsize=(7.4, 11))
    def lab(t):
        n = t.name or ""
        return ("■ " + n) if n in (A, B) else (("▪ " + n) if n in group else n)
    def col(x):
        x = str(x or "")
        return "#D55E00" if any(g in x for g in group) else "#444444"
    Phylo.draw(tree, axes=ax, do_show=False, label_func=lab, label_colors=col)
    ys = [t.get_position()[1] for t in ax.texts
          if t.get_text().strip() and any(g in t.get_text() for g in group)]
    if ys:
        y0, y1 = min(ys) - 0.6, max(ys) + 0.6
        ax.axhspan(y0, y1, color="#D55E00", alpha=0.09, zorder=0)
        cl = tree.common_ancestor(sp2[spA])
        ax.annotate("{} clade (bootstrap {})".format(spA, cl.confidence),
                    xy=(ax.get_xlim()[0], (y0 + y1) / 2),
                    xytext=(ax.get_xlim()[1] * 0.40, (y0 + y1) / 2),
                    fontsize=6.4, color="#B34700", va="center",
                    arrowprops=dict(arrowstyle="->", color="#B34700", lw=0.8))
    ax.set_title("{} loci concatenated tree ({} strains / {} species)".format(
        len([g for g in GENES if os.path.exists(os.path.join(a.workdir, "panel_{}.trim.fa".format(g)))]),
        len(leaves), len(sp2)), fontsize=7.6, loc="left", weight="bold")
    ax.set_xlabel("Patristic distance (substitutions per site)", fontsize=7)
    ax.tick_params(axis="y", labelsize=5.6)
    for sp_ in ("top", "right"):
        ax.spines[sp_].set_visible(False)
    if A and B:
        ax.legend(handles=[
            Line2D([], [], color="#D55E00", lw=2, label="{} clade".format(spA)),
            Line2D([], [], marker="s", ls="", color="#D55E00", label="■ " + core[0]),
            Line2D([], [], marker="s", ls="", color="#D55E00", alpha=0.55, label="■ " + core[1]),
        ], loc="lower right", fontsize=6, frameon=False)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(figdir, "tree_verdict." + ext), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("[plot] 树图 → figures/tree_verdict.png / .pdf")

    # ---- 图二：种内 vs 种间分布 ----
    intra = {sp: [tree.distance(x, y) for x, y in it.combinations(ls, 2)]
             for sp, ls in sp2.items() if len(ls) >= 2}
    inter = [tree.distance(x, y) for s1, s2 in it.combinations(sorted(sp2), 2)
             for x in sp2[s1] for y in sp2[s2]]
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.2), gridspec_kw={"width_ratios": [1.35, 1]})
    ax = axes[0]
    order = [sp for sp, _ in sorted(intra.items(), key=lambda x: -stt.median(x[1]))]
    bp = ax.boxplot([intra[sp] for sp in order], widths=0.62, patch_artist=True,
                    medianprops=dict(color="black", lw=1.0),
                    flierprops=dict(marker="o", ms=1.6, mfc="#888888", mec="none"))
    for i, box in enumerate(bp["boxes"]):
        box.set(facecolor="#D55E00" if order[i] == spA else "#BBBBBB", alpha=0.75, lw=0.6)
    if dAB:
        ax.axhline(dAB, color="#0072B2", lw=1.4, ls="--", zorder=3)
        ax.text(len(order) + 0.45, dAB, "contested\n({:.4f})".format(dAB),
                color="#0072B2", fontsize=6.6, va="center")
    ax.set_xticks(range(1, len(order) + 1))
    ax.set_xticklabels(order, rotation=55, ha="right", fontsize=6.4)
    ax.set_ylabel("Patristic distance", fontsize=7)
    ax.set_title("a  Intra-species distance", fontsize=8, loc="left", weight="bold")
    ax.grid(axis="y", lw=0.3, alpha=0.4); ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    ax2 = axes[1]
    alli = [d for v in intra.values() for d in v]
    ax2.hist(inter, bins=40, color="#999999", alpha=0.85, label="Inter-species (n={})".format(len(inter)))
    ax2.hist(alli, bins=40, color="#D55E00", alpha=0.85, label="Intra-species (n={})".format(len(alli)))
    if dAB:
        ax2.axvline(dAB, color="#0072B2", lw=1.6, ls="--", label="Contested pair")
    ax2.set_xlabel("Patristic distance", fontsize=7); ax2.set_ylabel("Pairs", fontsize=7)
    ax2.set_title("b  Intra vs inter", fontsize=8, loc="left", weight="bold")
    ax2.legend(fontsize=6, frameon=False)
    ax2.grid(axis="y", lw=0.3, alpha=0.4); ax2.spines["top"].set_visible(False); ax2.spines["right"].set_visible(False)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(figdir, "distance_dist." + ext), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("[plot] 分布图 → figures/distance_dist.png / .pdf")

PHYLO_BIN = "/home/yelu/miniforge3/envs/phylo/bin"
STAGES = dict(fetch=stage_fetch, prep=stage_prep, align=stage_align,
              tree=stage_tree, verdict=stage_verdict, plot=stage_plot)

def main():
    global PHYLO_BIN          # 必须在任何使用之前声明
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=list(STAGES) + ["all"])
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--acc-range", default="JX009400:JX010500")
    ap.add_argument("--pair", nargs=2, default=None, help="两个保藏号，如 ICMP:18578 ICMP:18642")
    ap.add_argument("--min-genes", type=int, default=8)
    ap.add_argument("--min-cols", type=int, default=200, help="每基因修剪后最小列数（闸门）")
    ap.add_argument("--species-map", default=None,
                    help="物种名纠正表（两列 TSV：保藏号\\t纠正后种名）；不按种名映射，避免误伤同名的其它株")
    ap.add_argument("--email", default=os.environ.get("NCBI_EMAIL", "researcher@example.org"))
    ap.add_argument("--sleep", type=float, default=1.2, help="NCBI 请求间隔秒")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--phylo-bin", default=PHYLO_BIN)
    a = ap.parse_args()
    PHYLO_BIN = a.phylo_bin
    os.makedirs(a.workdir, exist_ok=True)
    todo = list(STAGES) if a.stage == "all" else [a.stage]
    for s in todo:
        print("\n" + "#" * 92)
        print("# 阶段：{}".format(s))
        print("#" * 92)
        STAGES[s](a)

if __name__ == "__main__":
    main()
