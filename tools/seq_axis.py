#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""seq_axis.py —— 序列轴：种名与它自己的序列是否相符（按株最近邻，不建树）。

【为什么要这一步】
    号层/名字层核的是"号码与记录对不对得上"，**没有一个字涉及序列**。
    审稿人最可能问的一句是：「你说这号标错了种，但那株菌测的 ITS 到底属哪个种？」
    这一步就是回答它。

【做法（O(n)，不建树）】
    每个种取【它在「模式菌」列的第一个号】作为代表株 →
    取该号在 GenBank 里的 ITS → 全体用 MAFFT 比对 →
    按株算：d_in = 到同种其他株的最大距离；d_out = 到异种最近株的最小距离。
        d_out < d_in ⇒ 这株落到别的种里 ⇒ 并名 / 错配候选

【结论强度 —— 必须随结果一起转述】
    单基因判「不一致」= 强信号（可信）
    单基因判「一致」= 弱信号（近缘复合体内 ITS 常分不开）
    所以本步骤只用来【挑候选】，不能作为「种名正确」的证据。

用法
    python tools/seq_axis.py --out out_audit --genus colletotrichum [--locus ITS] [--max-species 200]
"""
import argparse, io, os, re, subprocess, sys, time

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "src"))

PHYLO_BIN = "/home/yelu/miniforge3/envs/phylo/bin"
CACHE = os.path.join(REPO, ".cache", "seq_axis")


def wsl(cmd, timeout=600):
    r = subprocess.run(["wsl", "-e", "bash", "-lc", cmd], capture_output=True,
                       text=True, timeout=timeout, encoding="utf-8", errors="replace")
    return r.returncode, (r.stdout or "") + (r.stderr or "")


def its_of_number(code, number, email, log=print):
    """取该号在 GenBank 里的 ITS 序列（有缓存）。返回 (acc, seq) 或 (None, None)。"""
    from Bio import Entrez
    Entrez.email = email
    tag = re.sub(r"[^A-Za-z0-9]", "_", "{}_{}".format(code, number))
    cp = os.path.join(CACHE, tag + ".json")
    if os.path.exists(cp):
        import json
        d = json.load(io.open(cp, encoding="utf-8"))
        return d.get("acc"), d.get("seq")
    term = '"{} {}"[All Fields] AND Colletotrichum[Organism]'.format(code, number)
    acc = seq = None
    try:
        h = Entrez.esearch(db="nuccore", term=term, retmax=50)
        ids = Entrez.read(h)["IdList"]; h.close()
        time.sleep(0.4)
        if ids:
            h = Entrez.efetch(db="nuccore", id=",".join(ids), rettype="gb", retmode="text")
            txt = h.read(); h.close()
            # 逐条记录找 ITS：定义行里带 internal transcribed spacer / ITS
            for rec in re.split(r"\n//\n", txt):
                if not rec.strip():
                    continue
                head = rec.split("ORIGIN")[0]
                if not re.search(r"internal transcribed spacer|ITS[12]?\b", head, re.I):
                    continue
                macc = re.search(r"^ACCESSION\s+(\S+)", head, re.M)
                morg = re.search(r'/(?:organism)="?Colletotrichum', head)
                if not (macc and morg):
                    continue
                body = rec.split("ORIGIN")[-1]
                s = re.sub(r"[^ACGTNacgtn]", "", body).upper()
                # ⚠️ 必须设上限：实测查询「(strain) 104-T」命中一条【全基因组】
                #    （5,815,657 bp），它的注释里含 "internal transcribed spacer"，
                #    于是把整条基因组当成 ITS 交了出去 → MAFFT 比对 15 分钟跑不完。
                #    ITS（含两侧保守区）一般 350–1200 bp，超过 2500 一律不是 ITS。
                if 350 <= len(s) <= 2500:
                    acc, seq = macc.group(1), s
                    break
    except Exception as e:
        log("    [!] {} {} 取序列失败: {}".format(code, number, str(e)[:60]))
    try:
        import json
        os.makedirs(CACHE, exist_ok=True)
        json.dump({"acc": acc, "seq": seq}, io.open(cp, "w", encoding="utf-8"))
    except Exception:
        pass
    return acc, seq


# ---------------------------------------------------------------- 一级：同号多名
def scan_same_number(out_dir, genus, email, log=print, max_numbers=400, sleep=0.5):
    """扫「同一保藏号在 GenBank 里被标成几个不同的种」。

    【为什么这一级最有价值】
        它有【直接证据】而不只是统计推断：同一管菌的若干条记录，organism 字段
        分别写着 Colletotrichum boninense 和 Colletotrichum karstii ——
        这不是"距离有点近"，是"数据库里自己写了两个名字"。
        实测 CBS 127535：MH877934 标 boninense，JQ005199…JQ005720 标 karstii。

    【它同时修掉了一个死标记】
        pipeline 里的 FLAG_MULTINAME「一株多名」定义了却从没被 append 过；
        本函数就是它的真实实现。
    """
    import json
    from Bio import Entrez
    Entrez.email = email
    tab = os.path.join(out_dir, "{}_三列表.tsv".format(genus.lower()))
    if not os.path.exists(tab):
        log("序列轴：没有三列表，跳过"); return []
    # 收集所有号（去重）
    nums = []
    seen = set()
    for line in io.open(tab, encoding="utf-8-sig").read().splitlines()[1:]:
        f = line.split("\t")
        sp = f[0].strip()
        for col in (1, 2):
            if len(f) <= col:
                continue
            for cell in f[col].split(" ｜ "):
                m = re.match(r"\s*([A-Za-z][A-Za-z\.\-]{0,12})\s*[:#]?\s*([\w\.\-]+)", cell)
                if not m:
                    continue
                # 过滤伪号：必须含数字，且代码不像 "SA G" 这种解析残留
                if not re.search(r"\d", m.group(2)):
                    continue
                if not re.match(r"^[A-Za-z][A-Za-z\.\-]{1,12}$", m.group(1)):
                    continue
                key = (m.group(1).upper(), m.group(2))
                if key in seen:
                    continue
                seen.add(key)
                nums.append((sp, key[0], key[1]))
    nums = nums[:max_numbers]
    log("序列轴·同号多名：检查 {} 个号".format(len(nums)))
    cp = os.path.join(CACHE, "same_number.json")
    cache = {}
    if os.path.exists(cp):
        try:
            cache = json.load(io.open(cp, encoding="utf-8"))
        except Exception:
            cache = {}
    rows = []
    for i, (sp, code, num) in enumerate(nums, 1):
        tag = "{}_{}".format(code, num)
        if tag in cache:
            orgs = cache[tag]
        else:
            orgs = []
            try:
                term = '"{} {}"'.format(code, num)
                h = Entrez.esearch(db="nuccore", term=term, retmax=40)
                ids = Entrez.read(h)["IdList"]; h.close()
                time.sleep(sleep)
                if ids:
                    h = Entrez.efetch(db="nuccore", id=",".join(ids), rettype="gb", retmode="text")
                    txt = h.read(); h.close()
                    for rec in re.split(r"\n//\n", txt):
                        mo = re.search(r'/organism="([^"]+)"', rec)
                        ma = re.search(r"^ACCESSION\s+(\S+)", rec, re.M)
                        if mo and ma:
                            orgs.append([mo.group(1).strip(), ma.group(1)])
                    time.sleep(sleep)
            except Exception as e:
                log("    [!] {} {} 查询失败: {}".format(code, num, str(e)[:50]))
                orgs = []
            cache[tag] = orgs
            if i % 20 == 0:
                log("  [{}/{}]".format(i, len(nums)))
                try:
                    os.makedirs(CACHE, exist_ok=True)
                    json.dump(cache, io.open(cp, "w", encoding="utf-8"))
                except Exception:
                    pass
        # 只保留【本属】的 organism —— 否则 "C 824" 会匹配到 Homo sapiens、
        # "SA G" 会匹配到 Pseudomonas（实测噪声）。跨属同名不算"一株多名"。
        genus_lc = genus.lower() + " "
        orgs = [o for o in orgs if o[0].lower().startswith(genus_lc)]
        names = sorted({o[0] for o in orgs})
        if len(names) >= 2:
            # 区分三种情形（不加区分会把噪声当结论）：
            #   ① 两个不同的【双名】        → 真·一株多名（错值）
            #   ② 只有一个双名 + 种下单元   → f. sp./var. 之类，非错
            #   ③ 含 Colletotrichum sp.（未命名）→ 忽略，不算第二个种
            def _base(x):
                # ⚠️ sp./cf./aff. 不是种加词 —— 不过滤会把「Colletotrichum sp. CBS 129823」
                #    当成第二个种（实测误判 2 条）。
                m = re.match(r"^(Colletotrichum [a-z][a-z\-]+)", x.strip(), re.I)
                if not m:
                    return ""
                epi = m.group(1).split()[-1].lower()
                return "" if epi in ("sp", "spp", "cf", "aff") else m.group(1).lower()
            bases = sorted({_base(x) for x in names} - {""})
            named = [x for x in names if _base(x)]
            if len(bases) >= 2:
                kind = "两个不同的种"
            elif len(bases) == 1 and len(named) >= 2:
                kind = "同一个种 + 种下单元"
            else:
                continue            # 只有未命名 sp.，不算
            rows.append({"species_claimed": sp, "code": code, "number": num,
                         "organisms": names, "accs": orgs, "kind": kind})
    try:
        os.makedirs(CACHE, exist_ok=True)
        json.dump(cache, io.open(cp, "w", encoding="utf-8"))
    except Exception:
        pass
    res = os.path.join(out_dir, "序列轴_同号多名.tsv")
    with io.open(res, "w", encoding="utf-8", newline="\n") as f:
        f.write("种\t号\tGenBank里出现的种名\t判定\t记录数\t登录号\n")
        for r in rows:
            f.write("{}\t{} {}\t{}\t{}\t{}\t{}\n".format(
                r["species_claimed"], r["code"], r["number"],
                " ｜ ".join(r["organisms"]), r["kind"], len(r["accs"]),
                ",".join(a for _, a in r["accs"][:8])))
    log("序列轴·同号多名：命中 {} 个号 → {}".format(len(rows), os.path.relpath(res, REPO)))
    for r in rows[:10]:
        log("  {:<30} {} {}  →  {}".format(
            r["species_claimed"].replace("Colletotrichum ", ""),
            r["code"], r["number"], " ｜ ".join(r["organisms"])))
    return rows


# ---------------------------------------------------------------- 二级：NCBI 分类库的一号多模式
def scan_type_collision(out_dir, genus, log=print, max_species=200):
    """反查 NCBI Taxonomy 的 type material：同一号被几个种登记为模式材料。

    【与一级的区别】
        一级看的是 GenBank【序列记录】的 organism 字段（提交者写的种名）；
        这一级看的是 NCBI【分类库】的 Type Material 字段（策展者登记的模式材料）。
        两者是不同的证据来源，可能一致也可能不一致 —— 都值得报。

    【为什么这算"一株多名"】
        一管菌不可能同时是两个种的正模。若 NCBI 分类库把同一个号登记在
        两个种名下，至少有一处登记是错的。
    """
    import sys as _sys
    _sys.path.insert(0, os.path.join(REPO, "src"))
    from typeaudit import sources as _S
    srcs = _S.load(os.path.join(REPO, "config", "sources.json"))
    tx = srcs.get("taxonomy")
    if tx is None or not hasattr(tx, "type_material_of"):
        log("分类库：没有 taxonomy 适配器，跳过")
        return []
    tab = os.path.join(out_dir, "{}_三列表.tsv".format(genus.lower()))
    if not os.path.exists(tab):
        return []
    sps = []
    for line in io.open(tab, encoding="utf-8-sig").read().splitlines()[1:]:
        sp = line.split("\t")[0].strip()
        if sp and sp not in sps:
            sps.append(sp)
    sps = sps[:max_species]
    log("分类库·一号多模式：查 {} 个种的 Type Material".format(len(sps)))

    owner = {}          # 规范化号 → {种名}
    for i, sp in enumerate(sps, 1):
        try:
            _, names = tx.type_material_of(sp)
        except Exception as e:
            log("    [!] {} 查询失败: {}".format(sp, str(e)[:50]))
            continue
        for nm in (names or []):
            m = re.match(r"^\s*([A-Za-z][A-Za-z\.\-]{0,12})\s*[:#]?\s*([\w\.\-]+)\s*$", str(nm))
            if not m:
                continue
            key = (m.group(1).upper(), m.group(2))
            owner.setdefault(key, set()).add(sp)
        if i % 20 == 0:
            log("  [{}/{}]".format(i, len(sps)))

    rows = [(k, v) for k, v in owner.items() if len(v) >= 2]
    res = os.path.join(out_dir, "分类库_一号多模式.tsv")
    with io.open(res, "w", encoding="utf-8", newline="\n") as f:
        f.write("号\t被几个种登记为正模\t种名\n")
        for (code, num), spset in sorted(rows, key=lambda x: -len(x[1])):
            f.write("{} {}\t{}\t{}\n".format(code, num, len(spset), " ｜ ".join(sorted(spset))))
    log("分类库·一号多模式：命中 {} 个号 → {}".format(len(rows), os.path.relpath(res, REPO)))
    for (code, num), spset in rows[:10]:
        log("  {} {}  →  {}".format(code, num, " ｜ ".join(sorted(x.split()[-1] for x in spset))))
    return rows


# ---------------------------------------------------------------- 三级：NCBI 分类库登记信息（预计算）
def dump_taxon_info(out_dir, genus, log=print, max_species=200):
    """把每个种的【NCBI 分类库登记信息】预计算成一张表。

    ⚠️ 为什么必须预计算：如果在网页渲染每张卡时现查 NCBI，
       250 张卡 × 2 次查询 = 首页要等好几分钟，服务会直接把连接断掉（实测）。
       所以改成"跑批时算一次、网页只读表"。
    输出：<out>/分类库_种名.tsv（种名 / taxid / 当前学名 / 定名人年份）
    """
    import sys as _sys
    _sys.path.insert(0, os.path.join(REPO, "src"))
    from typeaudit import sources as _S
    srcs = _S.load(os.path.join(REPO, "config", "sources.json"))
    tx = srcs.get("taxonomy")
    if tx is None or not hasattr(tx, "taxon_info"):
        log("分类库：没有 taxonomy 适配器，跳过"); return []
    tab = os.path.join(out_dir, "{}_三列表.tsv".format(genus.lower()))
    if not os.path.exists(tab):
        return []
    sps = []
    for line in io.open(tab, encoding="utf-8-sig").read().splitlines()[1:]:
        sp = line.split("\t")[0].strip()
        if sp and sp not in sps:
            sps.append(sp)
    sps = sps[:max_species]
    log("分类库·种名登记：查 {} 个种".format(len(sps)))
    out = os.path.join(out_dir, "分类库_种名.tsv")
    n_ok = 0
    with io.open(out, "w", encoding="utf-8", newline="\n") as f:
        f.write("种名\ttaxid\t当前学名\t定名人年份\t定名人全串\n")
        for i, sp in enumerate(sps, 1):
            try:
                d = tx.taxon_info(sp) or {}
            except Exception as e:
                log("    [!] {} 失败: {}".format(sp, str(e)[:50])); d = {}
            yr = ""
            m = re.search(r"(1[89]\d{2}|20\d{2})\s*$", (d.get("authority") or "").strip())
            if m:
                yr = m.group(1)
            if d.get("taxid"):
                n_ok += 1
            f.write("{}\t{}\t{}\t{}\t{}\n".format(
                sp, d.get("taxid", ""), d.get("scientific_name", ""), yr,
                (d.get("authority") or "").replace("&amp;", "&")[:120]))
            if i % 25 == 0:
                log("  [{}/{}]".format(i, len(sps)))
    log("分类库·种名登记：{} 个种里 {} 个查到 → {}".format(len(sps), n_ok, os.path.relpath(out, REPO)))
    return sps


# ---------------------------------------------------------------- 四级：序列自证（按模式株）
def _its_fasta(acc, email, cache_dir):
    """按登录号取 FASTA（带缓存）。"""
    import json
    from Bio import Entrez
    Entrez.email = email
    cp = os.path.join(cache_dir, "fa_{}.txt".format(acc))
    if os.path.exists(cp):
        return io.open(cp, encoding="utf-8").read().strip()
    try:
        h = Entrez.efetch(db="nuccore", id=acc, rettype="fasta", retmode="text")
        s = h.read(); h.close()
        s = "\n".join(l for l in s.splitlines() if not l.startswith(">")).replace("\n", "")
        s = re.sub(r"[^ACGTNacgtn]", "", s).upper()
        os.makedirs(cache_dir, exist_ok=True)
        io.open(cp, "w", encoding="utf-8").write(s)
        return s
    except Exception:
        return ""


def _type_strain(species, tx, log=print):
    """该种的模式株号（NCBI Taxonomy 的 type material），取第一个带数字的。"""
    try:
        _, names = tx.type_material_of(species)
    except Exception:
        return ""
    for nm in (names or []):
        m = re.match(r"^\s*([A-Za-z][A-Za-z\.\-]{0,12})\s*[:#]?\s*([\w\.\-]*\d[\w\.\-]*)\s*$", str(nm))
        if m:
            return "{}:{}".format(m.group(1).upper(), m.group(2))
    return ""


def self_check(out_dir, genus, email, log=print, max_cases=30):
    """序列自证：这管菌的序列，最像它声明的种，还是别的种？

    做法：只用【模式株】当参考 —— 不能拿随便一条同种序列当参考，那是循环论证。
    判定三档：落在自己种内 / 落到别的种 / 打平（复合群）。
    """
    import sys as _sys
    import json as _json
    _sys.path.insert(0, os.path.join(REPO, "src"))
    from typeaudit import sources as _S
    srcs = _S.load(os.path.join(REPO, "config", "sources.json"))
    tx = srcs.get("taxonomy")

    # ① 读「一株多名」清单，得到争议对
    p = os.path.join(out_dir, "序列轴_同号多名.tsv")
    if not os.path.exists(p):
        log("序列自证：没有同号多名清单，先跑 --mode same_number"); return []
    L = io.open(p, encoding="utf-8").read().splitlines()
    h = L[0].split("\t")
    cases = []
    for l in L[1:]:
        d = dict(zip(h, l.split("\t")))
        if (d.get("判定") or "") != "两个不同的种":
            continue
        sps = [x.strip() for x in (d.get("GenBank里出现的种名") or "").split("｜") if x.strip()]
        decl = (d.get("种") or "").strip()
        other = next((x for x in sps if x != decl), "")
        code_num = (d.get("号") or "").strip().split(" ", 1)
        if len(code_num) == 2 and other:
            cases.append({"species": decl, "other": other,
                          "code": code_num[0], "num": code_num[1]})
    cases = cases[:max_cases]
    log("序列自证：{} 条争议对".format(len(cases)))
    if not cases:
        return []

    # ② 参考集：涉及到的每个种的模式株 ITS
    species = sorted({c["species"] for c in cases} | {c["other"] for c in cases})
    log("  参考集：{} 个种（各取模式株）".format(len(species)))
    refs = {}
    for sp in species:
        st = _type_strain(sp, tx)
        if not st:
            log("    {} ：NCBI 没有类型材料登记".format(sp.replace("Colletotrichum ", ""))); continue
        code, num = st.split(":", 1)
        acc, seq = its_of_number(code, num, email, log=lambda *a: None)
        if seq:
            refs[sp] = (st, acc, seq)
            log("    {:<30} 模式株 {:<14} → {}（{} bp）".format(
                sp.replace("Colletotrichum ", ""), st, acc, len(seq)))
        else:
            log("    {:<30} 模式株 {:<14} → 无 ITS".format(sp.replace("Colletotrichum ", ""), st))

    # ③ 待检株：争议清单里的号
    recs = []
    for sp, (st, acc, seq) in refs.items():
        recs.append(("{}_type".format(sp.split()[-1]), seq))
    ok = 0
    for c in cases:
        acc, seq = its_of_number(c["code"], c["num"], email, log=lambda *a: None)
        if seq:
            recs.append(("{}_{}".format(c["species"].split()[-1],
                                        re.sub(r"[^A-Za-z0-9]", "", c["num"])), seq))
            ok += 1
    log("  待检株取到 {}/{} 条 ITS".format(ok, len(cases)))
    if len(recs) < 6:
        log("❌ 序列太少（{} 条），无法比较".format(len(recs))); return []

    # ④ 写 FASTA + MAFFT 比对
    fa = os.path.join(out_dir, "序列自证.fa")
    with io.open(fa, "w", encoding="utf-8", newline="\n") as f:
        for hd, sq in recs:
            f.write(">{}\n{}\n".format(hd, sq))
    aln = os.path.join(out_dir, "序列自证.aln.fa")
    def _w(pth):
        return pth.replace("\\", "/").replace("C:", "/mnt/c").replace("c:", "/mnt/c")
    rc, out = wsl("'{}/mafft' --auto --quiet '{}' > '{}'".format(
        PHYLO_BIN, _w(fa), _w(aln)), timeout=1200)
    if rc != 0 or not os.path.exists(aln):
        log("❌ MAFFT 失败：{}".format(out[-300:])); return []
    log("  比对完成 → {}".format(os.path.relpath(aln, REPO)))

    # ⑤ 最近邻检验
    res = os.path.join(out_dir, "序列自证.tsv")
    subprocess.run([sys.executable, os.path.join(REPO, "tools", "nearest_neighbor_check.py"),
                    "--aln", aln, "--out", res, "--locus", "ITS"],
                   cwd=REPO, capture_output=True, text=True, encoding="utf-8", errors="replace")
    if os.path.exists(res):
        rows = io.open(res, encoding="utf-8").read().splitlines()
        log("  结果 → {}（{} 行）".format(os.path.relpath(res, REPO), len(rows) - 1))
        for l in rows[:8]:
            log("    " + l[:150])
    return []


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="out_audit")
    ap.add_argument("--genus", default="colletotrichum")
    ap.add_argument("--locus", default="ITS")
    ap.add_argument("--max-species", type=int, default=200)
    ap.add_argument("--email", default=os.environ.get("NCBI_EMAIL") or "audit@example.org")
    ap.add_argument("--phylo-bin", default=PHYLO_BIN)
    ap.add_argument("--mode", default="same_number",
                    choices=["same_number", "type_collision", "taxon_info", "self_check", "nn"],
                    help="same_number=同号多名（便宜、证据直接）；nn=按株最近邻（要比对）")
    a = ap.parse_args()

    out_dir = a.out if os.path.isabs(a.out) else os.path.join(REPO, a.out)
    if a.mode == "self_check":
        self_check(out_dir, a.genus, a.email, print, max_cases=a.max_species)
        return 0
    if a.mode == "taxon_info":
        dump_taxon_info(out_dir, a.genus, print, max_species=a.max_species)
        return 0
    if a.mode == "type_collision":
        scan_type_collision(out_dir, a.genus, print, max_species=a.max_species)
        return 0
    if a.mode == "same_number":
        rows = scan_same_number(out_dir, a.genus, a.email, print, max_numbers=a.max_species * 6)
        return 0 if rows is not None else 2
    tab = os.path.join(out_dir, "{}_三列表.tsv".format(a.genus.lower()))
    if not os.path.exists(tab):
        print("❌ 没有三列表：{}".format(tab)); return 2

    # 每个种取「模式菌」列的第一个号
    reps = []
    for line in io.open(tab, encoding="utf-8-sig").read().splitlines()[1:]:
        f = line.split("\t")
        if len(f) < 3:
            continue
        sp = f[0].strip()
        m = re.search(r"([A-Za-z][A-Za-z\.\-]{0,12})\s*[:#]?\s*([\w\.\-]+)",
                      f[2].split(" ｜ ")[0] if f[2] else "")
        if m:
            reps.append((sp, m.group(1).upper(), m.group(2)))
    reps = reps[:a.max_species]
    print("取 {} 个种的代表株（每种 1 个模式号）".format(len(reps)))

    # 取序列
    recs, miss = [], []
    for i, (sp, code, num) in enumerate(reps, 1):
        acc, seq = its_of_number(code, num, a.email)
        if seq:
            epi = sp.split()[-1]
            recs.append((epi, code, num, acc, seq))
        else:
            miss.append((sp, code, num))
        if i % 20 == 0:
            print("  [{}/{}] 已取到 {} 条".format(i, len(reps), len(recs)), flush=True)
    print("取到序列 {} 条 ｜ 取不到 {} 条".format(len(recs), len(miss)))
    if len(recs) < 10:
        print("❌ 序列太少（{} 条），无法做组内/组间比较".format(len(recs))); return 2

    # 写 FASTA：header = <种加词>_<号>
    fa = os.path.join(out_dir, "序列轴_{}.fa".format(a.locus))
    with io.open(fa, "w", encoding="utf-8", newline="\n") as f:
        for epi, code, num, acc, seq in recs:
            f.write(">{}_{}__{}\n{}\n".format(epi, re.sub(r"[^A-Za-z0-9]", "", num), acc, seq))

    # MAFFT 比对（WSL）
    aln = os.path.join(out_dir, "序列轴_{}.aln.fa".format(a.locus))
    cmd = "'{}/mafft' --auto --quiet '{}' > '{}'".format(
        a.phylo_bin, fa.replace("\\", "/").replace("C:", "/mnt/c").replace("c:", "/mnt/c"),
        aln.replace("\\", "/").replace("C:", "/mnt/c").replace("c:", "/mnt/c"))
    rc, out = wsl(cmd, timeout=900)
    if rc != 0 or not os.path.exists(aln):
        print("❌ MAFFT 比对失败：{}".format(out[-400:])); return 3
    print("比对完成 → {}".format(os.path.relpath(aln, REPO)))

    # 最近邻检验
    res = os.path.join(out_dir, "序列轴.tsv")
    rc, out = subprocess.run(
        [sys.executable, os.path.join(REPO, "tools", "nearest_neighbor_check.py"),
         "--aln", aln, "--out", res, "--locus", a.locus],
        cwd=REPO, capture_output=True, text=True, encoding="utf-8", errors="replace").returncode, ""
    print("最近邻检验 rc = {}".format(rc))
    if os.path.exists(res):
        rows = io.open(res, encoding="utf-8").read().splitlines()
        print("结果 → {}（{} 行）".format(os.path.relpath(res, REPO), len(rows) - 1))
        # 打印候选（d_out < d_in）
        for l in rows[1:]:
            if "候选" in l or "错配" in l:
                print("  " + l[:150])
    return 0


if __name__ == "__main__":
    sys.exit(main())
