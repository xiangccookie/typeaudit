#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""verify_pairs.py —— 对「跨种共号」逐条做序列比对，给出方向。

输入：out_final/核对3_跨种共号.csv
做法：每一条（号 X，声明在种 A 与种 B 下）：
  ① 取 X 自己的 ITS
  ② 取 A 的模式株 ITS、B 的模式株 ITS
  ③ 算 X 到两个模式株的距离
  ④ 判：更接近谁

结论强度（必须随结果转述）：
  单基因判「不像」= 强信号；判「像」= 弱信号（近缘复合群 ITS 分不开）
"""
import csv, io, os, re, subprocess, sys, time

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "src"))
PHYLO = "/home/yelu/miniforge3/envs/phylo/bin"


def wsl(cmd, timeout=300):
    r = subprocess.run(["wsl", "-e", "bash", "-lc", cmd], capture_output=True,
                       text=True, timeout=timeout, encoding="utf-8", errors="replace")
    return r.returncode, (r.stdout or "") + (r.stderr or "")


def main():
    from typeaudit import sources as S
    srcs = S.load(os.path.join(REPO, "config", "sources.json"))
    tx = srcs.get("taxonomy")
    email = os.environ.get("NCBI_EMAIL") or "audit@example.org"

    # 复用 seq_axis 里的取序列函数
    sys.path.insert(0, os.path.join(REPO, "tools"))
    import importlib.util
    spec = importlib.util.spec_from_file_location("_sa", os.path.join(REPO, "tools", "seq_axis.py"))
    sa = importlib.util.module_from_spec(spec); spec.loader.exec_module(sa)

    p = os.path.join(REPO, "out_final", "核对3_跨种共号.csv")
    rows = list(csv.DictReader(io.open(p, encoding="utf-8-sig")))
    print("跨种共号 {} 条".format(len(rows)))

    out_rows = []
    recs = {}
    for r in rows:
        num = (r.get("号") or "").strip()
        sps = [x.strip() for x in (r.get("涉及种") or "").split(" ｜ ") if x.strip()]
        if len(sps) < 2:
            continue
        code, _, n = num.partition(" ")
        acc, seq = sa.its_of_number(code.upper(), n, email, log=lambda *a: None)
        if not seq:
            print("  {:<14} 无 ITS，跳过".format(num)); continue
        recs["q_" + re.sub(r"[^A-Za-z0-9]", "", num)] = seq
        tinfo = []
        for sp in sps[:2]:
            full = "Colletotrichum " + sp.replace("Colletotrichum ", "")
            st = sa._type_strain(full, tx, log=lambda *a: None)
            if not st:
                tinfo.append((sp, "", None)); continue
            c2, _, n2 = st.partition(":")
            a2, s2 = sa.its_of_number(c2.upper(), n2, email, log=lambda *a: None)
            key = "t_" + sp.replace("Colletotrichum ", "")[:18]
            if s2:
                recs[key] = s2
            tinfo.append((sp, st, key if s2 else None))
        out_rows.append({"num": num, "sps": sps[:2], "tinfo": tinfo,
                         "q": "q_" + re.sub(r"[^A-Za-z0-9]", "", num)})
        print("  {:<14} {} → {}".format(num, sps[0], sps[1]))

    if not out_rows:
        print("没有可用数据"); return 2
    fa = os.path.join(REPO, "out_final", "跨种共号_序列.fa")
    # 再保险一次：任何超过 2500 bp 的序列都不是 ITS（整条基因组曾混进来）
    _bad = {k: len(v) for k, v in recs.items() if len(v) > 2500}
    if _bad:
        print("  ⚠️ 丢弃超长序列:", _bad)
        recs = {k: v for k, v in recs.items() if len(v) <= 2500}
    with io.open(fa, "w", encoding="utf-8", newline="\n") as f:
        for k, v in recs.items():
            f.write(">{}\n{}\n".format(k, v))
    aln = os.path.join(REPO, "out_final", "跨种共号_序列.aln.fa")
    def _w(x):
        return x.replace("\\", "/").replace("C:", "/mnt/c").replace("c:", "/mnt/c")
    rc, o = wsl("'{}/mafft' --auto --quiet '{}' > '{}'".format(PHYLO, _w(fa), _w(aln)))
    if rc != 0:
        print("MAFFT 失败:", o[-300:]); return 3

    # 读比对结果算距离
    seqs, name = {}, None
    for l in io.open(aln, encoding="utf-8").read().splitlines():
        if l.startswith(">"):
            name = l[1:].strip(); seqs[name] = []
        elif name:
            seqs[name].append(l.strip())
    seqs = {k: "".join(v).upper() for k, v in seqs.items()}

    def pd(a, b):
        n = d = 0
        for x, y in zip(a, b):
            if x in "-~.?" or y in "-~.?":
                continue
            n += 1
            if x != y: d += 1
        return d / n if n else None

    res = []
    for o in out_rows:
        q = o["q"]
        if q not in seqs:
            continue
        a = seqs.get(o["tinfo"][0][2]) if o["tinfo"][0][2] else None
        b = seqs.get(o["tinfo"][1][2]) if o["tinfo"][1][2] else None
        da = pd(seqs[q], a) if a else None
        db = pd(seqs[q], b) if b else None
        res.append({"号": o["num"], "种1": o["sps"][0], "种2": o["sps"][1],
                    "到种1模式株": "{:.5f}".format(da) if da is not None else "缺参考",
                    "到种2模式株": "{:.5f}".format(db) if db is not None else "缺参考",
                    "更像": ("种1（{}）".format(o["sps"][0]) if da is not None and db is not None and da < db
                             else "种2（{}）".format(o["sps"][1]) if da is not None and db is not None and db < da
                             else "打平/缺参考")})
    outp = os.path.join(REPO, "out_final", "核对3_跨种共号_序列比对结果.csv")
    with io.open(outp, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["号", "种1", "种2", "到种1模式株", "到种2模式株", "更像"])
        w.writeheader(); w.writerows(res)
    print("\n结果 → {}（{} 条）".format(os.path.relpath(outp, REPO), len(res)))
    for r in res:
        print("  {:<14} {} vs {} ｜ {} / {} ｜ 更像 {}".format(
            r["号"], r["种1"][:14], r["种2"][:14], r["到种1模式株"], r["到种2模式株"], r["更像"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
