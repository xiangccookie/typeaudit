#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""chain_genus.py —— 任取一个属，一条链跑完：种清单 → 分块审计 → 合并 → 补全链。

为什么要它：炭疽菌有 chain_after.py（写死了炭疽菌），换属就得改代码。
本脚本把属名做成参数，后续换属只改命令行。

用法：
    python tools/chain_genus.py --genus Diaporthe --chunk 40
    python tools/chain_genus.py --genus Diaporthe --skip-species   # 已有清单时

断点续跑：按【种名集合】判断每块是否完整，已完成的块不重做。
"""
import argparse, datetime, io, os, subprocess, sys, time

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ENV = dict(os.environ, PYTHONPATH=os.path.join(REPO, "src"),
           NCBI_EMAIL=os.environ.get("NCBI_EMAIL") or "audit@example.org",
           PYTHONUNBUFFERED="1", PYTHONIOENCODING="utf-8")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genus", required=True)
    ap.add_argument("--chunk", type=int, default=40)
    ap.add_argument("--skip-species", action="store_true")
    a = ap.parse_args()

    g, gl = a.genus, a.genus.lower()
    root = os.path.join(REPO, "out_{}".format(gl))
    os.makedirs(root, exist_ok=True)
    LOG = os.path.join(root, "_chain.log")

    def log(msg):
        line = "[{}] {}".format(time.strftime("%H:%M:%S"), msg)
        print(line, flush=True)
        with io.open(LOG, "a", encoding="utf-8", newline="\n") as f:
            f.write(line + "\n")

    log("=" * 60)
    log("{} 全链开始".format(g))

    # ① 种清单
    spf = os.path.join(REPO, "data", "species_{}_all.txt".format(gl))
    if not a.skip_species or not os.path.exists(spf):
        log("① 取种清单（可能要几十分钟）")
        r = subprocess.run([sys.executable, os.path.join(REPO, "tools", "make_species_list.py"),
                            "--genus", g, "--out", spf],
                           cwd=REPO, env=ENV, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                           text=True, encoding="utf-8", errors="replace")
        for line in (r.stdout or "").splitlines()[-5:]:
            log("  " + line[:160])
    sps = [x.strip() for x in io.open(spf, encoding="utf-8").read().splitlines() if x.strip()]
    log("① 种清单：{} 个双名".format(len(sps)))
    if not sps:
        log("❌ 没取到种，终止"); return 2

    # ② 分块审计
    chunks = [sps[i:i + a.chunk] for i in range(0, len(sps), a.chunk)]
    log("② 分块：{} 块（每块 {}）".format(len(chunks), a.chunk))
    for idx, ch in enumerate(chunks, 1):
        d = os.path.join(root, "c{:02d}".format(idx))
        os.makedirs(d, exist_ok=True)
        tf = os.path.join(d, "_sp.txt")
        io.open(tf, "w", encoding="utf-8", newline="\n").write("\n".join(ch) + "\n")
        tbl = os.path.join(d, "{}_三列表.tsv".format(gl))
        have = set()
        if os.path.exists(tbl):
            have = {l.split("\t")[0].strip()
                    for l in io.open(tbl, encoding="utf-8-sig").read().splitlines()[1:] if l.strip()}
        missing = [x for x in ch if x not in have]
        if not missing:
            log("  跳过 c{:02d}（{} 个种已全在表里）".format(idx, len(ch))); continue
        if have:
            try: os.remove(tbl)
            except Exception: pass
        log("  开始 c{:02d}（{} 个种）".format(idx, len(ch)))
        t0 = time.time()
        r = subprocess.run([sys.executable, "-m", "typeaudit.pipeline", "--genus", g.capitalize(),
                            "--out", d, "--species-list", tf,
                            "--config", os.path.join(REPO, "config", "sources.json"),
                            "--ncbi-table", "off", "--email", ENV["NCBI_EMAIL"]],
                           cwd=os.path.join(REPO, "src"), env=ENV,
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
                           encoding="utf-8", errors="replace")
        log("  c{:02d} 结束 rc={} ｜ {:.0f}s".format(idx, r.returncode, time.time() - t0))

    # ③ 合并
    rows, order = {}, []
    cov, seen_cov = [], set()
    for idx in range(1, len(chunks) + 1):
        d = os.path.join(root, "c{:02d}".format(idx))
        tbl = os.path.join(d, "{}_三列表.tsv".format(gl))
        if os.path.exists(tbl):
            for l in io.open(tbl, encoding="utf-8-sig").read().splitlines()[1:]:
                if not l.strip(): continue
                c = l.split("\t")
                if c[0].strip() in rows: continue
                rows[c[0].strip()] = c; order.append(c[0].strip())
        cf = os.path.join(d, "{}_覆盖声明.tsv".format(gl))
        if os.path.exists(cf):
            for l in io.open(cf, encoding="utf-8-sig").read().splitlines()[1:]:
                if l.strip() and l not in seen_cov:
                    seen_cov.add(l); cov.append(l)
    merged = os.path.join(root, "{}_三列表.tsv".format(gl))
    with io.open(merged, "w", encoding="utf-8-sig", newline="") as f:
        f.write("种名\t菌株号（含问题）\t模式菌（含问题）\n")
        for sp in order:
            c = rows[sp]
            f.write("{}\t{}\t{}\n".format(c[0], c[1] if len(c) > 1 else "", c[2] if len(c) > 2 else ""))
    log("③ 合并：{} 个种 → {}".format(len(order), merged))
    if cov:
        cp = os.path.join(root, "{}_覆盖声明.tsv".format(gl))
        with io.open(cp, "w", encoding="utf-8-sig", newline="") as f:
            f.write("种名\t覆盖情况\n")
            for l in cov: f.write(l + "\n")
        log("③ 覆盖声明：{} 个种记录数超限".format(len(cov)))

    # ④ 补全链（含四方核验）
    log("④ 补全链开始（四方核验 → 反向审计 → 序列轴 → 分类库 → 汇总）")
    r = subprocess.run([sys.executable, os.path.join(REPO, "tools", "enrich.py"),
                        "--out", "out_{}".format(gl), "--genus", gl,
                        "--tree", "0", "--seq-max", "600", "--email", ENV["NCBI_EMAIL"]],
                       cwd=REPO, env=ENV, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                       text=True, encoding="utf-8", errors="replace")
    for line in (r.stdout or "").splitlines():
        if len(line.strip()) >= 6:
            log("  " + line[:170])
    log("④ 补全链结束 rc={}".format(r.returncode))
    log("🏁 {} 全链完成".format(g))
    return 0


if __name__ == "__main__":
    sys.exit(main())
