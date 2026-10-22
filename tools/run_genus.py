#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""run_genus.py —— 通用属扫描驱动（分块 + 可续跑）。

为什么不复用 run_chunks.py：那个脚本把属名与目录名写死成炭疽菌了。
本脚本把属名、种清单、块大小都做成参数，方便换属（用户要求：炭疽菌跑通后换一个属）。

续跑判定：按【种名集合】比 —— 表里已有的种不再跑，缺的才补。
⚠️ 不能按行数比：曾经用行数，结果"表里 17 行"被当成"本块 8 个种已完成"，
   那 8 个种一个都没跑（实测踩过）。
"""
import argparse, io, os, subprocess, sys, time

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ENV = dict(os.environ, PYTHONUNBUFFERED="1", PYTHONPATH=os.path.join(REPO, "src"),
           NCBI_EMAIL=os.environ.get("NCBI_EMAIL") or "audit@example.org",
           PYTHONIOENCODING="utf-8")


def log(msg, logf):
    line = "[{}] {}".format(time.strftime("%H:%M:%S"), msg)
    print(line, flush=True)
    with io.open(logf, "a", encoding="utf-8", newline="\n") as f:
        f.write(line + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genus", required=True, help="属名，如 Trichoderma")
    ap.add_argument("--species", required=True, help="种清单文件")
    ap.add_argument("--chunk", type=int, default=40)
    ap.add_argument("--start", type=int, default=1, help="从第几块开始（人工纠偏用）")
    a = ap.parse_args()

    g = a.genus
    gl = g.lower()
    spf = a.species if os.path.isabs(a.species) else os.path.join(REPO, a.species)
    sps = [x.strip() for x in io.open(spf, encoding="utf-8").read().splitlines() if x.strip()]
    chunks = [sps[i:i + a.chunk] for i in range(0, len(sps), a.chunk)]
    outroot = os.path.join(REPO, "out_{}".format(gl))
    os.makedirs(outroot, exist_ok=True)
    logf = os.path.join(outroot, "_chunks.log")

    log("=" * 60, logf)
    log("{} 扫描开始 ｜ {} 个种 ｜ {} 块（每块 {}）".format(g, len(sps), len(chunks), a.chunk), logf)

    for idx, ch in enumerate(chunks, 1):
        if idx < a.start:
            continue
        out = os.path.join(outroot, "c{:02d}".format(idx))
        tbl = os.path.join(out, "{}_三列表.tsv".format(gl))
        have = set()
        if os.path.exists(tbl):
            for l in io.open(tbl, encoding="utf-8-sig").read().splitlines()[1:]:
                if l.strip():
                    have.add(l.split("\t")[0].strip())
        missing = [s for s in ch if s not in have]
        if not missing:
            log("跳过 c{:02d}（{} 个种已全在表里）".format(idx, len(ch)), logf)
            continue
        if have:
            log("c{:02d} 缺 {} 个种，删表重跑".format(idx, len(missing)), logf)
            try:
                os.remove(tbl)
            except Exception:
                pass
        tf = os.path.join(out, "_sp.txt")
        os.makedirs(out, exist_ok=True)
        io.open(tf, "w", encoding="utf-8", newline="\n").write("\n".join(ch) + "\n")
        log("开始 c{:02d}（{} 个种）".format(idx, len(ch)), logf)
        t0 = time.time()
        r = subprocess.run(
            [sys.executable, "-m", "typeaudit.pipeline", "--genus", g, "--out", out,
             "--species-list", tf, "--config", os.path.join(REPO, "config", "sources.json"),
             "--ncbi-table", "off", "--email", ENV["NCBI_EMAIL"]],
            cwd=os.path.join(REPO, "src"), env=ENV,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
            encoding="utf-8", errors="replace", bufsize=1)
        for line in r.stdout:
            s = line.rstrip()
            if len(s.strip()) >= 6:
                log("  " + s[:180], logf)
        log("c{:02d} 结束 rc={} ｜ {:.0f}s".format(idx, r.returncode, time.time() - t0), logf)
    log("{} 全部块处理完毕".format(g), logf)
    return 0


if __name__ == "__main__":
    sys.exit(main())
