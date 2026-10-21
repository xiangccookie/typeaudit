#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# --- console utf8 bootstrap（Windows 控制台默认 GBK，输出 ↔ ⚠️ 🔴 会崩；详见 src/typeaudit/_console.py）
import sys as _sys
for _s in (_sys.stdout, _sys.stderr):
    try:
        _s.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass
"""run_full.py —— 全属跑批（审计 → 补全链），带日志与断点。

为什么单独一个脚本：
  · 全属 408 种约 4.5 小时（+序列轴约 2 小时），不能在一次工具调用里等；
  · 必须能【后台跑】、能【看进度】、能【中断后接着跑】。

用法
    python tools/run_full.py --out out_full --genus Colletotrichum --species-list data/species_colletotrichum_all.txt
进度
    <out>/_run.log        逐行日志
    <out>/<属>_进度.json  审计阶段的逐种进度（pipeline 自己写）
"""
import argparse, io, os, subprocess, sys, time

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PY = sys.executable
# PYTHONUNBUFFERED=1：不加的话子进程的 print 会攒成一大块才吐出来，
# 日志在跑的过程中几乎是空的（实测：跑了 1 小时只有 4 行）。
ENV = dict(os.environ, PYTHONUNBUFFERED="1", PYTHONPATH=os.path.join(REPO, "src"),
           NCBI_EMAIL=os.environ.get("NCBI_EMAIL") or "audit@example.org",
           PYTHONIOENCODING="utf-8")


def log(f, msg):
    line = "[{}] {}".format(time.strftime("%H:%M:%S"), msg)
    print(line, flush=True)
    f.write(line + "\n"); f.flush()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="out_full")
    ap.add_argument("--genus", default="Colletotrichum")
    ap.add_argument("--species-list", default="data/species_colletotrichum_all.txt")
    ap.add_argument("--tree", type=int, default=0, help="建树上限（0=不建树，留给手动并名裁决）")
    ap.add_argument("--seq-max", type=int, default=500, help="序列轴最多查多少个种")
    ap.add_argument("--skip-audit", action="store_true")
    a = ap.parse_args()

    out = a.out if os.path.isabs(a.out) else os.path.join(REPO, a.out)
    os.makedirs(out, exist_ok=True)
    spf = a.species_list if os.path.isabs(a.species_list) else os.path.join(REPO, a.species_list)
    sps = [x.strip() for x in io.open(spf, encoding="utf-8").read().splitlines() if x.strip()]
    lp = os.path.join(out, "_run.log")

    with io.open(lp, "a", encoding="utf-8", newline="\n") as f:
        log(f, "=" * 70)
        log(f, "全属跑批开始 ｜ {} 个种 ｜ 输出 {}".format(len(sps), out))
        log(f, "=" * 70)
        t0 = time.time()

        if not a.skip_audit:
            log(f, "① 审计（逐种落盘，可中断续跑）…")
            r = subprocess.run(
                [PY, "-m", "typeaudit.pipeline", "--genus", a.genus, "--out", out,
                 "--species-list", spf, "--config", os.path.join(REPO, "config", "sources.json"),
                 "--ncbi-table", "off", "--email", ENV["NCBI_EMAIL"]],
                cwd=os.path.join(REPO, "src"), env=ENV,
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                text=True, encoding="utf-8", errors="replace", bufsize=1)
            for line in r.stdout:
                line = line.rstrip()
                if line:
                    log(f, "  " + line[:200])
            log(f, "① 审计结束 rc={} ｜ {:.0f}s".format(r.returncode, time.time() - t0))
            if r.returncode != 0:
                log(f, "❌ 审计失败，停止（不继续补全）")
                return 1

        log(f, "② 补全链：四方核验 → 反向审计 → 汇总 → 序列轴 → 分类库 → 汇总")
        r2 = subprocess.run(
            [PY, os.path.join(REPO, "tools", "enrich.py"), "--out", out, "--genus",
             a.genus, "--tree", str(a.tree), "--seq-max", str(a.seq_max),
             "--email", ENV["NCBI_EMAIL"]],
            cwd=REPO, env=ENV, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, encoding="utf-8", errors="replace", bufsize=1)
        for line in r2.stdout:
            line = line.rstrip()
            if line:
                log(f, "  " + line[:200])
        log(f, "② 补全结束 rc={} ｜ 总耗时 {:.0f}s = {:.1f}h".format(
            r2.returncode, time.time() - t0, (time.time() - t0) / 3600))
        log(f, "🏁 全部完成 ｜ 产物目录 {}".format(out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
