#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""run_chunks.py —— 自愈式分批审计。

【为什么要它】
    实测两个事实：
      ① 非分离启动的子进程，工具调用一结束就被杀（out_c01 跑到 18 行就没了）；
      ② 分离启动的进程能活，但约 60 分钟后被系统杀掉（out_full 活了 60.3 分钟）。
    所以：每块单独落盘 + 启动时检查"这块是否已完成" ⇒ 被杀多少次都能接着跑。

用法
    python tools/run_chunks.py           # 自动跳过已完成的块
"""
import io, os, subprocess, sys, time

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PY = sys.executable
LOG = os.path.join(REPO, "out_chunks", "_chunks.log")
ENV = dict(os.environ, PYTHONUNBUFFERED="1", PYTHONPATH=os.path.join(REPO, "src"),
           NCBI_EMAIL="audit@example.org", PYTHONIOENCODING="utf-8")


def log(msg):
    line = "[{}] {}".format(time.strftime("%H:%M:%S"), msg)
    print(line, flush=True)
    os.makedirs(os.path.dirname(LOG), exist_ok=True)
    with io.open(LOG, "a", encoding="utf-8", newline="\n") as f:
        f.write(line + "\n")


def main():
    log("=" * 60)
    log("分批审计启动")
    try:
        os.remove(os.path.join(REPO, "out_chunks", "_all_done.flag"))
    except Exception:
        pass
    chunks = sorted(f for f in os.listdir(os.path.join(REPO, "out_chunks"))
                    if f.startswith("c") and f.endswith(".txt"))
    for cf in chunks:
        name = os.path.splitext(cf)[0]
        spf = os.path.join(REPO, "out_chunks", cf)
        want = [x.strip() for x in io.open(spf, encoding="utf-8").read().splitlines() if x.strip()]
        out = os.path.join(REPO, "out_n" + name[1:])
        tbl = os.path.join(out, "colletotrichum_三列表.tsv")
        have_set = set()
        if os.path.exists(tbl):
            have_set = {l.split("\t")[0].strip()
                        for l in io.open(tbl, encoding="utf-8-sig").read().splitlines()[1:] if l.strip()}
        # ⚠️ 必须按【种名集合】判断，不能数行数 ——
        #    旧版按行数比，导致"表里有 17 行"被当成"本块 8 个种已完成"，
        #    那 8 个种其实一个都没跑（实测踩过）。
        missing = [s for s in want if s not in have_set]
        if not missing:
            log("跳过 {}（{} 个种已全在表里）".format(name, len(want)))
            continue
        if have_set:
            log("{} 缺 {} 个种，删表重跑".format(name, len(missing)))
            try:
                os.remove(tbl)
            except Exception:
                pass
        log("开始 {}（{} 个种）".format(name, len(want)))
        t0 = time.time()
        r = subprocess.run(
            [PY, "-m", "typeaudit.pipeline", "--genus", "Colletotrichum", "--out", out,
             "--species-list", spf, "--config", os.path.join(REPO, "config", "sources.json"),
             "--ncbi-table", "off", "--email", ENV["NCBI_EMAIL"]],
            cwd=os.path.join(REPO, "src"), env=ENV,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
            encoding="utf-8", errors="replace", bufsize=1)
        for line in r.stdout:
            s = line.rstrip()
            # 跳过碎片：主程序最后那段提示会被逐字吐出来（实测日志里一堆「材」「料」单字行）
            if len(s.strip()) < 6:
                continue
            if s.strip():
                log("  " + s[:180])
        log("{} 结束 rc={} ｜ {:.0f}s".format(name, r.returncode, time.time() - t0))
    log("所有分块处理完毕")
    # 写完成标记：外面的「继续跑.bat」靠它决定要不要再启动一轮
    io.open(os.path.join(REPO, "out_chunks", "_all_done.flag"), "w",
            encoding="utf-8").write(time.strftime("%Y-%m-%d %H:%M:%S"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
