# -*- coding: utf-8 -*-
"""
rerun_all.py —— 一键重跑：全属三列表 → 语料回归 → 守恒断言。

【为什么需要这个脚本】
    代码改动（跨库映射、检索式修正、空载守卫）会影响结论，必须整体重跑并回归。
    手敲三条命令容易漏、容易把参数接错、也不会留下可追溯的记录。
    这个脚本把三步串起来，每步计时、全过程记日志、异常时给出非零退出码。

【三步分别是什么】
    第一步  跑管线：--genus Colletotrichum 的 69 个种（种名默认从上一轮的三列表第 1 列取）
    第二步  语料回归：data/known_errors.tsv 里 20 条已知案例，工具必须判对
    第三步  守恒断言：三列表各档数量与指标文件是否自洽

【怎么用】
    python tools/rerun_all.py                     # 全跑（默认写到 out_rerun）
    python tools/rerun_all.py --limit 3           # 只跑前 3 个种，快速自检
    python tools/rerun_all.py --skip-corpus --skip-conservation   # 只跑第一步
    python tools/rerun_all.py --out out_rerun2 --email you@x.org

【退出码】
    0 = 三步全过   1 = 管线失败   2 = 语料回归失败   3 = 守恒断言失败

【注意·速度】
    接入跨库映射后，每个「无目录」号会触发 CBS 实时查询（约 4–9 秒）。
    首次跑明显变慢；查过的结果落在 .cache/cbsmap/，重跑是秒级。
    想关掉映射：把 config/sources.json 里的 "cbsmap" 那一整块删掉即可。
"""
import argparse
import io
import os
import subprocess
import sys
import time

# 脚本在 <仓库>/tools/ 下，仓库根 = 上一级
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(REPO, "src")          # 包 typeaudit 所在目录


def species_from_previous(out_dir, genus):
    """从上一轮的三列表第 1 列取种名清单。

    参数
        out_dir : 上一轮的输出目录（形如 out_audit）
        genus   : 属名，用来确认取到的是全名（如 Colletotrichum arxii）
    返回
        种名列表（字符串）。取不到就返回空列表，由调用方报错。
    """
    tsv = os.path.join(out_dir, "{}_三列表.tsv".format(genus.lower()))
    if not os.path.exists(tsv):
        return []
    lines = io.open(tsv, encoding="utf-8-sig").read().splitlines()
    out = []
    for ln in lines[1:]:                  # 跳过表头
        name = ln.split("\t")[0].strip()
        if name:
            out.append(name)
    return out


def step(title, cmd, cwd, log):
    """跑一条子命令，计时并写日志；返回 (退出码, 耗时秒)。"""
    log.write("\n" + "=" * 78 + "\n")
    log.write("[{}] {}\n".format(time.strftime("%H:%M:%S"), title))
    log.write("命令: {}\n".format(" ".join(cmd)))
    log.flush()
    print("\n" + "=" * 78)
    print("[{}] {}".format(time.strftime("%H:%M:%S"), title))
    t0 = time.time()
    # env 里带上 PYTHONPATH，子进程才能 import typeaudit
    env = dict(os.environ)
    env["PYTHONPATH"] = SRC + (os.pathsep + env["PYTHONPATH"] if env.get("PYTHONPATH") else "")
    proc = subprocess.run(cmd, cwd=cwd, env=env, capture_output=True,
                          text=True, encoding="utf-8", errors="replace")
    dt = time.time() - t0
    out = (proc.stdout or "") + (proc.stderr or "")
    log.write(out)
    log.write("\n→ 退出码 {} ｜ 耗时 {:.1f} 秒\n".format(proc.returncode, dt))
    log.flush()
    # 只把最后几行打到屏幕，避免刷屏
    tail = [x for x in out.splitlines() if x.strip()][-8:]
    for x in tail:
        print("   " + x[:150])
    print("   → 退出码 {} ｜ 耗时 {:.1f} 秒".format(proc.returncode, dt))
    return proc.returncode, dt


def main():
    ap = argparse.ArgumentParser(description="一键重跑 + 回归")
    ap.add_argument("--genus", default="Colletotrichum")
    ap.add_argument("--out", default="out_rerun", help="第一步的输出目录（相对仓库根或绝对路径）")
    ap.add_argument("--species-from", default="out_audit", help="从哪个目录的三列表取种名清单")
    ap.add_argument("--species-list", default=None, help="也可以直接给一个种名清单文件")
    ap.add_argument("--limit", type=int, default=None, help="只跑前 N 个种（快速自检用）")
    ap.add_argument("--email", default=os.environ.get("NCBI_EMAIL", "audit@example.org"))
    ap.add_argument("--config", default=os.path.join(REPO, "config", "sources.json"))
    ap.add_argument("--skip-corpus", action="store_true")
    ap.add_argument("--skip-conservation", action="store_true")
    a = ap.parse_args()

    out_dir = a.out if os.path.isabs(a.out) else os.path.join(REPO, a.out)
    os.makedirs(out_dir, exist_ok=True)
    log_path = os.path.join(out_dir, "_重跑日志.txt")
    log = io.open(log_path, "w", encoding="utf-8", newline="\n")
    log.write("重跑开始 {}\n仓库 {}\n".format(time.strftime("%Y-%m-%d %H:%M:%S"), REPO))

    # ---------- 准备种名清单 ----------
    if a.species_list:
        sp_list_path = a.species_list
        names = [x.strip() for x in io.open(a.species_list, encoding="utf-8-sig")
                 if x.strip() and not x.startswith("#")]
    else:
        names = species_from_previous(os.path.join(REPO, a.species_from), a.genus)
        if not names:
            print("❌ 取不到种名清单：{} 下没有 {}_三列表.tsv".format(a.species_from, a.genus.lower()))
            return 1
        sp_list_path = os.path.join(out_dir, "_species.txt")
        io.open(sp_list_path, "w", encoding="utf-8", newline="\n").write("\n".join(names) + "\n")
    if a.limit:
        names = names[:a.limit]
        sp_list_path = os.path.join(out_dir, "_species_limited.txt")
        io.open(sp_list_path, "w", encoding="utf-8", newline="\n").write("\n".join(names) + "\n")
    print("种名清单 {} 个 → {}".format(len(names), sp_list_path))

    results = []

    # ---------- 第一步：管线 ----------
    rc, dt = step("第一步 · 跑管线（{} 个种）".format(len(names)),
                  [sys.executable, "-m", "typeaudit.pipeline",
                   "--genus", a.genus, "--out", out_dir,
                   "--species-list", sp_list_path, "--config", a.config,
                   "--ncbi-table", "off", "--email", a.email],
                  SRC, log)
    results.append(("管线", rc, dt))
    if rc != 0:
        print("\n❌ 管线失败，后续步骤跳过（日志见 {}）".format(log_path))
        log.close()
        return 1

    # ---------- 第二步：语料回归 ----------
    if not a.skip_corpus:
        rc, dt = step("第二步 · 语料回归（data/known_errors.tsv）",
                      [sys.executable, os.path.join(REPO, "tools", "check_known_errors.py"),
                       "--out", os.path.join(REPO, "out_known_rerun"), "--email", a.email],
                      REPO, log)
        results.append(("语料回归", rc, dt))
        if rc != 0:
            print("\n❌ 语料回归未通过（这是回归，不是「新发现」——先查是不是改动破坏了判定）")

    # ---------- 第三步：守恒断言 ----------
    if not a.skip_conservation:
        rc, dt = step("第三步 · 守恒断言",
                      [sys.executable, os.path.join(REPO, "tools", "check_conservation.py"),
                       "--out", out_dir, "--genus", a.genus],
                      REPO, log)
        results.append(("守恒断言", rc, dt))

    # ---------- 汇总 ----------
    log.write("\n" + "=" * 78 + "\n汇总\n")
    print("\n" + "=" * 78)
    print("汇总")
    bad = 0
    for name, rc, dt in results:
        flag = "✅ 通过" if rc == 0 else "❌ 失败({})".format(rc)
        print("  {:<10} {}  {:.1f} 秒".format(name, flag, dt))
        log.write("  {:<10} {}  {:.1f} 秒\n".format(name, flag, dt))
        if rc != 0:
            bad = bad or rc
    log.write("\n日志结束 {}\n".format(time.strftime("%Y-%m-%d %H:%M:%S")))
    log.close()
    print("\n日志 → {}".format(log_path))
    print("产物 → {}".format(out_dir))
    return bad


if __name__ == "__main__":
    sys.exit(main())
