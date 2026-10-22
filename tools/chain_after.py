#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""chain_after.py —— 等前一个扫描结束后，自动接着跑下一轮。

【为什么需要】
    两个属同时扫会争抢 NCBI 配额（无 API key 时限 3 请求/秒），
    争抢的后果不是报错，而是**部分种静默变成"未核"** —— 这正是最难查的一类问题。
    所以让它们串行：前一个退出后才起下一个。

【做什么】
    1) 轮询，直到指定 PID 不存在（或超时）
    2) 跑炭疽菌全属（写入 out_v2/cNN）
    3) 合并成 out_v2/colletotrichum_三列表.tsv
    4) 跑补全链
"""
import io, os, subprocess, sys, time

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LOG = os.path.join(REPO, "out_v2", "_chain.log")


def log(msg):
    line = "[{}] {}".format(time.strftime("%H:%M:%S"), msg)
    print(line, flush=True)
    os.makedirs(os.path.dirname(LOG), exist_ok=True)
    with io.open(LOG, "a", encoding="utf-8", newline="\n") as f:
        f.write(line + "\n")


def alive(pid):
    r = subprocess.run(
        ["powershell", "-NoProfile", "-Command",
         "Get-Process -Id {} -ErrorAction SilentlyContinue | Select-Object -First 1 -ExpandProperty Id".format(pid)],
        capture_output=True, text=True, encoding="utf-8", errors="replace")
    return str(pid) in (r.stdout or "")


def main():
    watch = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    outroot = os.path.join(REPO, "out_v2")
    os.makedirs(outroot, exist_ok=True)
    log("=" * 60)
    log("链式任务启动 ｜ 等待 PID {} 结束".format(watch) if watch else "链式任务启动")

    # ① 等待
    waited = 0
    while watch and alive(watch):
        time.sleep(60)
        waited += 60
        if waited % 600 == 0:
            log("  仍在等待（已等 {} 分钟）".format(waited // 60))
        if waited > 8 * 3600:
            log("⚠️ 等待超过 8 小时，放弃等待直接开始")
            break
    log("等待结束（共 {} 分钟）".format(waited // 60))

    env = dict(os.environ, PYTHONPATH=os.path.join(REPO, "src"),
               NCBI_EMAIL=os.environ.get("NCBI_EMAIL") or "audit@example.org",
               PYTHONUNBUFFERED="1", PYTHONIOENCODING="utf-8")

    # ② 炭疽菌全属（分块，写到 out_v2/cNN）
    sps = [x.strip() for x in io.open(os.path.join(REPO, "data", "species_colletotrichum_all.txt"),
                                      encoding="utf-8").read().splitlines() if x.strip()]
    chunks = [sps[i:i + 40] for i in range(0, len(sps), 40)]
    log("炭疽菌 {} 个种 → {} 块".format(len(sps), len(chunks)))
    for idx, ch in enumerate(chunks, 1):
        d = os.path.join(outroot, "c{:02d}".format(idx))
        os.makedirs(d, exist_ok=True)
        tf = os.path.join(d, "_sp.txt")
        io.open(tf, "w", encoding="utf-8", newline="\n").write("\n".join(ch) + "\n")
        tbl = os.path.join(d, "colletotrichum_三列表.tsv")
        if os.path.exists(tbl):
            have = {l.split("\t")[0].strip()
                    for l in io.open(tbl, encoding="utf-8-sig").read().splitlines()[1:] if l.strip()}
            if all(x in have for x in ch):
                log("跳过 c{:02d}（已完整）".format(idx))
                continue
        log("开始 c{:02d}（{} 个种）".format(idx, len(ch)))
        t0 = time.time()
        r = subprocess.run([sys.executable, "-m", "typeaudit.pipeline", "--genus", "Colletotrichum",
                            "--out", d, "--species-list", tf,
                            "--config", os.path.join(REPO, "config", "sources.json"),
                            "--ncbi-table", "off", "--email", env["NCBI_EMAIL"]],
                           cwd=os.path.join(REPO, "src"), env=env,
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
                           encoding="utf-8", errors="replace")
        log("c{:02d} 结束 rc={} ｜ {:.0f}s".format(idx, r.returncode, time.time() - t0))

    # ③ 合并
    rows, order = {}, []
    for idx in range(1, len(chunks) + 1):
        tbl = os.path.join(outroot, "c{:02d}".format(idx), "colletotrichum_三列表.tsv")
        if not os.path.exists(tbl):
            continue
        for l in io.open(tbl, encoding="utf-8-sig").read().splitlines()[1:]:
            if not l.strip():
                continue
            c = l.split("\t")
            sp = c[0].strip()
            if sp in rows:
                continue
            rows[sp] = c
            order.append(sp)
    merged = os.path.join(outroot, "colletotrichum_三列表.tsv")
    with io.open(merged, "w", encoding="utf-8-sig", newline="") as f:
        f.write("种名\t菌株号（含问题）\t模式菌（含问题）\n")
        for sp in order:
            c = rows[sp]
            f.write("{}\t{}\t{}\n".format(c[0], c[1] if len(c) > 1 else "", c[2] if len(c) > 2 else ""))
    log("合并完成：{} 个种 → {}".format(len(order), merged))

    # ④ 补全链
    log("补全链开始")
    r = subprocess.run([sys.executable, os.path.join(REPO, "tools", "enrich.py"),
                        "--out", "out_v2", "--genus", "colletotrichum",
                        "--tree", "0", "--seq-max", "600", "--email", env["NCBI_EMAIL"]],
                       cwd=REPO, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                       text=True, encoding="utf-8", errors="replace")
    for line in (r.stdout or "").splitlines():
        if len(line.strip()) >= 6:
            log("  " + line[:170])
    log("补全链结束 rc={}".format(r.returncode))
    log("🏁 链式任务全部完成")
    return 0


if __name__ == "__main__":
    sys.exit(main())
