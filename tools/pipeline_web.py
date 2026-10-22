#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""网页版 · 属级类型材料审计 + 审稿系统（零依赖：只用标准库 http.server）

启动
----
    python tools/pipeline_web.py               # 默认 http://127.0.0.1:8765
    python tools/pipeline_web.py --port 9000   # 换端口

为什么用 http.server 而不是 Flask
---------------------------------
审稿人不该为了看一页表面去装 web 框架。标准库够用：单进程、局域网内单人用。

两个区块
--------
① 跑流水线：填属名 → 开始 → 实时日志（轮询 /api/log）→ 出三列表
② 审稿：把三列表里的每个「候选问题」列成卡片，逐条裁定
        确认 / 驳回 / 存疑 + 备注，**追加写 review_log.jsonl 留痕**（append-only）

为什么必须有审稿区（项目一贯原则）
----------------------------------
机械检查的候选**假阳性高**，只可用于缩小候选集：
本项目实测 J-07「一号两主」12 条候选里绝大多数是块污染、J-08「位数孤例」14 条几乎全是假阳性。
所以：**机器只出候选，人裁定；裁定结果永不静默删除，一律留痕。**

留痕格式（每行一条，append-only）
--------------------------------
    {"ts": "...", "genus": "...", "species": "...", "item": "CBS 147945",
     "system_flag": "目录记为别种", "verdict": "确认|驳回|存疑",
     "reviewer": "...", "note": "..."}
"""
import argparse
import io
import json
import os
import queue
import re
import subprocess
import sys
import threading
import time
from collections import deque
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import parse_qs, urlparse

try:
    sys.stdout.reconfigure(encoding="utf-8")
    sys.stderr.reconfigure(encoding="utf-8")
except Exception:
    pass

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
SRC = os.path.join(REPO, "src")
if os.path.isdir(SRC) and SRC not in sys.path:
    sys.path.insert(0, SRC)

LOG = deque(maxlen=4000)          # 日志环形缓冲
LOCK = threading.Lock()
STATE = {"running": False, "stop": False, "genus": "", "out_dir": "",
         "rc": None, "t0": None}


def load_refs():
    """读 find_refdata.py 写的 registry（~/.dsh/typeaudit_refs.json），
    用于**自动填**三个参考数据框——省得用户每次手打长路径（还容易丢 `C:`）。"""
    p = os.path.join(os.path.expanduser("~"), ".dsh", "typeaudit_refs.json")
    try:
        if os.path.exists(p):
            import json as _json
            d = _json.load(io.open(p, encoding="utf-8"))
            return {k: (v or {}).get("path") for k, v in d.items() if (v or {}).get("path")}
    except Exception:
        pass
    return {}


REFS = load_refs()


def _mtime(path):
    """表的生成时间——用来提醒"这是上次的结果"，别把旧数据当本次的。"""
    try:
        if path and os.path.exists(path):
            return time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(os.path.getmtime(path)))
    except Exception:
        pass
    return None


def deps():
    """依赖体检：网页顶部要显示的东西。"""
    out = {"python": sys.executable}
    try:
        import tkinter as _tk
        out["tkinter"] = "OK {}".format(_tk.TkVersion)
    except Exception:
        out["tkinter"] = "缺失"
    try:
        import Bio
        out["biopython"] = "OK {}".format(Bio.__version__)
    except Exception:
        out["biopython"] = None
    return out


DEPS = deps()

# 启动时装载一次外部源（CBS 目录 / MycoBank），供审稿区做名字对照。
# 目录或异名表缺失时 SEARCHES 为空 → 名字对照显示「未配」，不报错。
_SRC_CACHE = {}


def sources_for(out_dir):
    """按输出目录找运行时生成的 _sources.json（网页版把它写在 out_dir 下）。
    找不到再退回仓库的 config/sources.json。结果按目录缓存，避免反复读 61 MB 的异名表。"""
    key = os.path.abspath(out_dir or ".")
    if key in _SRC_CACHE:
        return _SRC_CACHE[key]
    got = {}
    for cfgp in (os.path.join(key, "_sources.json"),
                 os.path.join(REPO, "config", "sources.json")):
        if os.path.exists(cfgp):
            try:
                from typeaudit import sources as S
                got = S.load(cfgp)
                break
            except Exception:
                got = {}
    _SRC_CACHE[key] = got
    return got


def logln(s):
    with LOCK:
        for line in str(s).rstrip("\n").split("\n"):
            LOG.append(line)


# ---------------------------------------------------------------- 流水线线程
def worker(genus, out_dir, limit, framework, species_list, cbs, email, cfg,
           mycobank="", biocol="", skip_done=False):
    STATE.update(running=True, stop=False, genus=genus, out_dir=out_dir,
                 rc=None, t0=time.time())
    LOG.clear()
    logln("=" * 60)
    logln("属名：{}    输出：{}    限制种数：{}".format(genus, out_dir, limit or "不限"))
    if cbs:
        logln("CBS 目录：{}".format(cbs))
    else:
        logln("⚠️ 未提供 CBS 目录 → 核验列将全部标「无目录」（不会有错误标注）")
    logln("=" * 60)
    # ⚠️ 邮箱必须在任何 S.load / 任何请求【之前】设好：
    #    sources._entrez() 读的是 os.environ["NCBI_EMAIL"]（因为 load() 不往适配器传 email）。
    if email:
        os.environ["NCBI_EMAIL"] = email
        logln("NCBI_EMAIL 已设：{}".format(email))
    else:
        logln("⚠️ 未填邮箱 → 联网请求会被 NCBI 拒（日志会报 SystemExit）")

    import contextlib

    class _Tee:
        """把流水线里的 print 也导进网页日志。
        ⚠️ 不重定向的话，pipeline.run() 的所有 print 只进服务进程的控制台，
           网页上什么都看不见 —— 实测就因此让人误以为「卡住」。"""
        def write(self, s):
            if s and s.strip():
                logln(s)
            return len(s)
        def flush(self):
            pass

    try:
        from typeaudit import pipeline as P

        # 运行时配置：以仓库 config/sources.json 为底，注入界面上填的路径。
        # ⚠️ 条件是 `not cfgpath`，**不能带 and cbs**：
        #    原来写成 `if not cfgpath and cbs` —— 于是【CBS 目录留空】时不写这个文件，
        #    cfgpath 保持 ""，P.run(config="") → load("") → **一个源都不装**，
        #    全表被标成「无目录」。这与界面上「留空 = 用自带」的承诺正好相反，
        #    而且日志只会淡淡说一句「可用源：（无外部源，仅联网检索）」，不会报错。
        cfgpath = cfg
        if not cfgpath:
            cfgpath = os.path.join(out_dir, "_sources.json")
            os.makedirs(out_dir, exist_ok=True)
            # ⚠️ 三个数据集【无论界面上填没填都写进配置】—— 填了用填的，
            #    留空写空串，由 sources.load() 回退到自动探测到的位置。
            #    原来只在【填了】才写，导致留空时该源根本不出现 ⇒
            #    机构登记表留空 → 非 CBS 机构全变「无目录」；
            #    MycoBank 留空 → 「目录记为别种」无法被解释，非错的号留在错值栏。
            # 【以仓库的 config/sources.json 为底，只覆盖界面上填的路径】
            # ⚠️ 不要在这里硬写源清单。原来硬写了 5 个，于是网页跑的流水线
            #    缺了 CGMCC 和 ICMP —— 三源交叉验证与「书写错误」纠错在网页端静默失效，
            #    而命令行（读 config/sources.json）却是正常的。**两套配置 = 两套行为**，
            #    而且新加适配器时没人会记得回来改这里。
            _cfg = json.loads(io.open(os.path.join(REPO, "config", "sources.json"),
                                      encoding="utf-8").read())
            srcs = {k: dict(v) for k, v in (_cfg.get("sources") or {}).items()}
            # 界面上填了的优先；留空写空串，由 sources.load() 回退到自动探测的位置
            for key, val in (("cbs", cbs), ("mycobank", mycobank),
                             ("biocollections", biocol)):
                if key in srcs:
                    if val and not os.path.isabs(val):
                        # 界面上填的相对路径按【仓库根】解析。原来交给 load() 的
                        # 两级相对解析，而那里的 base 是配置所在目录（= 输出目录），
                        # 于是 `data/xxx` 会被找成 `<out_dir>/data/xxx`。
                        val = os.path.normpath(os.path.join(REPO, val))
                    srcs[key]["path"] = val or ""
            # 【缓存根统一到仓库级 .cache】，不再按输出目录分家。
            # ⚠️ 原来写的是 os.path.join(out_dir, ".cache")：out_dir 是界面上的
            #    相对值（如 "out_audit"），load() 又按配置所在目录（也是 out_dir）
            #    拼一次 → 变成 out_audit/out_audit/.cache 双层目录（实测 277.41 MiB），
            #    而且每换一个输出目录就把属级 8.38 MiB 重下一遍。
            #    现在写【绝对路径】，与 CWD 和 out_dir 都无关。
            for k in ("genbank", "taxonomy", "cgmcc", "icmp"):
                if k in srcs:
                    srcs[k]["path"] = os.path.join(REPO, ".cache")
            logln("数据源（来自 config/sources.json）：{}".format("、".join(sorted(srcs))))
            # ⚠️ 网页上写着「机构登记表本仓库已内置，留空即用自带精简版」——
            #    但原来只在【填了】才写入配置，留空就【根本不写这一项】，
            #    于是所有非 CBS 机构的号被标成「无目录」而不是「未知机构代码」。
            #    本机实测同一批号两次跑出不同结果（一次 未知机构代码×7，一次 无目录×17），
            json.dump({"sources": srcs}, open(cfgpath, "w", encoding="utf-8"),
                      ensure_ascii=False, indent=1)

        # 让「停止」在种之间生效
        orig = P.collect

        def guard(*a, **kw):
            if STATE["stop"]:
                raise KeyboardInterrupt("用户停止")
            return orig(*a, **kw)

        P.collect = guard
        # 【跳过已完成的种】—— 让"服务重启后接着跑"变成一键操作。
        # 为什么需要：跑批在服务进程的线程里，重启服务/关机就中断；
        # 但每完成一个种已落盘到三列表，所以启动前把已在表里的种剔掉即可。
        if skip_done and species_list:
            try:
                txt = (io.open(species_list, encoding="utf-8-sig").read()
                       if os.path.exists(species_list) else species_list)
                _want = [ln.strip() for ln in txt.splitlines() if ln.strip()]
                _tbl = os.path.join(out_dir, "{}_三列表.tsv".format(genus.lower()))
                _done = set()
                if os.path.exists(_tbl):
                    for _l in io.open(_tbl, encoding="utf-8-sig").read().splitlines()[1:]:
                        if _l.strip():
                            _done.add(_l.split("\t")[0].strip())
                _rest = [s for s in _want if s not in _done]
                logln("⏭ 跳过已完成 {} 个种，本次跑 {} 个".format(len(_want) - len(_rest), len(_rest)))
                if not _rest:
                    logln("✅ 清单里的种全都已完成，无需重跑")
                    STATE["rc"] = 0
                    return
                _tf = os.path.join(out_dir, "_species_resume.txt")
                io.open(_tf, "w", encoding="utf-8", newline="\n").write("\n".join(_rest) + "\n")
                species_list = _tf
            except Exception as _e:
                logln("⚠️ 跳过已完成失败（照原清单跑）：{}: {}".format(type(_e).__name__, _e))

        with contextlib.redirect_stdout(_Tee()), contextlib.redirect_stderr(_Tee()):
            P.run(genus=genus, out_dir=out_dir, config=cfgpath,
                  framework=framework or None, species_list=species_list or None,
                  limit=(int(limit) if limit else None), email=email or None)
        STATE["rc"] = 0
        # 【守恒断言】跑完当场对账：有没有种从所有选项卡消失。
        # 放在这里而不是 pipeline 里，是因为"五个去向"的定义就在本文件
        # （read_candidates / read_trusted / read_undecidable）——
        # 换一份实现去算，就又变成"两套实现 = 两套 bug"。
        try:
            _rep = conservation_check(out_dir, genus)
            if _rep:
                logln("\n" + conservation_line(_rep))
        except Exception as _e:
            logln("\n[match] ⚠️ 对账未跑成（{}: {}）".format(type(_e).__name__, _e))
        logln("\n✅ 完成。")
        p = os.path.join(out_dir, "{}_三列表.tsv".format(genus.lower()))
        if os.path.exists(p):
            logln("📄 结果表：{}".format(p))
            logln("→ 切到「审稿」区块开始裁定。")
    except KeyboardInterrupt as e:
        STATE["rc"] = 130
        logln("\n⛔ 已停止：{}".format(e))
    except BaseException as e:
        # ⚠️ 必须抓 BaseException：sources._entrez() 在缺邮箱时 raise 的是
        #    **SystemExit**（不是 Exception），用 `except Exception` 抓不到 ——
        #    实测线程会静默死掉，界面上只剩一行「耗时 N 秒」，看起来像卡住。
        import traceback
        STATE["rc"] = 1
        logln("\n❌ 出错（{}）：{}".format(type(e).__name__, e))
        logln(traceback.format_exc()[-1500:])
    finally:
        STATE["running"] = False
        logln("（耗时 {:.1f}s）".format(time.time() - (STATE["t0"] or time.time())))
        # ★ 自动补全：主审计正常结束就自动做四方核验 + 反向审计 + 汇总。
        #   这三件按「种」计费，开销小；号连接按「号」计费，留给按钮手动触发。
        if STATE.get("rc") == 0:
            try:
                threading.Thread(target=enrich_worker, daemon=True,
                                 kwargs=dict(out_dir=out_dir, genus=genus, link=False)).start()
                logln("\n▶ 已自动启动补全：四方核验 + 反向审计 + 汇总")
                logln("   （号连接按号计费，需要时点「📋 全部发现」里的『一键补全（含号连接）』）")
            except Exception as e:
                logln("⚠️ 自动补全未启动：{}".format(e))


# ---------------------------------------------------------------- 跑树（并名裁决）
TREE_STATE = {"running": False, "stop": False, "rc": None, "pair": ""}
TLOG = deque(maxlen=3000)


# ---------------------------------------------------------------- CBS 目录下载
# 走 fetch_cbs_live.py 的子进程，逐行读它的 stdout —— 和手工在命令行跑完全同一条路径，
# 行为一致、崩溃隔离，而且已拉的页有磁盘缓存，中断后重跑会接着走。
CBS_STATE = {"running": False, "rc": None, "proc": None, "mode": "", "out_tsv": None}
CLOG = deque(maxlen=4000)
CLOCK = threading.Lock()


def cbs_log(s):
    with CLOCK:
        for line in str(s).rstrip("\n").split("\n"):
            CLOG.append(line)


def cbs_worker(genus, out_dir, email, page, sleep, no_syn, syn_min):
    """后台跑一次 CBS 目录下载。genus 为空 = 拉全量。"""
    CBS_STATE.update(running=True, rc=None, mode=(genus or "全量"), out_tsv=None)
    with CLOCK:
        CLOG.clear()
    script = os.path.join(REPO, "tools", "fetch_cbs_live.py")
    args = [sys.executable, "-u", script]
    if genus:
        args += ["--genus", genus]
    if no_syn:
        args += ["--no-synonyms"]
    if syn_min:
        args += ["--syn-min", str(syn_min)]
    if email:
        args += ["--email", email]
    if page:
        args += ["--page", str(page)]
    if sleep:
        args += ["--sleep", str(sleep)]
    if out_dir:
        args += ["--out-dir", out_dir]
    cbs_log("本次范围：{}".format(genus or "全量（约 63,000 条 / 30 次请求 / 10–40 分钟）"))
    cbs_log("$ " + " ".join(args))
    cbs_log("-" * 60)
    proc = None
    try:
        proc = subprocess.Popen(args, cwd=REPO, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, bufsize=1,
                                universal_newlines=True, encoding="utf-8", errors="replace")
        CBS_STATE["proc"] = proc
        for line in proc.stdout:
            cbs_log(line.rstrip("\n"))
            # 抓产物路径回填给前端 —— 用户点完下载就不用再手输路径了。
            # 只认主表（带 .tsv 且不含 diff），否则会抓到差异表的路径。
            if "已写" in line and ".tsv" in line and "diff" not in line:
                m = re.search(r"已写\s+(.+?\.tsv)", line)
                if m:
                    p = m.group(1).strip()
                    # 兜底转绝对路径：子进程的 cwd 是仓库根，相对路径按它解析。
                    # 前端直接把它填进输入框，用户随后跑流水线时会写进 _sources.json，
                    # 而 load() 解析相对路径是按配置目录找的 —— 不对齐就会静默用错文件。
                    if not os.path.isabs(p):
                        p = os.path.abspath(os.path.join(REPO, p))
                    CBS_STATE["out_tsv"] = p
        proc.wait()
        CBS_STATE["rc"] = proc.returncode
        cbs_log("-" * 60)
        cbs_log("结束，退出码 {}。已拉的页都在 cache/ 里，重跑会自动续。".format(proc.returncode))
    except BaseException as e:
        import traceback
        CBS_STATE["rc"] = 1
        cbs_log("出错（{}）：{}".format(type(e).__name__, e))
        cbs_log(traceback.format_exc()[-1200:])
    finally:
        CBS_STATE["running"] = False
        CBS_STATE["proc"] = None


def tree_log(s):
    with LOCK:
        for line in str(s).rstrip("\n").split("\n"):
            TLOG.append(line)


# ---------------------------------------------------------------- 反向审计
# 做什么：CBS 目录说「某一种的模式菌是 CBS X」→ 反过来问 GenBank 有没有声明这个号。
# 为什么要接进网页：这是「漏登」的唯一检测口，以前只能命令行单独跑，网页上看不到。
REV_STATE = {"running": False, "rc": None}
REV_LOG = []
REV_OUT = os.path.join(REPO, "out_reverse")


def rev_log(s):
    """写一行进度到内存缓冲（前端轮询 /api/revlog 取）。"""
    REV_LOG.append("[{}] {}".format(time.strftime("%H:%M:%S"), s))
    if len(REV_LOG) > 600:
        del REV_LOG[:150]


# ---------------------------------------------------------------- 号连接（连通分组）
# 做什么：把一批菌株号按「同一管实物」聚成组。两座桥：CBS.other_coll + ICMP.Duplicates。
LINK_STATE = {"running": False, "rc": None}
LINK_LOG = []
LINK_RESULT = {}


def link_log(s):
    LINK_LOG.append("[{}] {}".format(time.strftime("%H:%M:%S"), str(s)[:220]))
    if len(LINK_LOG) > 400:
        del LINK_LOG[:100]


def _numbers_from_table(out_dir, genus="colletotrichum"):
    """从三列表两列里抽出全部号（去标记、去重）。"""
    p = os.path.join(out_dir, "{}_三列表.tsv".format(genus))
    if not os.path.exists(p):
        return []
    seen, out = set(), []
    for line in io.open(p, encoding="utf-8-sig").read().splitlines()[1:]:
        f = line.split("\t")
        for ci in (1, 2):
            if len(f) <= ci or f[ci] in ("—", ""):
                continue
            for cell in f[ci].split(" ｜ "):
                m = re.match(r"^(.*?)\s*\[", cell)
                disp = (m.group(1) if m else cell).strip()
                if not disp or disp.startswith("(strain)"):
                    continue
                k = re.sub(r"\s+", "", disp).upper()
                if k in seen:
                    continue
                seen.add(k)
                out.append(disp)
    return out


# ---------------------------------------------------------------- 一键补全（四方核验 / 号连接 / 反向审计 / 汇总）
ENRICH_STATE = {"running": False, "rc": None}
ENRICH_LOG = []


def enrich_log(s):
    ENRICH_LOG.append("[{}] {}".format(time.strftime("%H:%M:%S"), str(s)[:240]))
    if len(ENRICH_LOG) > 800:
        del ENRICH_LOG[:200]


def enrich_worker(out_dir, genus, link=False):
    """跑 tools/enrich.py。

    快档（link=False）：四方核验 + 反向审计 + 汇总 —— 按【种】计费，开销小，
                      所以主审计一结束就自动跑（用户不需要再点任何按钮）。
    全档（link=True） ：再加号连接 —— 按【号】计费（1.5–4 秒/号），必须手动触发。
    """
    ENRICH_STATE.update(running=True, rc=None)
    ENRICH_LOG.clear()
    cmd = [sys.executable, os.path.join(HERE, "enrich.py"),
           "--out", out_dir, "--genus", genus]
    if link:
        cmd.append("--link")
    enrich_log("▶ 补全开始（{}；建树上限 3 对）".format("含号连接" if link else "不含号连接"))
    env = dict(os.environ)
    env["PYTHONPATH"] = os.path.join(REPO, "src") + os.pathsep + env.get("PYTHONPATH", "")
    env["PYTHONIOENCODING"] = "utf-8"
    if not env.get("NCBI_EMAIL"):
        env["NCBI_EMAIL"] = "researcher@example.org"
    # 进程内调用不读 env 字典，必须写进真实环境变量
    if not os.environ.get("NCBI_EMAIL"):
        os.environ["NCBI_EMAIL"] = env["NCBI_EMAIL"]
    # 建树：上限 3 对（每对 1–5 分钟）。写死在这里，让"自动"有明确边界。
    args = ["--out", out_dir, "--genus", genus, "--tree", "3"] + (["--link"] if link else [])
    try:
        rc = run_py_script("enrich.py", args, enrich_log, "补全")
        ENRICH_STATE.update(running=False, rc=rc)
        enrich_log("■ 补全结束（rc={}）".format(rc))
    except Exception as e:
        enrich_log("❌ {}".format(e))
        ENRICH_STATE.update(running=False, rc=1)


# ---------------------------------------------------------------- 全部发现（汇总）
FIND_STATE = {"running": False, "rc": None}
FIND_LOG = []


def find_log(s):
    FIND_LOG.append("[{}] {}".format(time.strftime("%H:%M:%S"), str(s)[:230]))
    if len(FIND_LOG) > 200:
        del FIND_LOG[:60]


def read_findings(out_dir, genus):
    """读 <out>/<genus>_发现汇总.tsv。返回 (路径, 行列表, 汇总)。"""
    p = os.path.join(out_dir, "{}_发现汇总.tsv".format((genus or "colletotrichum").lower()))
    if not os.path.exists(p):
        return p, [], {}
    lines = io.open(p, encoding="utf-8").read().splitlines()
    if not lines:
        return p, [], {}
    hdr = lines[0].split("\t")
    rows = [dict(zip(hdr, l.split("\t"))) for l in lines[1:] if l.strip()]
    cnt = {}
    for r in rows:
        k = (r.get("层级", ""), r.get("归类", ""))
        cnt["{}·{}".format(*k)] = cnt.get("{}·{}".format(*k), 0) + 1
    return p, rows, cnt


def collect_worker(out_dir, genus):
    """后台跑 tools/collect_findings.py。"""
    FIND_STATE.update(running=True, rc=None)
    FIND_LOG.clear()
    script = os.path.join(HERE, "collect_findings.py")
    find_log("▶ 汇总 {} / {}".format(out_dir, genus))
    env = dict(os.environ)
    env["PYTHONPATH"] = os.path.join(REPO, "src") + os.pathsep + env.get("PYTHONPATH", "")
    env["PYTHONIOENCODING"] = "utf-8"
    try:
        rc = run_py_script("collect_findings.py", ["--out", out_dir, "--genus", genus],
                           find_log, "汇总")
        FIND_STATE.update(running=False, rc=rc)
    except Exception as e:
        find_log("❌ {}".format(e))
        FIND_STATE.update(running=False, rc=1)


# ---------------------------------------------------------------- 自动三级裁决
AUTO_STATE = {"running": False, "rc": None}
AUTO_LOG = []
AUTO_ROWS = []


def auto_log(s):
    AUTO_LOG.append("[{}] {}".format(time.strftime("%H:%M:%S"), str(s)[:230]))
    if len(AUTO_LOG) > 600:
        del AUTO_LOG[:150]


def auto_worker(pairs, out_root, min_cols, per_name):
    """后台跑 tools/auto_verdict.py（三级降级：8 → 5 → ITS 单基因）。

    为什么走子进程而不是 import：它是带 argparse 的命令行工具，
    import 进来还得自己拼 Namespace，反而更脆。子进程还能把每行日志实时抛出来。
    """
    AUTO_STATE.update(running=True, rc=None)
    AUTO_LOG.clear(); AUTO_ROWS.clear()
    pairs_file = os.path.join(out_root, "_待裁决对.tsv")
    try:
        os.makedirs(out_root, exist_ok=True)
        with io.open(pairs_file, "w", encoding="utf-8", newline="") as f:
            for a, b in pairs:
                f.write("{}\t{}\n".format(a, b))
    except Exception as e:
        auto_log("❌ 写待裁决对失败: {}".format(e))
        AUTO_STATE.update(running=False, rc=1); return

    script = os.path.join(HERE, "auto_verdict.py")
    cmd = [sys.executable, script, "--pairs-file", pairs_file, "--out", out_root,
           "--min-cols", str(min_cols), "--per-name", str(per_name)]
    auto_log("▶ 开始：{} 对，命令 {}".format(len(pairs), " ".join(cmd[1:])))
    env = dict(os.environ)
    env["PYTHONPATH"] = os.path.join(REPO, "src") + os.pathsep + env.get("PYTHONPATH", "")
    env.setdefault("NCBI_EMAIL", "researcher@example.org")
    env["PYTHONIOENCODING"] = "utf-8"
    try:
        args = ["--pairs-file", pairs_file, "--out", out_root,
                "--min-cols", str(min_cols), "--per-name", str(per_name)]
        rc = run_py_script("auto_verdict.py", args, auto_log, "自动三级裁决")
        AUTO_STATE.update(running=False, rc=rc)
        auto_log("■ 子进程结束（rc={}）".format(rc))
    except Exception as e:
        auto_log("❌ 启动失败 {}: {}".format(type(e).__name__, e))
        AUTO_STATE.update(running=False, rc=1)
        return
    # 读汇总表
    try:
        sp = os.path.join(out_root, "自动裁决汇总.tsv")
        if os.path.exists(sp):
            lines = io.open(sp, encoding="utf-8").read().splitlines()
            hdr = lines[0].split("\t")
            for l in lines[1:]:
                f = l.split("\t")
                AUTO_ROWS.append(dict(zip(hdr, f)))
            auto_log("■ 汇总 {} 行".format(len(AUTO_ROWS)))
    except Exception as e:
        auto_log("⚠️ 读汇总失败: {}".format(e))


# ---------------------------------------------------------------- 四方核验
SPEC_STATE = {"running": False, "rc": None}
SPEC_LOG = []
SPEC_ROWS = []


def spec_log(s):
    SPEC_LOG.append("[{}] {}".format(time.strftime("%H:%M:%S"), str(s)[:230]))
    if len(SPEC_LOG) > 400:
        del SPEC_LOG[:100]


def spec_worker(names, out_dir):
    """后台跑四方核验：MycoBank / NCBI / CBS 目录 / CGMCC 官网，再对账。

    规则（用户定）：某家"空白"不计错；只有【双方都正面声明且互斥】才算矛盾；
    双方各指定不同模式号 → 判「待核」，不判错。
    """
    SPEC_STATE.update(running=True, rc=None)
    SPEC_LOG.clear(); SPEC_ROWS.clear()
    try:
        from typeaudit import sources as S
    except Exception as e:
        spec_log("❌ 导入失败: {}".format(e)); SPEC_STATE.update(running=False, rc=1); return
    srcs = sources_for(out_dir)
    spec_log("▶ 四方核验 {} 个种".format(len(names)))
    for i, nm in enumerate(names, 1):
        try:
            st = S.species_status(nm, srcs)
        except Exception as e:
            spec_log("  [{}] {} ❌ {}: {}".format(i, nm, type(e).__name__, str(e)[:60]))
            continue
        nc = st.get("number_check") or {}
        SPEC_ROWS.append({
            "种": nm,
            "MycoBank接受名": ((st.get("myco") or {}).get("current") or "")[:80],
            "NCBI_taxid": str((st.get("ncbi") or {}).get("taxid") or ""),
            "NCBI类型材料号": "; ".join((nc.get("ncbi_type_numbers") or [])[:5]),
            "CBS类型号": "; ".join((nc.get("cbs_type_numbers") or [])[:4]),
            "CGMCC该名号": "; ".join((nc.get("cgmcc_all") or [])[:4]),
            "CGMCC模式号": "; ".join((nc.get("cgmcc_type") or [])[:3]),
            "矛盾数": len(st.get("conflicts") or []),
            "待核数": len(st.get("pending") or []),
            "提示数": len(st.get("notes") or []),
            "矛盾明细": " ｜ ".join((st.get("conflicts") or [])[:2])[:260],
            "待核明细": " ｜ ".join((st.get("pending") or [])[:2])[:260],
        })
        spec_log("  [{}/{}] {:<34} 矛盾{} 待核{} 提示{}".format(
            i, len(names), nm.replace("Colletotrichum ", ""),
            len(st.get("conflicts") or []), len(st.get("pending") or []),
            len(st.get("notes") or [])))
    # 落盘
    try:
        p = os.path.join(out_dir, "四方核验.tsv")
        cols = ["种", "MycoBank接受名", "NCBI_taxid", "NCBI类型材料号", "CBS类型号",
                "CGMCC该名号", "CGMCC模式号", "矛盾数", "待核数", "提示数", "矛盾明细", "待核明细"]
        with io.open(p, "w", encoding="utf-8", newline="") as f:
            f.write("\t".join(cols) + "\n")
            for r in SPEC_ROWS:
                f.write("\t".join(str(r.get(c, "")).replace("\t", " ") for c in cols) + "\n")
        spec_log("■ 结果表 → {}".format(p))
    except Exception as e:
        spec_log("⚠️ 写表失败: {}".format(e))
    SPEC_STATE.update(running=False, rc=0)


def annotate_worker(out_dir, icmp_budget):
    """后台跑「连接反哺」：把连通组里干净 CBS 锚点的结论借给无目录/未核的号。

    安全边界（写在 sources.inherit_from_groups 里，不在这里放宽）：
    只借 CBS 锚点、锚点必须干净（白名单）、绝不继承「号不存在」。
    """
    LINK_STATE.update(running=True, rc=None)
    LINK_LOG.clear()
    try:
        from typeaudit import pipeline as P
    except Exception as e:
        link_log("❌ 导入 pipeline 失败: {}".format(e))
        LINK_STATE.update(running=False, rc=1)
        return
    link_log("▶ 连接反哺：{}".format(out_dir))
    try:
        r = P.annotate_links(out_dir, sources=sources_for(out_dir), icmp_budget=icmp_budget,
                             log=link_log)
    except Exception as e:
        link_log("❌ 失败 {}: {}".format(type(e).__name__, e))
        LINK_STATE.update(running=False, rc=1)
        return
    LINK_RESULT.update({"groups": [], "singletons": [], "edges": [],
                        "stats": {"annotate": {k: v for k, v in r.items() if k != "stats"}}})
    link_log("■ " + str({k: v for k, v in r.items() if k not in ("stats", "dst")}))
    LINK_STATE.update(running=False, rc=0 if r.get("ok") else 1)


def link_worker(numbers, out_dir, icmp_budget):
    """后台跑连通分组，结果写 TSV 并放进 LINK_RESULT。"""
    LINK_STATE.update(running=True, rc=None)
    LINK_LOG.clear()
    LINK_RESULT.clear()
    try:
        from typeaudit import sources as S
    except Exception as e:
        link_log("❌ 导入 typeaudit 失败: {}".format(e))
        LINK_STATE.update(running=False, rc=1)
        return
    srcs = sources_for(out_dir)
    link_log("▶ 开始：{} 个号（ICMP 预算 {}）".format(len(numbers), icmp_budget))
    link_log("   数据源：cbsmap={} icmp={}".format(
        "有" if srcs.get("cbsmap") else "无", "有" if srcs.get("icmp") else "无"))
    try:
        res = S.link_numbers(numbers, srcs, icmp_budget=icmp_budget,
                             progress=lambda n: link_log("   ICMP 已查 {} 个".format(n)))
    except Exception as e:
        link_log("❌ 失败 {}: {}".format(type(e).__name__, e))
        LINK_STATE.update(running=False, rc=1)
        return
    LINK_RESULT.update(res)
    st = res.get("stats") or {}
    link_log("■ 完成：{} 组，{} 个孤立号".format(len(res.get("groups") or []),
                                          len(res.get("singletons") or [])))
    link_log("   统计 {}".format(st))
    # 落盘
    try:
        p = os.path.join(out_dir, "号连接组.tsv")
        with io.open(p, "w", encoding="utf-8", newline="") as f:
            f.write("组号\t成员数\t成员\t机构代码\t连接来源\n")
            for i, g in enumerate(res.get("groups") or [], 1):
                codes = sorted({re.match(r"^([A-Za-z]+)", x).group(1).upper()
                                for x in g if re.match(r"^([A-Za-z]+)", x)})
                vias = sorted({(e.get("via") or "") for e in (res.get("edges") or [])
                               if e.get("a") in g and e.get("b") in g})
                f.write("{}\t{}\t{}\t{}\t{}\n".format(
                    i, len(g), "; ".join(g), "; ".join(codes), "; ".join(vias)))
        link_log("   结果表 → {}".format(p))
    except Exception as e:
        link_log("⚠️ 写表失败: {}".format(e))
    LINK_STATE.update(running=False, rc=0)


def reverse_worker(email):
    """后台跑 tools/reverse_audit.py，把 stdout 逐行转进日志缓冲。

    ⚠️ 用【子进程】而不是 import：那个脚本的 main() 里有嵌套的缓存函数，
       import 拿不到，且它是命令行工具、有自己的退出码。
    """
    REV_STATE.update(running=True, rc=None)
    REV_LOG.clear()
    script = os.path.join(HERE, "reverse_audit.py")
    if not os.path.exists(script):
        rev_log("❌ 找不到 {}".format(script))
        REV_STATE.update(running=False, rc=1)
        return
    rev_log("▶ 开始反向审计：CBS 类型号 → GenBank 是否声明")
    rev_log("   输出目录 {}".format(REV_OUT))
    env = dict(os.environ)
    env["PYTHONPATH"] = os.path.join(REPO, "src")
    try:
        rc = run_py_script("reverse_audit.py", ["--out", REV_OUT,
                                                "--email", email or "audit@example.org"],
                           rev_log, "反向审计")
    except Exception as e:
        rev_log("❌ 启动失败 {}: {}".format(type(e).__name__, e))
        rc = 1
    rev_log("■ 结束，退出码 {}".format(rc))
    REV_STATE.update(running=False, rc=rc)


def tree_worker(workdir, pair, acc_range, email, phylo_bin):
    """依次跑 type_synonym_verdict.py 的五个阶段（fetch→prep→align→tree→verdict）。

    为什么要分阶段跑而不是一条 all：**可断点续跑**——fetch/prep 都有缓存，
    断了再点一次不会重复下载 GenBank。
    """
    TREE_STATE.update(running=True, stop=False, rc=None, pair=" vs ".join(pair))
    TLOG.clear()
    script = os.path.join(HERE, "type_synonym_verdict.py")
    if not os.path.exists(script):
        tree_log("❌ 找不到 {}".format(script))
        TREE_STATE.update(running=False, rc=1)
        return
    tree_log("=" * 60)
    tree_log("并名裁决：{}  vs  {}".format(pair[0], pair[1]))
    tree_log("工作目录：{}".format(workdir))
    tree_log("号段：{}".format(acc_range))
    tree_log("=" * 60)
    try:
        stages = [["fetch", "--acc-range", acc_range],
                  ["prep", "--min-genes", "8"],
                  ["align", "--min-cols", "200"],
                  ["tree"],
                  ["verdict", "--pair", pair[0], pair[1]]]
        for st in stages:
            if TREE_STATE["stop"]:
                tree_log("⛔ 已停止（当前阶段跑完即停）")
                TREE_STATE["rc"] = 130
                break
            cmd = [sys.executable, script, st[0], "--workdir", workdir,
                   "--phylo-bin", phylo_bin, "--email", email] + st[1:]
            tree_log("\n$ " + " ".join(st))
            env = dict(os.environ)
            env["PYTHONPATH"] = SRC          # 让脚本能 import typeaudit
            env["NCBI_EMAIL"] = email
            p = subprocess.run(cmd, cwd=REPO, capture_output=True, text=True,
                               encoding="utf-8", errors="replace", env=env, timeout=None)
            for line in (p.stdout or "").splitlines():
                tree_log("  " + line)
            if p.returncode != 0:
                tree_log("  ⚠️ {} 退出码 {}".format(st[0], p.returncode))
                for line in (p.stderr or "").splitlines()[-12:]:
                    tree_log("  ! " + line)
                if st[0] in ("fetch", "prep"):
                    TREE_STATE["rc"] = 1
                    break
        else:
            TREE_STATE["rc"] = 0
            tree_log("\n✅ 五阶段跑完。裁决见上方 verdict 段（并名-1/2/3/4）。")
            vp = os.path.join(workdir, "verdict_{}_{}.tsv".format(
                re.sub(r"\W+", "", pair[0]), re.sub(r"\W+", "", pair[1])))
            if os.path.exists(vp):
                tree_log("📄 裁决表：" + vp)
    except Exception as e:
        import traceback
        TREE_STATE["rc"] = 1
        tree_log("❌ 出错（{}）：{}".format(type(e).__name__, e))
        tree_log(traceback.format_exc()[-1000:])
    finally:
        TREE_STATE["running"] = False


# ---------------------------------------------------------------- 读三列表 → 候选卡
BASIONYM_GENERA = {"vermicularia", "gloeosporium", "glomerella"}

# 【有性型属名】一菌两名：有性世代（teleomorph）用的属名，与无性世代属名指同一批物种。
# 2013 年「一菌一名」后，下列后者是接受名、前者并入。
# ⚠️ 为什么要表：原来只把 BASIONYM_GENERA 写死成炭疽菌的三个基名属，
#    换到 Trichoderma 后 Hypocrea 不在表里 → 40 条「Hypocrea xxx」被判成
#    「别的属 → 错值」。而它们只是同一物种的有性世代名，**不是错**。
TELEOMORPH_OF = {
    "hypocrea": "trichoderma",
    "glomerella": "colletotrichum",
    "nectria": "fusarium",
    "haematonectria": "fusarium",
    "gibberella": "fusarium",
    "calonectria": "cylindrocladium",
    "mycosphaerella": "pseudocercospora",
}


def _epithet_stem(e):
    """种加词词干：去掉随属名性别变化的词尾（-us/-um/-a/-is/-e）。

    实测例子：Hypocrea cerace**a** ↔ Trichoderma cerace**um** —— 同一个种，
    只因属名性别不同，词尾一个 a 一个 um。直接比字符串会判成"不同种"。
    """
    e = (e or "").lower()
    for suf in ("um", "us", "is", "ae", "a", "e"):
        if e.endswith(suf) and len(e) - len(suf) >= 4:
            return e[: -len(suf)]
    return e


def name_relation(number, claimed_species, sources):
    """给一个号，回答三件事：目录怎么记的、MycoBank 的接受名是什么、两者与声明名什么关系。

    审稿人看到「目录记为别种」时必须能立刻分清两种情形（否则等于把判断丢回给人瞎猜）：
      · 异名      —— 目录记的是**已并入**的名字（如 hymenocallidis → siamense）⇒ **非错**
      · 基名属    —— 目录记的是**基名属**（Vermicularia / Gloeosporium / Glomerella）⇒ **非错**
      · 真别种    —— 两个都是**现行接受名**且不同 ⇒ **可疑，须查原始出处**
    """
    out = {"catalog_name": None, "current_name": None, "verdict": None}
    cbs = sources.get("cbs") if sources else None
    myb = sources.get("mycobank") if sources else None
    m = re.match(r"^\s*([A-Za-z][A-Za-z\.\-]{0,12})\s*[:#]?\s*([\w\.\-]+)\s*$", number or "")
    if cbs and m and _cbs_ok(cbs):
        try:
            st = cbs.status_of(m.group(1), m.group(2))
            if st and st.get("name"):
                out["catalog_name"] = st["name"]
        except Exception:
            pass
    if out["catalog_name"] and myb:
        try:
            out["current_name"] = myb.current_name(out["catalog_name"])
        except Exception:
            pass
    cat = (out["catalog_name"] or "").strip()
    parts = cat.split()
    gen = (parts[0].lower() if parts else "")
    cat_epi = (parts[1].lower() if len(parts) > 1 else "")
    sp_epi = (claimed_species.split()[-1].lower() if claimed_species else "")
    decl_gen = (claimed_species.split()[0].lower() if claimed_species else "")
    cur = (out["current_name"] or "").strip()
    mcur = re.match(r"^[A-Z][a-z]+\s+([a-z][a-z\-]+)", cur)
    cur_epi = (mcur.group(1).lower() if mcur else "")
    # 种下单元（f. / var. / subsp.）：名字里的"末词"是命名人，不能拿来比种加词
    # ⚠️ 末尾不能再加 \b：`f.` 结尾是句点（非单词字符），`\.` 后面紧跟空格时 \b 不成立，
    #    正则整体失配 —— 实测 `Colletotrichum gloeosporioides f. alatae` 因此被误判成
    #    「并名冲突」，而不是「种下单元」。
    is_infra = bool(re.search(r"\b(f\.|f\.sp\.|var\.|subsp\.|forma)\s", cat + " ", re.I))

    # ⚠️ 必须先判「相同」。原来是先判「声称种出现在接受名里 → 异名」，
    #    于是目录名与声明种完全一致的号（84/120）全被标成「异名（非错）」——
    #    把"一致"报成了"异名"，方向正好反了。
    if not cat:
        pass
    elif gen in BASIONYM_GENERA:
        out["verdict"] = "基名属（非错）"
    elif cat_epi and sp_epi and cat_epi == sp_epi:
        out["verdict"] = "一致（机构记名与声明种相同）"
    elif is_infra:
        out["verdict"] = "CBS 记为种下单元（变型/变种，非错）"
    elif gen and decl_gen and gen != decl_gen:
        # 属名不同，三种可能必须分开（混在一起就会把"有性型"当错）：
        #   ① 有性型属名（Hypocrea ↔ Trichoderma、Glomerella ↔ Colletotrichum …）→ 非错
        #   ② 真·别的属 → 须查
        # 词尾按词干比（ceracea ↔ ceraceum 是同一种）。
        if TELEOMORPH_OF.get(gen) == decl_gen:
            if _epithet_stem(cat_epi) and _epithet_stem(cat_epi) == _epithet_stem(sp_epi):
                out["verdict"] = "有性型属名（非错）"
            else:
                out["verdict"] = "别的属（须查原始出处）"
        elif cur_epi and sp_epi and cur_epi == sp_epi:
            out["verdict"] = "异名（非错）"
        else:
            out["verdict"] = "别的属（须查原始出处）"
    elif cur_epi and cur_epi == cat_epi:
        out["verdict"] = "并名冲突（须查原始出处）"
    elif cur_epi and sp_epi and cur_epi == sp_epi:
        out["verdict"] = "异名（非错）"
    else:
        out["verdict"] = "并名冲突（须查原始出处）"
    # 若【声明名本身】在 MycoBank 是拼写变体（如 karstii → karsti），
    # 那"两个都接受"的判法是错的 —— 它只是拼写不同，属非错。
    try:
        myb0 = sources.get("mycobank") if sources else None
        if myb0 is not None and hasattr(myb0, "status_of"):
            # ⚠️ 要查【机构记的那个名字】的状态，不是声明种 ——
            #    karstii 是 karsti 的拼写变体，而声明种 boninense 是正常的 Legitimate。
            #    查错对象 → 这条永远判不出来（实测踩过）。
            for _n in (cat, claimed_species):
                if not _n:
                    continue
                st0 = (myb0.status_of(_n) or "").strip()
                if not st0:
                    continue
                out["name_status"] = st0
                out["name_status_of"] = _n
                # ⚠️ 只有在【接受名的种加词 == 声明种】时才叫拼写变体。
                #    否则会误判：karstii 是 karsti 的拼写变体没错，但声明种是 boninense ——
                #    这仍是"boninense vs karsti"两个不同种的冲突，不能降级为非错。
                #    （实测：不加这个条件，CBS 127536/127552 徽章写并名冲突、名字层写拼写变体，自相矛盾）
                if st0.lower().startswith("orthographic") and cur_epi and sp_epi and cur_epi == sp_epi:
                    out["verdict"] = "拼写变体（非错）"
                    break
    except Exception:
        pass
    return out


def _cbs_ok(obj):
    if obj is None:
        return False
    if hasattr(obj, "n_rows"):
        try:
            return obj.n_rows() > 0
        except Exception:
            return False
    return True


# ---------------------------------------------------------------- 标记释义（详情按钮用）
# 每个标记： (归类, 一句话定义, 判别依据) —— 归类决定它进哪一栏
FLAG_EXPLAIN = {
    "号不存在": ("错值候选", "该机构公开的机读数据库里查不到这个号。",
                 "查了该机构的公开数据库，返回空。可能是号写错、或该号未公开。"),
    "号有误（应为 X）": ("错值候选", "该号查不到，但把机构代码换成相近的、且确实存在的那家后能查到。",
                        "代码不在机构登记表里 + 近似代码在登记表里 + 用同一个号在新代码下确实查到 ⇒ 判定为写错。"),
    "CBS 记为别种": ("错值候选", "该机构把这个号记成【另一个种/属】，且 MycoBank 也解释不成声明种的异名。",
                     "机构记的名 → MycoBank 接受名 ≠ 声明种，且两者都是现行接受名 ⇒ 须查该菌株原始出处。"),
    "异名（非错）": ("非错", "该机构记的是声明种的【老名（异名）】，两者指同一种。",
                     "机构记的名 → MycoBank 接受名 == 声明种 ⇒ 非错。"),
    "基名属": ("非错", "该机构记的是【异属的老属名】（如 Vermicularia→Colletotrichum、Glomerella→Colletotrichum）。",
               "机构记的名 → MycoBank 接受名落在声明属内，且机构名的属 ≠ 声明属 ⇒ 非错。"),
    "书写错误（应为 X）": ("错值候选",
                          "机构代码前缀写错了 —— 更正后的代码 + 原号码确实查得到。",
                          "代码不在登记表 + 形近代码属于 CBS/CGMCC/ICMP + 用同一个号在更正后代码下查到 ⇒ 书写错误。"
                          "更正号已并入受信菌株表。"),
    "书写错误（形近 X，但该号在 X 查无）": ("错值候选",
                                          "机构代码不在登记表里，且形近 CBS/CGMCC/ICMP 之一；但那个号在更正后代码下也查不到。",
                                          "代码确非登记机构 ⇒ 判书写错误；给不出更正号 ⇒ 只报错不更正，须人工核原始出处。"),
    "未知机构代码": ("错值候选", "机构代码不在 NCBI Biocollections 登记表内。",
                     "代码查登记表返回空。**但登记表未必收录齐全，且本工具只对 CBS/CGMCC/ICMP 下确定结论**，故须人判。"),
    "疑似非机构编号": ("不确定", "这个代码既不在机构登记表里、也不像写错的机构代码 —— "
                                "多半是实验室内部编号或宿主标签被当成了机构代码。",
                       "判据：含小写字母 / 号码≤2位 / 登记表里无相近代码。**启发式，须人判**。"),
    "疑似应为 X（待人工核）": ("不确定", "机构代码不在登记表里，但登记表里有相近的机构代码。",
                              "只报「疑似」不报「应为」—— 因为用它下面的号去查并不能查到（否则会报「号有误」）。"),
    "标本号": ("不确定", "这是标本馆号（如 CBS H-21492），不该拿菌种保藏记录去核。",
               "号形如 H-xxxx 且机构是 CBS。标本要核的是标本馆档案，不是菌种保藏记录。"),
    "经映射": ("不确定",
               "这个号在别家机构没有公开可查的数据库，但 CBS 的记录里另一个号登记了它作为别名，"
               "于是借 CBS 那条记录来判。详见下方「映射到 / 该号在 CBS 记的名 / 模式地位」。",
               "只是换一家权威来看，**本身未经独立验证**；能证明「是同一管菌」，"
               "不能证明「它确实是模式材料」。"),
    "并名冲突": ("不确定",
                "该机构记的种与声明种不同，而两者都是现行接受名 —— 不是异名，也不一定是错。",
                "须查原始发表确认是同一菌还是两株。异名是\"该机构用了已并入的旧名\"，"
                "这一档是两个名字都还活着，性质不同。"),
    "机构旧代码": ("非错",
                   "这个代码是机构的历史名称（如 AS＝Academia Sinica，即 CGMCC 的前身），不是编造。",
                   "按现行代码能查到同一管菌即判非错；查不到才继续按疑似非机构编号处理。"),
    "模式地位冲突": ("错值",
                    "记录原文写的是正模（holotype/neotype/lectotype），而机构把这个号记成次要模式。",
                    "两边指的类型等级不同，至少一方写错了 —— 须查原始发表确认。"),
    "模式地位待核": ("不确定",
                    "机构把这个号记为次要模式，但 GenBank 记录里没有 /type_material 字段 —— 记录从未声明它是类型材料。",
                    "记录既没声明、机构却有类型标注，两者不构成冲突；只能说原来的『模式地位不符』不成立。"),
    "记录原文亦为 epitype": ("非错",
                             "GenBank 记录自己写的就是 epitype，机构也记 Epitype —— 两边一致。",
                             "这是措辞准确的情形，不是错。原判定只看『号出现在类型材料字段』就报冲突，属误报。"),
    "NCBI 记为模式材料": ("不确定",
                  "该机构对这个号的类型标注是【次要模式】——比如补充模式（Epitype，原模式无法鉴别时额外指定）、"
                  "等模、副模、权威菌株。这些都是**合法的类型材料**，但都不是正模（Holotype）。",
                  "记录侧只笼统写「模式材料」，没写是哪一级；机构侧明确是【次要模式】。"
                  "一个种可以同时有正模与补充模式，所以两者不一定冲突 —— 但等级不同，"
                  "读者会把它当正模引用。按本项目口径（会让读者误解即算错）记为错值；"
                  "要定案须核记录里 /type_material 的原文写的是哪一级。"),
    "无公开数据库": ("不确定", "机构真实存在（在登记表里），但它没有公开的机读数据库 ⇒ **无从核验**。",
               "**不等于号不存在，也不等于号有问题**。出路：自备该机构的数据 data/catalogue_<代码>.tsv。"),
    "未核": ("不确定", "号在该机构的记录里，但**没给它记模式地位标注** ⇒ 无法核验它是不是模式材料。",
             "**不等于非模式菌株**。实测 CBS 记录 62,740 行里只有 24.1% 有模式地位标注。"),
    "记为非模式菌株": ("错值候选", "记录把它当类型材料用，但该机构标注它是【非模式菌株】。",
                       "声明侧（GenBank/Taxonomy 的类型材料槽位）与该机构（如 CGMCC 官网 is_typestrain）矛盾。"),
}


def _norm_key(code, number):
    """惰性导入 norm_key。

    ⚠️ 不要在函数里裸写 norm_key(...)：pipeline_web 模块级【没有】导入它，
       而两处调用都包在 try/except 里 —— NameError 会被静默吞掉，
       表现为「详情里少了几行」，完全看不出原因（实测踩过）。
    """
    from typeaudit.sources import norm_key
    return norm_key(code, number)


RANK_CN = {
    "holotype": "正模（主模式）", "ex-holotype": "正模衍生培养物",
    "epitype": "补型（原模式材料不可用时指定的补充模式）",
    "ex-epitype": "补型的衍生培养物",
    "neotype": "新模式（原模式丢失后重新指定）",
    "isotype": "等模式（正模的复份）", "paratype": "副模式",
    "syntype": "全模式", "lectotype": "选模式（从全模中选定）",
    "type": "模式材料（该机构未细分类型）",
    "authentic strain": "权威菌株（作者鉴定过，非正式模式）",
    "ex-type": "模式衍生培养物",
}


def rank_cn(raw):
    """'Epitype of Colletotrichum acutatum' → '补型（…）'。取不到就原样返回。"""
    if not raw:
        return ""
    m = re.match(r"\s*([A-Za-z][A-Za-z\- ]*?)\s+of\b", raw, re.I)
    key = (m.group(1).strip().lower() if m else str(raw).strip().lower())
    return RANK_CN.get(key, key)


def mapped_detail(flag, sources):
    """把「经映射→CBS X」展开成通顺的叙述条目（徽章放最短事实，详情说清来龙去脉）。"""
    out = []
    m = re.search(r"经映射→\s*([A-Za-z]+)\s*([\w\.\-]+)", flag or "")
    if not m or not sources:
        return out
    code, num = m.group(1), m.group(2)
    rec = {}
    try:
        cbs = sources.get("cbs")
        if cbs is not None and hasattr(cbs, "status_of"):
            rec = cbs.status_of(code, num) or {}
    except Exception:
        rec = {}
    out.append(("该号在 CBS 记的名", (rec.get("name") or "").strip() or "—"))
    # ⚠️ parse_status() 只返回 rank / name 两个键，**没有 raw**。
    #    原文 status 要从目录行（bykey）里取，否则这里永远显示空。
    raw = ""
    try:
        cbs0 = sources.get("cbs")
        _r0 = (getattr(cbs0, "bykey", {}) or {}).get(_norm_key(code, num)) or {}
        raw = (_r0.get("status") or "").strip()
    except Exception:
        raw = ""
    rk = (rec.get("rank") or "").strip()
    if not raw and not rk:
        out.append(("CBS 记的模式地位", "为空（该记录没有模式地位标注）"))
    else:
        out.append(("CBS 记的模式地位",
                    "{}　原文：{}".format(rank_cn(raw or rk), raw or rk)))
    try:
        cbs = sources.get("cbs")
        r = (getattr(cbs, "bykey", {}) or {}).get(_norm_key(code, num)) or {}
        nc = (r.get("name_changes") or "").strip()
        mb = (r.get("mb_id") or "").strip()
        if nc:
            out.append(("该号曾用名（CBS 记录）", nc[:120]))
        if mb:
            out.append(("对应 MycoBank 编号", "MB#{}".format(mb)))
    except Exception:
        pass
    out.append(("这条能证明什么",
                "能证明「这个号在 CBS 的记录里对应 {} {}」；"
                "不能证明「它确实是模式材料」—— 映射只是换一家权威来看，本身未经独立验证。"
                .format(code, num)))
    return out


def cbs_extra(item, sources):
    """任何 CBS 号都试着带上「曾用名 / MycoBank 编号」——有才显示。

    这两列来自合并后的 CBS 目录（fetch_cbs_live.py --genus 抓的富字段）：
      · name_changes = 该号在 CBS 记录里曾经用过的名字 —— 并名/异名的直接线索
      · mb_id        = 对应的 MycoBank 编号 —— 可与 MycoBank 记录机械对齐
    覆盖不全（只有抓过富字段的那些行有），所以「没有」不等于「没变更过」。
    """
    out = []
    if not sources or not item:
        return out
    m = re.match(r"^\s*CBS\s*[:#]?\s*([\w\.\-]+)\s*$", item, re.I)
    if not m:
        return out
    try:
        cbs = sources.get("cbs")
        r = (getattr(cbs, "bykey", {}) or {}).get(_norm_key("CBS", m.group(1))) or {}
    except Exception:
        return out
    nc = (r.get("name_changes") or "").strip()
    mb = (r.get("mb_id") or "").strip()
    full = (r.get("taxon_full") or "").strip()
    if full:
        out.append(("CBS 全名（含命名人）", full[:110]))
    if nc:
        out.append(("该号曾用名（CBS 记录）", nc[:120]))
    if mb:
        out.append(("对应 MycoBank 编号", "MB#{}".format(mb)))
    return out


def read_candidates(out_dir, genus, sources=None):
    p = os.path.join(out_dir, "{}_三列表.tsv".format((genus or "").lower()))
    if not os.path.exists(p):
        return None, []
    cards = []
    with open(p, encoding="utf-8-sig") as f:
        head = f.readline().rstrip("\n").split("\t")
        for ln in f:
            fld = ln.rstrip("\n").split("\t")
            if len(fld) < 3:
                continue
            sp = fld[0]
            for col, colname in ((1, "菌株号"), (2, "模式菌")):
                for cell in fld[col].split("｜"):
                    cell = cell.strip()
                    if not cell or cell == "—":
                        continue
                    m = re.match(r"^(.*?)\s*\[(.*?)\]$", cell)
                    if not m:
                        continue
                    item, flags = m.group(1).strip(), m.group(2).strip()
                    if flags in ("—", ""):
                        continue
                    for flag in [x.strip() for x in flags.split(";") if x.strip()]:
                        # 徽章只放短名，完整内容留在 system_flag 里给详情用。
                        # 为什么：徽章容器窄，60 多字的标记会竖排成一列，读不了
                        # （实测「号不存在（声明:… ‖ 查无:…）」就是这种效果）。
                        _short = re.split(r"（|‖", flag, 1)[0].strip() or flag
                        # 从号里取机构代码（CBS / ICMP / CGMCC / MFLU …）。
                        # 为什么：界面上写「机构记的名」和写「目录记的名」一样模糊 ——
                        # 读者要知道到底是哪一家在说话（用户实测指出）。
                        _mi = re.match(r"^\s*([A-Za-z][A-Za-z\.\-']{0,14})\s*[:#]?\s*[\w\.\-]+\s*$", item)
                        _inst = _mi.group(1).upper() if _mi else ""
                        card = {"species": sp, "column": colname,
                                "item": item, "system_flag": flag,
                                "badge": _short, "inst": _inst,
                                "raw_flag": flag}
                        # 【统一判定】三列表里的标记与网页层现算的裁决可能打架 ——
                        # 实测同一号左边标「目录记为别种」(红)、右边裁决写「异名（非错）」(绿)。
                        # 原因是有两套判定器。这里统一到 name_relation 的裁决上，
                        # 并把原始标记留在 raw_flag 里备查（不丢信息）。
                        # 对所有卡都算 name_rel —— 它是本地字典查询，很便宜，
                        # 而详情页需要「目录记的名 / MycoBank 接受名」才能解释清楚。
                        if sources:
                            nrel = name_relation(item, sp, sources)
                            card["name_rel"] = nrel
                            # 统一判定：三列表标记与网页层裁决打架时，以裁决为准
                            # （实测同一号左边标「目录记为别种」红、右边裁决「异名（非错）」绿）
                            if flag in ("目录记为别种", "基名属") and nrel and nrel.get("verdict"):
                                card["system_flag"] = nrel["verdict"]
                        if _short == "无公开数据库" and _inst:
                            card["badge"] = "{} 无公开数据库".format(_inst)
                        cls, short, why = FLAG_EXPLAIN.get(
                            card["system_flag"], ("不确定", "（该标记尚无释义）", ""))
                        # 带参数的标记（如「经映射→CBS 112996」）没有精确条目，按前缀兜底；
                        # 不兜底详情页就写「该标记尚无释义」——等于没解释。
                        if short == "（该标记尚无释义）":
                            # ⚠️ 用【包含】不是【开头】：标记现在以机构名开头
                            #    （如「CBS 记为别种」），开头匹配会全部落空。
                            for _pre, _v in FLAG_EXPLAIN.items():
                                if _pre in card["system_flag"] and not _pre.endswith(("→", "（")):
                                    cls, short, why = _v
                                    break
                        card["cls"] = classify_flag(card["system_flag"])
                        card["detail"] = [
                            ("在哪里发现",
                             "{} 的「{}」列：{}".format(sp, colname, item)),
                            ("来源表", "{}（{}）".format(os.path.basename(p or ""),
                                                        _mtime(p) or "时间未知")),
                            ("系统标记", card["system_flag"] +
                             ("　（三列表原始标记：{}）".format(flag)
                              if card["system_flag"] != flag else "")),
                            ("归类", card["cls"]),
                            ("这个标记是什么意思", short),
                            ("为什么这么判", why),
                            ("{} 记的名".format(_inst or "机构"),
                             (card.get("name_rel") or {}).get("catalog_name") or "—"),
                            ("MycoBank 接受名", (card.get("name_rel") or {}).get("current_name") or "—"),
                            ("关系裁决", (card.get("name_rel") or {}).get("verdict") or "—"),
                            ("诚实边界", "本工具只核【号与机构公开记录是否一致】；"
                                        "不核序列是否属于该种、不核它是否真是模式菌株、不核论文声明。"),
                        ]
                        # 经映射的卡：把叙述补进详情（徽章只放最短事实，详情才说清来龙去脉）
                        for _lab, _val in mapped_detail(card["system_flag"], sources):
                            card["detail"].insert(-1, (_lab, _val))
                        # 任何 CBS 号都带上富字段（曾用名 / MyCoBank 编号），有才显示
                        for _lab, _val in cbs_extra(item, sources):
                            card["detail"].insert(-1, (_lab, _val))
                        cards.append(card)

    # ---------- 合并「经映射」：同一个 CBS 号只留一张卡 ----------
    # 一条事实被拆成多张卡（ATCC 56816、STE-U 5292 各有一张，CBS 112996 又一张）是重复；
    # 而且被映射的号往往还带着自己的标记卡，光删「经映射」那张只吃掉一半（实测踩过）。
    # 正确做法：**按号整组折叠** —— 凡是映射到 CBS X 的号，它的所有卡都并进 CBS X 那张，
    # 别名和它自己的标记都记在详情里。
    try:
        def _nkey(s):
            return str(s).upper().replace(" ", "")

        _items = {}
        for _c in cards:
            _items.setdefault((_c["species"], _nkey(_c["item"])), []).append(_c)
        # 每个号 → 它映射到的 CBS 号
        _map_tgt = {}
        for _c in cards:
            _m = re.match(r"经映射→\s*(CBS\s*[:#]?\s*[\w\.\-]+)",
                          str(_c.get("system_flag") or ""))
            if _m:
                _map_tgt[(_c["species"], _nkey(_c["item"]))] = " ".join(_m.group(1).split())
        _drop = set()
        for (_sp, _it), _tgt in _map_tgt.items():
            _mine = _items.get((_sp, _it)) or []
            if not _mine:
                continue
            _hosts = _items.get((_sp, _nkey(_tgt)))
            _alias = _mine[0]["item"]
            if not _hosts:
                # 该种下没有 CBS 号那张卡 —— 把「经映射」这张改挂到 CBS 号上
                _c0 = _mine[0]
                # ⚠️ 别把 item 直接改成目标号：那样卡片会显示成「CBS 123755」，
                #    看起来像是 123755 自己有问题，其实是另一个号映射过来的（实测被误读）。
                #    改成「原号 → 目标号」，来源一眼可见。
                _c0["item"] = "{} → {}".format(_alias, _tgt)
                _c0["badge"] = "经映射"
                _c0["inst"] = (_alias.split()[0].upper() if _alias.split() else "")
                _c0["detail"].insert(0, ("原始号",
                                         "{}（CBS 的记录里登记它为 {} 的别名）".format(_alias, _tgt)))
                for _c in _mine[1:]:
                    _drop.add(id(_c))
                continue
            _host = _hosts[0]
            _host.setdefault("aliases", []).append(_alias)
            _other = [str(x.get("system_flag")) for x in _mine
                      if not str(x.get("system_flag") or "").startswith("经映射")]
            if _other:
                _host.setdefault("alias_flags", []).append((_alias, _other))
            for _c in _mine:
                _drop.add(id(_c))
        if _drop:
            cards = [_c for _c in cards if id(_c) not in _drop]
        for _c in cards:
            if _c.get("aliases"):
                _c["aliases"] = sorted(set(_c["aliases"]))
                _c["detail"].insert(1, ("别名（经映射）",
                                        "{} —— 这些号在 CBS 记录里都指向本号".format(
                                            "、".join(_c["aliases"]))))
            for _alias, _fl in (_c.get("alias_flags") or []):
                _c["detail"].insert(-1, ("别名 {} 自己的标记".format(_alias),
                                         "；".join(_fl)[:160]))
    except Exception:
        pass

    # 【把「全部发现」里的错也变成卡片】
    # 为什么要这一步：原来只有三列表的【号层】标记会生成卡片，于是
    # 目录层的拼写错误、两家必有一方错的分歧，只躺在汇总表里，进不了审稿流程。
    # 用户实测反馈：「确定错」里只看到号不存在一种错误 —— 正是因为这些没变成卡。
    try:
        fp = os.path.join(out_dir, "{}_发现汇总.tsv".format((genus or "").lower()))
        if os.path.exists(fp):
            _L = io.open(fp, encoding="utf-8").read().splitlines()
            if _L:
                _h = _L[0].split("\t")
                for _l in _L[1:]:
                    d = dict(zip(_h, _l.split("\t")))
                    if (d.get("层级") or "号") == "号":
                        continue                      # 号层已由三列表生成，别重复
                    # 目录层【不生成卡片】：那几条是 CBS 目录【自身】的质量检查
                    # （拼写错误 95 处、占位符 2 处、撞车 1 处…），对象是"整个目录"，
                    # 不是某个物种/编号。放进逐条审稿流程没有意义 ——
                    # 「确认/驳回」对一个目录级统计无从下手（用户实测指出）。
                    # 它们留在「全部发现」表的目录层，报告里单列一节即可。
                    if (d.get("层级") or "") in ("目录", "源库"):
                        continue
                    _cls = d.get("归类") or ""
                    if _cls not in ("错值", "矛盾"):
                        continue                      # 只把"确定的错/矛盾"放进来
                    obj = d.get("对象") or ""
                    sp_ = obj.split(" / ")[0].strip() if " / " in obj else obj.strip()
                    item_ = obj.split(" / ")[-1].strip() if " / " in obj else "—"
                    flag_ = "{}（{}）".format(d.get("类别") or "", (d.get("层级") or "") + "层")
                    cards.append({
                        "species": sp_, "column": (d.get("层级") or "") + "层",
                        "item": item_, "system_flag": flag_,
                        "badge": d.get("类别") or flag_,
                        "raw_flag": d.get("发现") or "",
                        "cls": _cls, "from_findings": True,
                        "detail": [
                            ("在哪里发现", "{}（{}层）".format(obj or "—", d.get("层级") or "")),
                            ("来源表", os.path.basename(fp)),
                            ("系统标记", d.get("发现") or "—"),
                            ("归类", _cls),
                            ("这条的依据", d.get("证据") or "—"),
                            ("采集来源", d.get("来源") or "—"),
                            ("诚实边界", "本工具只核【号与机构公开记录是否一致】；"
                                        "不核序列是否属于该种、不核它是否真是模式菌株、不核论文声明。"),
                        ],
                    })
    except Exception:
        pass
    return p, cards


# ---------------------------------------------------------------- 留痕
def review_path(out_dir):
    os.makedirs(out_dir, exist_ok=True)
    return os.path.join(out_dir, "review_log.jsonl")


def append_review(out_dir, rec):
    rec["ts"] = time.strftime("%Y-%m-%d %H:%M:%S")
    with open(review_path(out_dir), "a", encoding="utf-8") as f:
        f.write(json.dumps(rec, ensure_ascii=False) + "\n")


# 「错值类」标记：只有这些才算候选错值。其余（无目录/未核/标本号/基名属）
# 是"无法核验"或"已排除的非错"，列进去只会淹没真问题。
ERROR_FLAGS = ("号不存在", "未知机构代码", "目录记为别种", "别的属",
               "属异名", "一株多名", "位数孤例")


# ---------------------------------------------------------------- 三档归类（诚实边界的机器表达）
#   错值   : 语料鉴定过、且判为出错  → 进【审稿区】，这是唯一可以下结论的地方
#   非错   : 语料鉴定过、且判为不出错 → 进【受信菌株表】
#   不确定 : **语料没判过的一律进这里** —— 不知道就说不确定，绝不猜
# 归类以 data/known_errors.tsv 的「归类」列为准；文件里没列的标记用下表兜底；
# 仍然没有的，默认「不确定」。
# 带动态值的标记用【前缀】匹配（如 `号有误（应为 CGMCC 3.14）`）
FLAG_CLASS_PREFIX = [
    ("书写错误", "错值"),          # 机构代码确非登记机构 ⇒ 真错值（含「应为 X」已证实那种）
    ("号有误（应为", "错值"),      # 旧写法，兼容已有产物
    # ⚠️ 必须加这条：`号不存在（形近 CBS 516.97 但目录记的是 X）` 这类【带后缀】的标记
    #    在兜底表里是精确匹配不到的。不加的话 classify_flag 会走到最后一档
    #    「不确定」——把明确的错值降级成不确定，判档就错了。
    ("号不存在", "错值"),
    # —— 两家名录说法相反 ——
    # 为什么归「错值」而不是「不确定」：两个名录对同一问题给出相反答案，
    # 这个矛盾本身客观、可复核 —— 至少一家在这件事上是错的（或过时）。
    # 但「哪一家错」需第三方证据（原始发表 / 模式菌号 / 序列），故标记写「方向待裁决」。
    ("权威分歧", "错值"),
    # NCBI 声明是模式号、但机构公开目录查无 —— **归不确定**（不是错）：
    # 实测 14 条这类号里 8 条属于 2022 年新发表的种，多半是「已保藏、目录未公开上线」。
    ("NCBI 声明为模式号，机构无此记录", "不确定"),
    # —— 跨库映射（无目录 → CBS 号）——
    # 归类为什么是「不确定」而不是「非错」：映射只证明 CBS 记录里登记过这个别名，
    # **不等于该机构自己的目录里有这个号**。按本项目纪律「没经语料鉴定的不进受信表」，
    # 所以它只能进「暂无法判定」——但格子里会显示映射到的 CBS 号及其判定，信息比「无目录」多。
    # 注意顺序：这条排在上面三条之后，所以「经映射 + 号不存在」这种组合仍会判成错值。
    ("经映射", "不确定"),
    # —— 并名冲突 ——
    # 为什么归「错值」（用户 2026-10-20 裁定）：目录记的种与声明种不同，而 MycoBank
    # 两个都还在接受名里 ⇒ **必有一方错**（记录标错，或机构目录滞后）。
    # 「错是确定的，方向未定」——所以进确定错，但标记里写明须查原始出处。
    ("并名冲突", "错值"),
    # 机构旧代码（如 AS = CGMCC 前身）：不是编造，是机构改名/合并留下的历史代码 → 非错
    ("机构旧代码", "非错"),
    # 有性型属名（Hypocrea→Trichoderma 等）：同一物种的另一个世代名，不是错
    ("有性型属名", "非错"),
    # —— 类型材料等级（2026-10-21 改进：核对记录原文）——
    # 只有「记录称正模、机构记次要模式」才算冲突（错值）；
    # 记录原文也是 epitype 的本就一致 → 非错；记录没有该字段 → 待核。
    ("模式地位冲突", "错值"),
    ("模式地位不精确", "不确定"),
    ("模式地位待核", "不确定"),
    ("记录原文亦为 epitype", "非错"),
]


def classify_flag(flag):
    """标记 → 归类。先精确匹配，再前缀匹配，最后保守判「不确定」。"""
    cls = _FLAG_CLASSES_CACHE.get("v")
    if cls is None:
        cls = flag_classes()
        _FLAG_CLASSES_CACHE["v"] = cls
    if flag in cls:
        v = cls[flag]
    else:
        v = "不确定"
        for pre, v2 in FLAG_CLASS_PREFIX:
            if str(flag).startswith(pre):
                v = v2
                break
    return _final_class(flag, v)


# 只有这两类留在「不确定」——它们的含义正是**我们核不了**：
#   未核   = 号在目录里，但目录没给它标模式地位
#   无目录  = 该机构根本没有公开目录
# 「模式地位待核」也属于"我们核不了"：机构给了类型标注，但 GenBank 记录里
# 根本没有 /type_material 字段，无法判断记录想说的是哪一级。
STAY_UNKNOWN = ("未核", "无公开数据库", "无目录", "模式地位待核")
# 明确判「非错」的不动 —— 那是"我们确定它没错"，塞进确定错会自相矛盾
KEEP_AS_IS = ("非错", "矛盾")


def _final_class(flag, v):
    """统一裁定（2026-10-20 用户定）：除「未核 / 无目录」外，一律进确定错。

    边界说明：这条规则会把「NCBI 声明为模式号但机构目录查无」「标本号」
    「经映射」这类原本归"不确定"的也变成错值 —— 用户的判据是"会让读者误解即算错"。
    但两类例外保持不变：
      · 非错（异名/基名属/种下单元）—— 我们确定它**不是**错
      · 矛盾（两家相反）—— 已进确定错，但保留"方向待裁决"的性质
    """
    if v in KEEP_AS_IS:
        return v
    # 标记文字里自己写明「非错」的（如「…（变型/变种，非错）」）也保持非错，
    # 否则会被下面"一律错值"的规则吃掉 —— 实测把两条非错判成了错值。
    if "非错" in str(flag or ""):
        return "非错"
    head = str(flag or "").split("（")[0].strip()
    if head in STAY_UNKNOWN or str(flag or "").startswith(STAY_UNKNOWN):
        return "不确定"
    return "错值"


_FLAG_CLASSES_CACHE = {}


FLAG_CLASS_FALLBACK = {
    "未知机构代码": "错值", "目录记为别种": "错值", "号不存在": "错值",
    "一株多名": "错值", "位数孤例": "错值", "属异名": "错值", "别的属": "错值",
    "异名（非错）": "非错", "基名属": "非错",
    "标本号": "不确定", "疑似非机构编号": "不确定", "无目录": "不确定", "未核": "不确定",
}


def flag_classes():
    """标记 → 归类。以语料为准，语料没写的用兜底表，再没有就是「不确定」。"""
    cls = dict(FLAG_CLASS_FALLBACK)
    p = os.path.join(REPO, "data", "known_errors.tsv")
    if os.path.exists(p):
        try:
            with io.open(p, encoding="utf-8") as f:
                head = f.readline().rstrip("\n").split("\t")
                if "期望标记" in head and "归类" in head:
                    i_f, i_c = head.index("期望标记"), head.index("归类")
                    for line in f:
                        q = line.rstrip("\n").split("\t")
                        if len(q) > max(i_f, i_c) and q[i_f].strip() and q[i_c].strip():
                            cls[q[i_f].strip()] = q[i_c].strip()
        except Exception:
            pass
    return cls


def read_trusted(out_dir, genus):
    """
    受信菌株表：三列表里【没有任何问题标记】的号。三列：种名 / 模式菌株 / 菌株号。

    ⚠️ 「受信」的含义**仅限于「号与外部目录/登记表一致」** ——
        不等于这株菌真是模式菌株，也不等于序列属于该种。界外的事本工具答不了。
    """
    p = os.path.join(out_dir, "{}_三列表.tsv".format((genus or "").lower()))
    rows = []
    if not os.path.exists(p):
        return None, rows
    def nk(x):
        return re.sub(r"[^A-Za-z0-9\.]", "", (x or "").upper())

    cls = flag_classes()
    with io.open(p, encoding="utf-8-sig") as f:
        f.readline()
        for line in f:
            fld = line.rstrip("\n").split("\t")
            if len(fld) < 3:
                continue
            cells = {"strain": [], "type": []}     # 每列：(展示号, 该格的标记列表)
            allf = {}                              # 归一化号 -> 全部标记（跨两列汇总）
            for col, k in ((1, "strain"), (2, "type")):
                for cell in fld[col].split("｜"):
                    cell = cell.strip()
                    if not cell or cell == "—":
                        continue
                    m = re.match(r"^(.*?)\s*\[(.*?)\]\s*$", cell)
                    if not m:
                        continue
                    raw = m.group(1).strip()
                    # ⚠️ (strain) 开头的是 GenBank /strain 的值（实验室编号），
                    #    没有机构前缀、【根本没被检查过】—— 照收就是把「没检查」说成
                    #    「检查通过」。一律排除。
                    if raw.startswith("(strain)"):
                        continue
                    fl = [x.strip() for x in m.group(2).split(";")
                          if x.strip() and x.strip() != "—"]
                    cells[k].append((raw, fl))
                    allf.setdefault(nk(raw), []).extend(fl)
            # ⚠️ 受信规则（本机实测修正过一次）：
            #    原来按「有没有标记」判定，于是把「非错」类标记（异名（非错）/
            #    基名属）也一并排除 —— 可它们的归类就是【非错】，按选项卡设计
            #    本该进受信表。实测 ICMP 19118 就是这样被漏掉的。
            #    正确规则：**该号在任一列的所有标记都属「非错」类**（或没有标记）才算受信；
            #    只要出现一个「错值」或「不确定」类标记，就出局。
            # 先扫一遍「书写错误（应为 X Y）」—— 更正已证实（同号在更正后代码下确实查到）
            # ⇒ 把【更正后的号】收进受信表。依据：更正形式与目录一致。
            fixed = {}
            for k in ("strain", "type"):
                for raw, fl_ in cells[k]:
                    for f_ in fl_:
                        m_ = re.match(r"^书写错误（应为\s*(\S+)\s+(\S+?)\s*）$", f_)
                        if m_:
                            fixed[nk(raw)] = (m_.group(1), m_.group(2))
            out, why = {}, {}
            for k in ("strain", "type"):
                seen, keep = set(), []
                for raw, fl_ in cells[k]:
                    key = nk(raw)
                    if key in seen:
                        continue
                    seen.add(key)
                    flags_here = set(allf.get(key, []))
                    bad = [x for x in flags_here if classify_flag(x) in ("错值", "不确定")]
                    if bad:
                        # 原号带错值/不确定标记 ⇒ 本身不进受信；
                        # 但若它是「书写错误且更正已证实」，把更正后的号放进来。
                        if key in fixed:
                            c_, n_ = fixed[key]
                            show = "{} {}".format(c_, n_)
                            keep.append(show)
                            why[show] = "由 {} 更正而来（该号在 {} 下确实查得到）".format(raw, c_)
                        continue
                    keep.append(raw)
                    good = [x for x in flags_here if classify_flag(x) == "非错"]
                    if good:
                        why[raw] = "；".join(sorted(good))
                out[k] = "; ".join(keep)
            if out["strain"] or out["type"]:
                # ⚠️ 种名一律用【结果表里的全名】。原来这里缩写成 "C. xxx"，
                #    而 read_candidates 用的是全名 —— 两个选项卡口径不一致，
                #    以后做「点种名联动筛选」就会匹配不上。缩写交给界面层。
                rows.append({"species": fld[0],
                             "type": out["type"], "strain": out["strain"],
                             "why_type": "; ".join(
                                 "{}:{}".format(k, v) for k, v in why.items()
                                 if k in out["type"]) if out["type"] else "",
                             "why_strain": "; ".join(
                                 "{}:{}".format(k, v) for k, v in why.items()
                                 if k in out["strain"]) if out["strain"] else ""})
    return p, rows


def read_undecidable(out_dir, genus):
    """
    暂无法判定：三列表里【一个可核验的号都没有】的种。

    为什么必须有这个选项卡（本机实测，3 个种的那一轮）：
      · `(strain)` 开头的是 GenBank `/strain` 值（论文内部编号，如 ZHKUCC 210077），
        没有机构前缀、【从来没被检查过】⇒ 收进受信表就是把「没检查」说成「检查通过」，
        所以设计上一律排除；
      · 某一列若为空，写表时落成占位符「—」，同样不算一条。
      ⇒ 「只有 (strain) + 空列」的种在【确定错 / 不确定 / 受信菌株表】三处都不出现，
        **凭空消失**。而它和「检查过、干净」的种外观完全一样 ——
        这是最坏的一类静默：你分不清「无从检查」与「检查通过」。
      ⇒ 单独立一个选项卡显式列出来，让五个去向的种数之和等于结果表的种数。
    """
    p = os.path.join(out_dir, "{}_三列表.tsv".format((genus or "").lower()))
    rows = []
    if not os.path.exists(p):
        return None, rows
    with io.open(p, encoding="utf-8-sig") as f:
        f.readline()
        for line in f:
            fld = line.rstrip("\n").split("\t")
            if len(fld) < 3:
                continue
            checked, strain_only, empty_cols = [], [], []
            for col, cn in ((1, "菌株号"), (2, "模式菌")):
                got = False
                for cell in fld[col].split("｜"):
                    cell = cell.strip()
                    if not cell or cell == "—":
                        continue
                    got = True
                    m = re.match(r"^(.*?)\s*\[(.*?)\]\s*$", cell)
                    raw = m.group(1).strip() if m else cell
                    if raw.startswith("(strain)"):
                        strain_only.append(raw.replace("(strain)", "").strip())
                    else:
                        checked.append(raw)
                if not got:
                    empty_cols.append(cn)
            # 只要有一个【可核验的号】，就交给别的选项卡（确定错 / 不确定 / 受信）
            if checked:
                continue
            why = []
            if strain_only:
                uniq = sorted(set(strain_only))
                why.append("只有 GenBank /strain 的论文内部编号 {} 个（{}）；"
                           "/strain 没有机构前缀，本工具从不对它做核验".format(
                               len(uniq), "、".join(uniq[:6])))
            if empty_cols:
                why.append("{}列为空".format("、".join(empty_cols)))
            if not why:
                why.append("两列都没有任何编号")
            rows.append({"species": fld[0],
                         "strain_only": "; ".join(sorted(set(strain_only))),
                         "n_strain_only": len(set(strain_only)),
                         "reason": "；".join(why)})
    return p, rows


def read_progress(out_dir, genus):
    """读进度边车文件 `<属>_进度.json`。没有文件就返回 None（旧产物）。

    为什么界面必须读它：一份跑了一半的三列表，格式与内容和「全属都查过了」完全一样。
    不显式提示，它就会被当成完整结论 —— 这是本项目的诚实边界原则在「部分结果」上的延伸。
    """
    p = os.path.join(out_dir, "{}_进度.json".format((genus or "").lower()))
    if not os.path.exists(p):
        return None
    try:
        d = json.load(io.open(p, encoding="utf-8"))
        d["_file"] = p
        return d
    except Exception:
        return None


def conservation_check(out_dir, genus):
    """【守恒断言】结果表里的种，有没有从所有选项卡里凭空消失的。

    规矩与由来：
      · 本机实测踩过 —— `Colletotrichum acaciae` 只有 GenBank `/strain` 内部编号、
        模式菌列为空，于是在【确定错/不确定/受信菌株表】三处都不出现，
        **静默消失**；而它和「检查过、干净」的种外观一模一样。
      · ⚠️ 五个去向**不是对「种」的划分**：一个种可以同时出现在多个选项卡里
        （如它有一个号干净、另一个号有问题：干净号进受信表，问题号进不确定）。
        所以断言只能是【覆盖】——每个种至少出现在一个去向里 ——
        **不能**写成「各选项卡种数之和 = 种数」，那样会误报。

    返回 dict；`uncovered` 非空即为 FAIL。
    """
    p = os.path.join(out_dir, "{}_三列表.tsv".format((genus or "").lower()))
    if not os.path.exists(p):
        return None
    sp_all = []
    with io.open(p, encoding="utf-8-sig") as f:
        f.readline()
        for line in f:
            c = line.rstrip("\n").split("\t")
            if c and c[0].strip():
                sp_all.append(c[0].strip())
    _, cands = read_candidates(out_dir, genus)
    _, trust = read_trusted(out_dir, genus)
    _, na = read_undecidable(out_dir, genus)
    sp_err, sp_unk = set(), set()
    for c in cands:
        sp = c.get("species")
        if str(c.get("cls") or "").strip() == "错值":
            sp_err.add(sp)
        else:
            sp_unk.add(sp)
    sp_trust = set(x.get("species") for x in trust)
    sp_na = set(x.get("species") for x in na)
    covered = sp_err | sp_unk | sp_trust | sp_na
    uncovered = [s for s in sp_all if s not in covered]
    return {
        "file": p, "genus": genus,
        "progress": read_progress(out_dir, genus),
        "n_species": len(sp_all),
        "n_err": len(sp_err), "n_unk": len(sp_unk),
        "n_trust": len(sp_trust), "n_na": len(sp_na),
        "uncovered": uncovered,
        "ok": not uncovered,
    }


def conservation_line(rep):
    """把守恒结果排成一行，给流水线日志用。"""
    if not rep:
        return "[match] 没有结果表，无法对账"
    head = ("[match] 种数 {n_species} ｜ 确定错 {n_err} ｜ 不确定 {n_unk} ｜ "
            "受信 {n_trust} ｜ 暂无法判定 {n_na} ｜ 未覆盖 {nunc}").format(
        nunc=len(rep["uncovered"]), **rep)
    _pg = rep.get("progress") or {}
    _tail = ""
    if _pg and _pg.get("finished") is False:
        _tail = "  ｜  ⚠️ 部分结果（{}/{} 个种，未完）".format(
            _pg.get("done"), _pg.get("total"))
    if rep["ok"]:
        return head + "  →  ✅ 无种遗漏" + _tail
    return (head + "  →  ❌ 有 {} 个种从所有选项卡消失了：\n       {}".format(
        len(rep["uncovered"]), "\n       ".join(rep["uncovered"][:12])) + _tail)


def read_errors(out_dir, genus, sources=None):
    """错值表：从三列表筛出真错值类候选，附机构记名/接受名对照。"""
    p = os.path.join(out_dir, "{}_三列表.tsv".format((genus or "").lower()))
    if not os.path.exists(p):
        return None, []
    rows = []
    with io.open(p, encoding="utf-8-sig") as f:
        head = f.readline().rstrip("\n").split("\t")
        for ln in f:
            fld = ln.rstrip("\n").split("\t")
            if len(fld) < 3:
                continue
            for col, colname in ((1, "菌株号"), (2, "模式菌")):
                for cell in fld[col].split("｜"):
                    cell = cell.strip()
                    m = re.match(r"^(.*?)\s*\[(.*?)\]$", cell)
                    if not m:
                        continue
                    item, flags = m.group(1).strip(), m.group(2).strip()
                    hits = [x.strip() for x in flags.split(";") if x.strip()]
                    if not any(any(h.startswith(e) for e in ERROR_FLAGS) for h in hits):
                        continue
                    rec = {"species": fld[0], "column": colname, "item": item,
                           "flags": "; ".join(hits)}
                    if sources:
                        rec.update(name_relation(item, fld[0], sources) or {})
                    rows.append(rec)
    return p, rows


def read_reviews(out_dir):
    p = review_path(out_dir)
    out = []
    if os.path.exists(p):
        for ln in open(p, encoding="utf-8"):
            ln = ln.strip()
            if ln:
                try:
                    out.append(json.loads(ln))
                except Exception:
                    pass
    return out


# ---------------------------------------------------------------- 页面
PAGE = r"""<!DOCTYPE html><html lang="zh"><head><meta charset="utf-8">
<title>类型材料审计 · 流水线 + 审稿</title>
<style>
 body{font:14px/1.5 -apple-system,"Segoe UI","Microsoft YaHei",sans-serif;margin:0;background:#f5f6f8;color:#222}
 header{background:#1f2937;color:#fff;padding:10px 16px;font-weight:600}
 .wrap{display:flex;gap:12px;padding:12px;align-items:flex-start;flex-wrap:wrap}
 .card{background:#fff;border:1px solid #e3e6ea;border-radius:8px;padding:12px;box-shadow:0 1px 2px rgba(0,0,0,.04)}
 .left{flex:1 1 380px;min-width:340px} .right{flex:2 1 620px;min-width:420px}
 label{display:block;font-size:12px;color:#555;margin:6px 0 2px}
 input[type=text],input[type=number]{width:100%;box-sizing:border-box;padding:5px 7px;border:1px solid #cbd2d9;border-radius:5px}
 .row{display:flex;gap:8px}.row>*{flex:1}
 button{padding:6px 14px;border:0;border-radius:5px;background:#2563eb;color:#fff;cursor:pointer;font-size:13px}
 button.gray{background:#6b7280} button.red{background:#dc2626} button:disabled{opacity:.45;cursor:not-allowed}
 pre#log{background:#0f172a;color:#d7dde8;padding:10px;border-radius:6px;height:320px;overflow:auto;font:12px/1.45 Consolas,monospace;white-space:pre-wrap}
 .alias{font-size:11px;color:#6b7280;font-weight:400;margin-left:4px}
.flag{display:inline-block;padding:1px 7px;border-radius:10px;font-size:11px;background:#fee2e2;color:#991b1b;margin-left:6px;max-width:none;white-space:normal;overflow:visible;text-overflow:clip;line-height:1.5;vertical-align:middle}  /* 不截断：宁可换行也要把话说完（用户要求）。 */
 .flag.ok{background:#dcfce7;color:#166534}.flag.warn{background:#fef3c7;color:#92400e}.flag.gray{background:#e5e7eb;color:#374151}
 .item{display:flex;align-items:center;gap:8px;padding:7px 9px;border:1px solid #eceff3;border-radius:6px;margin-bottom:6px}
 .item .sp{width:190px;font-size:12px;color:#4b5563;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}
 .item .nm{flex:1;font-family:Consolas,monospace;font-size:13px}
 .item .cl{width:54px;font-size:11px;color:#6b7280}
 .item .st{width:auto;min-width:56px;max-width:46%;flex:0 0 auto;font-size:11px;color:#6b7280}  /* 原来固定 64px：徽章一长就被压成竖排、再截断（实测「NCBI 记为模式材料，CBS 记为次要模式」被切成两行小方块）。改成按内容宽度，最多占 46%。 */
 .btns button{padding:3px 9px;font-size:12px;margin-left:4px}
 .b-ok{background:#16a34a}.b-no{background:#dc2626}.b-wait{background:#d97706}
 .done{opacity:.55}
 .tabb{background:#e5e7eb;color:#111;margin-right:6px}
 .tabb.on{background:#dc2626;color:#fff}
 table.et{width:100%;border-collapse:collapse;font-size:12px}
 table.et th{background:#f3f4f6;text-align:left;padding:5px 7px;border-bottom:1px solid #e5e7eb}
 table.et td{padding:5px 7px;border-bottom:1px solid #f1f3f5;vertical-align:top}
 table.et code{font-family:Consolas,monospace}
 .nmrel{font-size:11px;color:#4b5563;background:#f9fafb;border-left:2px solid #d1d5db;
        padding:2px 6px;margin:0 8px;flex:1 1 240px}
 .group{border:1px solid #e5e7eb;border-radius:8px;margin-bottom:10px;overflow:hidden}
 .ghead{background:#f3f4f6;padding:6px 10px;font-size:13px;display:flex;align-items:center;gap:8px}
 .ghead b{font-weight:600}
 .ghead .cnt{font-size:11px;color:#6b7280;background:#e5e7eb;border-radius:9px;padding:1px 7px}
 .ghead .todo{font-size:11px;color:#92400e;background:#fef3c7;border-radius:9px;padding:1px 7px}
 .gbody{padding:8px 10px}
 .item{margin-bottom:5px}
 .item:last-child{margin-bottom:0}
 .tiny{font-size:11px;color:#6b7280}
 .filters{margin:4px 0 8px;display:flex;flex-wrap:wrap;gap:4px;align-items:center}
  .filters button{padding:2px 8px;font-size:12px;white-space:nowrap}.filters button{background:#e5e7eb;color:#111;margin-right:5px}
 .filters button.on{background:#2563eb;color:#fff}
</style></head><body>
<header>类型材料审计 · 流水线 + 审稿系统　<span class="tiny" style="color:#9ca3af">机器只出候选，人裁定；裁定一律留痕</span></header>
  <div class="card" id="sumcard" style="margin:10px 16px 0;padding:7px 12px">
    <div style="display:flex;align-items:center;gap:10px;flex-wrap:wrap">
      <b style="font-size:13px;white-space:nowrap">三属审计结论</b>
      <span id="sumline" class="tiny" style="color:#374151">加载中…</span>
      <a href="javascript:void(0)" onclick="toggleSum()" id="sumtoggle"
         class="tiny" style="color:#2563eb;margin-left:auto;white-space:nowrap">展开详情 ▾</a>
    </div>
    <div id="sumbody" class="tiny" style="display:none;margin-top:8px;line-height:1.85;
         border-top:1px dashed #d1d5db;padding-top:8px"></div>
  </div>
<div id="depbar" style="display:none;background:#fef3c7;border-bottom:1px solid #fcd34d;
     padding:8px 16px;font-size:13px;color:#92400e"></div>
<div class="wrap">
 <div class="left card">
  <div><label>属名（必填）</label><input type="text" id="genus" value="Colletotrichum"></div>
  <div class="row">
    <div><label>限制种数（先填 3 试）</label><input type="number" id="limit" value="3"></div>
    <div><label>NCBI 邮箱</label><input type="text" id="email" placeholder="you@example.org"></div>
  </div>
  <details open style="margin:4px 0">
    <summary class="tiny" style="cursor:pointer;color:#6b7280">第三方数据源（只在自己重跑时需要填；看结论不必填）</summary>
  <label>CBS 目录（不填 → 不会有错误标注）</label>
  <div style="display:flex;gap:8px;align-items:center">
    <input type="text" id="cbs" placeholder="...\cbs_catalogue.tsv" style="flex:1">
    <button style="white-space:nowrap" onclick="dlCbs()">⇩ 下载目录</button>
  </div>
  <div class="tiny" id="cbsstat" style="margin-top:3px"></div>
  <label>MycoBank 异名表（填了才能判「CBS 记为别种」是异名还是真别种）</label>
  <input type="text" id="mycobank" placeholder="...\mycobank_full.tsv">
  <div class="tiny" style="margin-top:3px">
    没有？去 MycoBank 下载 <b>MBlist</b> 压缩包：
    <a href="https://www.mycobank.org/Basic%20names%20search" target="_blank"
       style="color:#2563eb">https://www.mycobank.org/Basic names search</a>
    ，然后解析：<code>python tools/import_mycobank.py &lt;下载的.zip&gt; -o data/mycobank_full.tsv</code>
  </div>
  <div class="tiny" style="margin-top:6px">
    机构登记表：<b>仓库已内置</b>（<code>data/biocollections_codes.csv</code>，9,672 行），自动使用，无需填。
  </div>
  <div class="tiny" id="refhint" style="margin-top:2px"></div>
  </details>
  <label>输出目录</label><input type="text" id="out" value="out_audit">
  <div class="tiny" style="margin-top:5px;line-height:1.9;padding:6px 8px;
       background:#ecfdf5;border:1px solid #a7f3d0;border-radius:6px">
    <b>不想跑流水线？</b> 直接载入仓库里已有的结果（不联网、不占 CPU）：
    <select id="loadres" onchange="loadResults(this.value)" style="font-size:12px;padding:2px 4px">
      <option value="">— 选择属 —</option>
      <option value="Colletotrichum|results/Colletotrichum">Colletotrichum（408 种 · 3,633 个号）</option>
      <option value="Trichoderma|results/Trichoderma">Trichoderma（445 种 · 4,200 个号）</option>
      <option value="Diaporthe|results/Diaporthe">Diaporthe（507 种 · 3,270 个号）</option>
    </select>
  </div>
  <label>框架论文 / 物种清单（可选）</label>
  <input type="text" id="framework" placeholder="框架论文 txt（可选）">
  <label>物种清单（可选）—— <b>每行一个种名</b>；填了就不联网取种名，直接跑这些种</label>
  <textarea id="species" rows="4" style="width:100%;box-sizing:border-box;padding:5px 7px;
    border:1px solid #cbd2d9;border-radius:5px;font-family:Consolas,monospace;font-size:12px"
    placeholder="Colletotrichum siamense&#10;Colletotrichum falcatum&#10;Colletotrichum arxii"></textarea>
      <label style="display:inline-flex;align-items:center;gap:5px;margin-top:4px"><input type="checkbox" id="skipdone" checked style="width:auto;margin:0"> <span>跳过已完成的种（接着上次跑）</span></label>
  <div style="margin-top:10px">
    <button id="run" onclick="run()">▶ 开始</button>
    <button id="stop" class="gray" onclick="stop()" disabled>■ 停止</button>
    <button id="openapi" class="gray" onclick="loadCands()">↻ 载入结果表</button>
  </div>
  <div class="tiny" style="margin-top:6px">停止 = 当前种跑完后停（联网请求无法硬中断）</div>
  <pre id="log"></pre>
 </div>
 <div class="right card">
  <div style="display:flex;align-items:center;gap:10px">
    <div id="warnbox"></div>
    <b>审稿区</b><span class="tiny" id="stat"></span>
    <span style="flex:1"></span>
    <button class="gray" onclick="loadCands()">↻ 刷新候选</button>
    <button class="gray" onclick="exportReviews()">⤓ 导出裁定记录</button>
  </div>
  <div style="margin:2px 0 8px">
    <button id="tabErr" class="tabb on" onclick="setTab('err')">确定错（审稿区）</button>
    <button id="tabUnk" class="tabb" onclick="setTab('unk')">不确定</button>
    <button id="tabTrust" class="tabb" onclick="setTab('trust')">受信菌株表</button>
    <button id="tabNA" class="tabb" onclick="setTab('na')">暂无法判定</button>
    <button id="tabCov" class="tabb" onclick="setTab('cov')">📋 覆盖声明</button>
    <button id="tabT" class="tabb" onclick="setTab('tree')">🌳 并名裁决</button>
    <span class="tiny" id="errcnt"></span>
  </div>
  <div id="trustview" style="display:none"></div>
  <div id="naview" style="display:none"></div>
  <div id="covview" style="display:none"></div>
  <div id="errview" style="display:none"></div>
  <div id="treeview" style="display:none">
    <div style="border:1px solid #d1d5db;border-radius:8px;padding:10px;margin-bottom:12px;background:#fbfdff">
      <div class="tiny" style="margin-bottom:6px">
        <b>自动三级（推荐）</b>：只填两个种加词 → 自动取两株的模式保藏号 → 自动判定能凑几个位点
        （<b>8 → 5 → ITS 单基因</b>逐级降级）→ 自动建树并给出裁决。
        <b>结论强度随层级下降</b>：8 位点「较强」、5 位点「中等」、ITS 单基因<b>只能「仅提示」</b>。
      </div>
      <div><label>争议对（每行一对：种加词 空格 种加词）</label>
        <textarea id="autopairs" rows="2" style="width:100%;font:12px/1.5 Consolas,monospace"
          placeholder="celtidis karsti&#10;analogum camelliae"></textarea></div>
      <div style="margin-top:6px">
        <button onclick="runAuto()">▶ 自动三级</button>
        <span class="tiny" id="autosum" style="margin-left:8px"></span>
      </div>
      <div id="autotable" style="margin-top:6px"></div>
      <pre id="autolog" style="background:#0f172a;color:#d7dde8;padding:8px;border-radius:6px;
        height:170px;overflow:auto;font:12px/1.45 Consolas,monospace;white-space:pre-wrap;margin-top:6px"></pre>
    </div>
    <div class="tiny" style="border-top:1px dashed #d1d5db;padding-top:8px;margin-bottom:6px">
      <b>手动模式（高级）</b>：需要自己给一个「同批提交的连号段」。
    </div>
    <div class="tiny" style="margin-bottom:6px">
      建树比较两个保藏号是否同种：<b>并名-2 异株·遗传支持</b> / <b>并名-3 遗传不支持</b> / <b>并名-4 未裁</b>。
      需要本机 WSL 的 phylo 环境（iqtree2 / mafft / trimal）。一次约 2–4 分钟。
    </div>
    <div class="row">
      <div><label>号 A</label><input type="text" id="ta" placeholder="ICMP:18578"></div>
      <div><label>号 B</label><input type="text" id="tb" placeholder="ICMP:18642"></div>
    </div>
    <div class="row">
      <div><label>工作目录</label><input type="text" id="twork" value="tree_work"></div>
      <div><label>号段（同批提交）</label><input type="text" id="trange" value="JX009400:JX010500"></div>
    </div>
    <div><label>phylo bin 目录（WSL 路径）</label>
      <input type="text" id="tphylo" value="/home/yelu/miniforge3/envs/phylo/bin"></div>
    <div style="margin-top:8px">
      <button onclick="runTree()">▶ 跑树</button>
      <button class="gray" onclick="fetch('/api/treestop',{method:'POST'})">■ 停止</button>
    </div>
    <pre id="tlog" style="background:#0f172a;color:#d7dde8;padding:10px;border-radius:6px;
      height:340px;overflow:auto;font:12px/1.45 Consolas,monospace;white-space:pre-wrap;margin-top:8px"></pre>
  </div>
  <div id="findview" style="display:none">
    <details style="margin-bottom:6px">
      <summary class="tiny" style="cursor:pointer;color:#2563eb">📖 这一页怎么用 · 能查出哪些类别（点开看）</summary>
    <div class="tiny" style="margin-top:4px;background:#f8fafc;border-left:3px solid #93c5fd;padding:5px 8px">
      <b>一次跑完就够了。</b>点「▶ 开始」会自动串完：号层 → 名字层（并名异名）→
      记录层（GenBank 声明）→ 汇总。<b>不需要你分层操作</b>；下面的「层级」列只是告诉你
      这条证据来自哪一步。<br>
      <b>本表默认不含号层。</b>号层（三列表逐号标记）在下方「逐号卡片」里已经一条一条列出、并且可以逐条裁定；本表只列卡片里没有的部分：名字层、序列层、记录层、源库层。想看全部，勾选上面的「含号层」。<br><b>发现不都在这一页。</b>「确定错」「不确定」两页是<b>逐号卡片</b>（可裁定）；
      本表是<b>四层汇总</b>（号 / 名字 / 记录 / 目录），含卡片里没有的
      <b>并名异名分歧</b>、<b>GenBank 无声明</b>、<b>CBS 目录自身问题</b>。
      找不齐时先看本表。
    </div>
    <div id="findkinds" class="tiny" style="margin:6px 0;line-height:1.9"></div>
    <div class="tiny" style="margin-bottom:6px">
      把<b>四个来源</b>的发现合成一张表：号层（三列表标记）｜名字层（并名/异名分歧）｜
      记录层（GenBank 声明缺失）｜目录层（CBS 目录自身问题）。
      每条都带<b>层级 / 归类 / 来源 / 证据</b>。<b>不确定</b>是诚实边界，不是没查。<br><br>
      <b>三档归类</b>　<b>确定错</b> = 语料鉴定过、判为出错（只有这一栏能下结论）；
      <b>不确定</b> = 语料没判过的一律进这里 —— 这就是我们的诚实边界；
      <b>受信菌株表</b> = 没有任何问题标记的号。界外的事（序列是否属于该种、
      是否真模式菌株、论文声明）本工具答不了。
    </div>
    </details>
    <div style="margin-bottom:6px">
      <button onclick="runFind()">↻ 重新汇总</button>
      <button onclick="runEnrich(false)">↻ 只重跑核验（不重跑审计）</button>
      <button class="gray" onclick="runEnrich(true)">↻ 只重跑核验 + 号连接</button>
      <span class="tiny" id="findsum" style="margin-left:8px"></span>
    </div>
    <details style="margin-bottom:6px">
      <summary class="tiny" style="cursor:pointer;color:#2563eb">📋 表格：卡片管不了的发现
        （异名/独立种分歧 · 一株多名 · GenBank 无声明 · 目录自身问题）—— 点开查看</summary>
    <div style="margin-top:6px" class="tiny">
      筛选：<button class="gray" onclick="setFindFilter('all')">全部</button>
      <button class="gray" onclick="setFindFilter('错值')">只看错值</button>
      <button class="gray" onclick="setFindFilter('矛盾')">只看矛盾</button>
      <button class="gray" onclick="setFindFilter('不确定')">只看不确定</button>
      <label style="margin-left:12px;cursor:pointer;color:#6b7280">
        <input type="checkbox" id="showhao" style="vertical-align:-2px"
               onchange="SHOW_HAO=this.checked;renderFindTable()"> 含号层（与下方卡片重复）
      </label>
    </div>
    <div id="findtable" style="max-height:520px;overflow:auto"></div>
    </details>
    <pre id="findlog" style="display:none;background:#0f172a;color:#d7dde8;padding:8px;
      border-radius:6px;height:120px;overflow:auto;font:12px/1.45 Consolas,monospace;
      white-space:pre-wrap;margin-top:8px"></pre>
  </div>


  <div class="filters">
    <button data-f="all" class="on" onclick="setFilter('all',this)">全部</button>
    <button data-f="todo" onclick="setFilter('todo',this)">待裁定</button>
    <button data-f="确认" onclick="setFilter('确认',this)">已确认</button>
    <button data-f="驳回" onclick="setFilter('驳回',this)">已驳回</button>
    <button data-f="存疑" onclick="setFilter('存疑',this)">存疑</button>
    <input type="text" id="kw" placeholder="搜索种名/号..." oninput="render()" style="width:150px;padding:2px 6px;font-size:12px">
  </div>
  <div id="cards"></div>
  <div id="dlg" style="display:none;position:fixed;left:0;top:0;right:0;bottom:0;
       background:rgba(0,0,0,.35);z-index:99" onclick="if(event.target===this)closeDetail()">
    <div style="background:#fff;max-width:760px;margin:6vh auto;border-radius:10px;
         padding:16px 18px;max-height:82vh;overflow:auto;box-shadow:0 8px 30px rgba(0,0,0,.2)">
      <div style="display:flex;align-items:center;gap:10px;margin-bottom:8px">
        <b>详情</b><span style="flex:1"></span>
        <button class="gray" onclick="closeDetail()">关闭</button>
      </div>
      <div id="dlgbody"></div>
    </div>
  </div>
  <div id="empty" class="tiny">还没有候选。先跑一次流水线，或点「刷新候选」（会去读输出目录里的 &lt;属&gt;_三列表.tsv）。</div>
 </div>
</div>
<script>
let CANDS=[], REVIEWS={}, FILTER='all', CUR='', SHOW_UNKNOWN=false;
let FLAGCLASS={}, TAB='err', TRUSTED=null, TSRC={};
const $=id=>document.getElementById(id);
function q(){return new URLSearchParams({genus:$('genus').value,out:$('out').value});}
async function run(){
  // ⚠️ 不要再引用 biocol —— 那个输入框已从界面移除（机构登记表现在由后端自动用内置的）。
  //    $('biocol') 返回 null，取 .value 会抛 TypeError，点「开始」直接中断。
  const body={genus:$('genus').value,limit:$('limit').value,email:$('email').value,
   cbs:$('cbs').value,mycobank:$('mycobank').value,biocol:'',
   out:$('out').value,
   framework:$('framework').value,species:$('species').value,skip_done:$('skipdone').checked};
  const r=await fetch('/api/run',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify(body)});
  const j=await r.json();
  if(!j.ok){alert(j.msg||'启动失败');return;}
  $('run').disabled=true;$('stop').disabled=false;poll();
}
async function stop(){await fetch('/api/stop',{method:'POST'});}
async function poll(){
  const r=await fetch('/api/log');const j=await r.json();
  $('log').textContent=j.log.join('\n');$('log').scrollTop=$('log').scrollHeight;
  $('stat').textContent=j.status;
  if(j.deps){renderDeps(j.deps);}
  if(j.running){setTimeout(poll,700);}
  else{$('run').disabled=false;$('stop').disabled=true;loadCands();}
}
async function loadCands(){
  CUR=$('genus').value;
  const r=await fetch('/api/candidates?'+q());const j=await r.json();
  CANDS=j.cards||[]; FLAGCLASS=j.classes||{};
  window.__PROG=j.progress||null;
  window.__SRC={file:j.file,mtime:j.mtime};
  const rv=await fetch('/api/reviews?'+q());const jv=await rv.json();
  REVIEWS={}; (jv.reviews||[]).forEach(x=>{REVIEWS[key(x)]=x;});
  render();
}
function key(c){return c.species+'|'+c.column+'|'+c.item+'|'+c.system_flag;}
function setFilter(f,btn){FILTER=f;document.querySelectorAll('.filters button[data-f]').forEach(b=>b.classList.remove('on'));btn.classList.add('on');render();}
async function fillRefs(){
  try{
    const r=await fetch('/api/refs'); const j=await r.json();
    // biocol 的输入框已移除（机构登记表由后端自动用内置的），不再往里填
    const pairs=[['cbs',j.cbs],['mycobank',j.mycobank]];
    let n=0;
    pairs.forEach(([id,v])=>{ const el=$(id); if(el&&v&&!el.value){el.value=v;n++;} });
    const hint=$('refhint');
    if(hint) hint.textContent = n? `（已从 registry 自动填入 ${n} 个路径）`
                                 : '（未找到 registry，请手动填或先跑 tools/find_refdata.py）';
  }catch(e){}
}
function renderDeps(d){
  const bar=$('depbar');
  if(d.biopython){ bar.style.display='none'; return; }
  bar.style.display='block';
  bar.innerHTML='<b>⚠️ 缺 biopython —— 点「开始」会失败。</b>　当前解释器：<code>'+d.python+'</code><br>'+
    '修法一（免安装）：用带 biopython 的解释器启动本服务，例如<br>'+
    '<code>"C:\\Users\\Administrator\\.dsh\\dsh-bio-genie\\python-env\\Scripts\\python.exe" tools\\pipeline_web.py</code><br>'+
    '修法二（装一次）：<code>"'+d.python+'" -m pip install biopython</code>';
}
function flagClass(f){
  if(f==='号不存在'||f.indexOf('记为别种')>=0)return '';
  if(f==='基名属'||f==='标本号')return 'ok';
  if(f==='无目录'||f==='未核')return 'warn';
  return 'gray';
}
function toggleUnknown(el){SHOW_UNKNOWN=!!el.checked;render();}
// 归类由后端解析好（带动态值的标记如「号有误（应为 …）」前端做字符串匹配会错）。
// 语料没判过的一律「不确定」—— 这个保守默认在后端。
function clsOf(c){ return (c && c.cls) || '不确定'; }
function curList(){ return CANDS.filter(c => clsOf(c) === ((TAB==='err')?'错值':'不确定')); }
// 判定词的缩写：徽章与短句用短形式，完整说法留在详情里。
function shortVerdict(v){
  v=String(v||'');
  if(v.indexOf('一致')===0) return '一致';
  if(v.indexOf('并名冲突')===0||v.indexOf('并名冲突')===0) return '并名冲突，须查发表';
  if(v.indexOf('种下单元')>=0) return '种下单元（非错）';
  if(v.indexOf('拼写变体')>=0) return '拼写变体（非错）';
  if(v.indexOf('基名属')>=0) return '基名属（非错）';
  if(v.indexOf('异名')===0) return '异名（非错）';
  if(v.indexOf('别的属')>=0) return '别的属，须查发表';
  return v;
}
// 从 "Colletotrichum arboricola M. Zapata ..., 2018" 里取年份
function _yr(s){const m=String(s||'').match(/(1[89]\d{2}|20\d{2})\s*$/);return m?m[1]:''; }
function shortBadge(b){
  b=String(b||'');
  if(b==='并名冲突'||b==='并名冲突') return '并名冲突';
  return b;
}
function render(){
  const kw=$('kw').value.trim().toLowerCase();
  const base=curList();
  let list=base.filter(c=>{
    const rv=REVIEWS[key(c)];
    if(FILTER==='todo'&&rv)return false;
    if(['确认','驳回','存疑'].includes(FILTER)&&(!rv||rv.verdict!==FILTER))return false;
    if(kw&&!((c.species+c.item+c.system_flag).toLowerCase().includes(kw)))return false;
    return true;});
  $('empty').style.display=base.length?'none':'block';
  const done=Object.keys(REVIEWS).length;
  const nErr=CANDS.filter(c=>clsOf(c)==='错值').length;
  const nUnk=CANDS.filter(c=>clsOf(c)==='不确定').length;
  // ⚠️ 部分结果必须自报家门：跑了一半的表和"全属都查过了"的表长得一样。
  //    不标注就会被当成完整结论（本项目的诚实边界原则）。
  let warn='';
  if(window.__PROG && window.__PROG.finished===false){
    warn=`<div class="tiny" style="margin-bottom:6px;background:#fef2f2;border-left:3px solid #dc2626;padding:7px 9px">
      ⚠️ <b>这是部分结果</b>：只跑了 <b>${window.__PROG.done}</b> / ${window.__PROG.total} 个种
      （${window.__PROG.updated||''}）。<b>不要把它当成全属结论</b> —— 没跑的种在这里根本没有候选。
      接着跑不会重复下载已完成的种。
    </div>`;
  }
  $('warnbox').innerHTML=warn;
  $('stat').innerHTML='<b>确定错 '+nErr+' 条</b> · 不确定 '+nUnk+' 条 · 已裁定 '+done
    +' 条 · 当前显示 '+list.length+' 条'+(CUR?' · '+CUR:'')+
    (window.__SRC&&window.__SRC.mtime?`<br><span class="tiny">表：${window.__SRC.file}（生成于 ${window.__SRC.mtime}）</span>`:'');
  // ---- 按【种】分组：每个种一个区块，组内逐号一行（保留逐号裁定）----
  const groups={};
  list.forEach(c=>{ (groups[c.species]=groups[c.species]||[]).push(c); });
  $('cards').innerHTML=Object.keys(groups).sort().map(sp=>{
    const items=groups[sp];
    const todo=items.filter(c=>!REVIEWS[key(c)]).length;
    const body=items.map(c=>{
      const k=encodeURIComponent(key(c));const rv=REVIEWS[key(c)]||{};
      return `<div class="item ${rv.verdict?'done':''}">
        <div class="cl">${c.column}</div>
        <div class="nm">${c.item}${(c.aliases&&c.aliases.length)?`<span class="alias">（别名 ${c.aliases.join('、')} 映射到此号）</span>`:''}</div>
        <div class="st"><span class="flag ${flagClass(c.system_flag)}">${shortBadge(c.badge||c.system_flag)}</span></div>
        ${c.name_rel && c.name_rel.catalog_name ?
          // 这一行是【名字层】判定：只回答"机构记的名字与声明种是不是同一个名字"，
          // 与【模式地位】判定（左侧徽章）是两件事，不能混为一谈。
          // 「一致」原来被染成黄色（warn），看着像有问题 —— 其实是最正常的状态，改绿色。
          // 写成一句通顺的话，不再用悬空小标签 —— 徽章讲「类型地位」，
          // 这一行只讲「名字对不对得上」，两者维度不同，用不同的句式区分开。
          // 三方并列：声明种 ｜ 机构记的名 ｜ MycoBank 接受名 → 裁决（缩写）。
          // 之前只写机构记的名，读者看不出"声明的是什么种"，读不懂裁决从哪来。
          `<div class="nmrel"><span style="color:#6b7280">名字</span>：NCBI 声明 <b>${c.species.split(' ').slice(-1)[0]}</b>` +
          `　｜　${c.inst||'机构'} 记 <b>${c.name_rel.catalog_name.split(' ').slice(0,2).join(' ')}</b>` +
          (c.name_rel.current_name?`　｜　MycoBank <b>${c.name_rel.current_name.split(' ').slice(0,2).join(' ')}</b>`:'') +
          (c.name_rel.verdict?`　→　<span style="color:${(c.name_rel.verdict.includes('非错')||c.name_rel.verdict.indexOf('一致')===0)?'#16a34a':'#b45309'}" title="${c.name_rel.verdict}（完整说法见详情）">${shortVerdict(c.name_rel.verdict)}</span>`:'') +
          `</div>`:''}
        <div style="flex:1"></div>
        <div class="btns">
          <button class="b-ok" onclick="judge('${k}','确认')">确认</button>
          <button class="b-no" onclick="judge('${k}','驳回')">驳回</button>
          <button class="b-wait" onclick="judge('${k}','存疑')">存疑</button>
          <button class="gray" onclick="showDetail('${k}')">详情</button>
        </div></div>`;}).join('');
    return `<div class="group">
      <div class="ghead"><b title="${sp}">${sp.replace('Colletotrichum ','C. ')}</b>
        <span class="cnt">${items.length} 条</span>
        ${todo?`<span class="todo">待裁定 ${todo}</span>`:`<span class="cnt">已裁定完</span>`}</div>
      <div class="gbody">${body}</div></div>`;}).join('');
}
async function judge(k,verdict){
  const c=CANDS.find(x=>encodeURIComponent(key(x))===k);if(!c)return;
  const note=prompt('裁定依据/备注（可留空）：', (REVIEWS[key(c)]||{}).note||'');
  if(note===null)return;
  const body={genus:CUR,out:$('out').value,species:c.species,item:c.item,column:c.column,
    system_flag:c.system_flag,verdict:verdict,note:note,reviewer:localStorage.getItem('reviewer')||''};
  const r=await fetch('/api/review',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify(body)});
  const j=await r.json();
  if(j.ok){REVIEWS[key(c)]=body;render();}else{alert('保存失败：'+(j.msg||''));}
}
function toggleSum(){
  const b=$('sumbody'), a=$('sumtoggle');
  const open=b.style.display==='none';
  b.style.display=open?'block':'none';
  a.textContent=open?'收起 ▴':'展开详情 ▾';
}
async function loadSummary(){
  const b=$('sumbody'), ln=$('sumline');
  try{
    const r=await fetch('/api/summary'); const d=await r.json();
    if(d.error){ ln.innerHTML='<span style="color:#b91c1c">未找到 results/summary.json</span>'; return; }
    ln.textContent=d.headline||'';
    let h='';
    for(const k of (d.blocks||[])){ h+='<div style="margin-bottom:6px"><b>'+k.title+'</b><br>'+k.lines.join('<br>')+'</div>'; }
    b.innerHTML=h||'（summary.json 里没有 blocks）';
  }catch(e){ ln.innerHTML='<span style="color:#b91c1c">读取失败：'+e+'</span>'; }
}
async function loadResults(v){
  if(!v) return;
  const p=v.split('|');
  $('genus').value=p[0]; $('out').value=p[1];
  CANDS=[];                       // 清掉上一属的卡片，避免串味
  try{ await loadCands(); }catch(e){}
  setTab('err');                  // 切到「确定错（审稿区）」，同时渲染全部发现
}
async function setTab(t){
  TAB=t;
  $('tabErr').classList.toggle('on',t==='err');
  $('tabUnk').classList.toggle('on',t==='unk');
  $('tabTrust').classList.toggle('on',t==='trust');
  $('tabNA').classList.toggle('on',t==='na');
  $('tabT').classList.toggle('on',t==='tree');
  $('tabCov').classList.toggle('on',t==='cov');
  $('trustview').style.display=(t==='trust')?'block':'none';
  $('naview').style.display=(t==='na')?'block':'none';
  $('treeview').style.display=(t==='tree')?'block':'none';
  $('covview').style.display=(t==='cov')?'block':'none';
  // 「全部发现」已并入审稿区：只在「确定错」这一页显示
  $('findview').style.display=(t==='err')?'block':'none';
  $('errview').style.display='none';
  const cardTab=(t==='err'||t==='unk');
  document.querySelector('.filters').style.display=cardTab?'flex':'none';
  $('cards').style.display=cardTab?'block':'none';
  $('empty').style.display='none';
  if(t==='trust')await renderTrust();
  else if(t==='na')await renderNA();
  else if(t==='cov')await renderCov();
  else if(t==='err')await renderFind();
  // ⚠️ 「确定错/不确定」靠 CANDS 渲染，而 CANDS 只在 loadCands() 里加载 ——
  //    页面加载时**不会自动拉**。于是刚刷新完页面切到「不确定」是空的，
  //    那个种看起来就还是「消失了」（本机实测踩到）。
  //    另两个选项卡（受信 / 暂无法判定）是按需 fetch 的，所以不缺。
  else if(cardTab){ if(CANDS.length)render(); else await loadCands(); }
}
async function runEnrich(link){
  const r=await fetch('/api/enrich',{method:'POST',
    headers:{'Content-Type':'application/json'},
    body:JSON.stringify({out:'out_audit', genus:'colletotrichum', link:!!link})});
  const j=await r.json();
  const lg=document.getElementById('findlog');
  lg.textContent = j.ok ? ('已启动补全（'+(j.link?'含号连接':'快档')+'）…') : ('失败：'+(j.msg||''));
  if(j.ok) pollEnrich();
}
async function pollEnrich(){
  const j=await fetch('/api/enrichresult',{cache:'no-store'}).then(x=>x.json());
  const lg=document.getElementById('findlog');
  lg.textContent=(j.log||[]).join('\n'); lg.scrollTop=1e9;
  if(j.running) setTimeout(pollEnrich,1500);
  else pollFind();
}
let FIND_ROWS = [], FIND_BY = {}, FIND_FILTER = 'all';
async function runFind(){
  const r=await fetch('/api/collect_findings',{method:'POST',
    headers:{'Content-Type':'application/json'},
    body:JSON.stringify({out:'out_audit', genus:'colletotrichum'})});
  const j=await r.json();
  document.getElementById('findlog').textContent = j.ok ? '已提交汇总…' : ('失败：'+(j.msg||''));
  if(j.ok) pollFind();
}
async function pollFind(){
  const j=await fetch('/api/findings?out=out_audit&genus=colletotrichum',{cache:'no-store'}).then(x=>x.json());
  const _fl=document.getElementById('findlog');
  _fl.textContent=(j.log||[]).join('\n');
  _fl.style.display = (j.log||[]).length ? 'block' : 'none';
  FIND_ROWS=j.rows||[]; FIND_BY=j.by||{};
  renderFindKinds();
  renderFindTable();
  if(j.running) setTimeout(pollFind,1500);
}
async function renderFind(){
  const j=await fetch('/api/findings?out=out_audit&genus=colletotrichum',{cache:'no-store'}).then(x=>x.json());
  FIND_ROWS=j.rows||[]; FIND_BY=j.by||{};
  renderFindKinds();
  const _fl2=document.getElementById('findlog');
  _fl2.textContent=(j.log||[]).join('\n');
  _fl2.style.display = (j.log||[]).length ? 'block' : 'none';
  renderFindTable();
}
function setFindFilter(k){ FIND_FILTER=k; renderFindTable(); }
function renderFindKinds(){
  const box=document.getElementById('findkinds');
  if(!box) return;
  const cnt={};
  // 归一：去掉括号里的参数、去掉「→CBS 112996」这类具体号，
  // 否则同一个类别会散成好几行（实测「经映射」被拆成 3 行）。
  const norm=k=>String(k||'—').replace(/（.*$/,'').replace(/→.*$/,'').trim();
  (FIND_ROWS||[]).forEach(r=>{const k=norm(r.类别);cnt[k]=(cnt[k]||0)+1;});
  const items=Object.entries(cnt).sort((a,b)=>b[1]-a[1]);
  const CLS={'号不存在':'#dc2626','书写错误':'#dc2626','并名冲突':'#dc2626',
             '异名/独立种分歧':'#dc2626','权威分歧':'#dc2626','CBS 记为别种':'#dc2626',
             '无目录':'#6b7280','未核':'#6b7280','标本号':'#6b7280',
             '疑似非机构编号':'#6b7280','NCBI 声明为模式号但机构目录查无':'#6b7280',
             'NCBI 记为模式材料':'#dc2626','经映射':'#6b7280',
             '异名（非错）':'#16a34a','基名属':'#16a34a','CBS 记为种下单元':'#16a34a'};
  box.innerHTML='<b>本工具能查出的类别（'+items.length+' 类，共 '+FIND_ROWS.length+' 条）：</b><br>'+
    items.map(([k,v])=>'<span style="display:inline-block;margin:2px 10px 2px 0">'+
      '<span style="color:'+(CLS[k]||'#111')+'">●</span> '+
      k.replace(/（.*$/,'')+' <b>'+v+'</b></span>').join('');
}
// 类别的白话解释：很多条目的「发现」列只是重复类别，读者看不出所以然。
// 只在【发现 === 类别】时替换为这段解释；数据本身一个字不改（统计仍按类别算）。
const FIND_WHY = {
  '无公开数据库':'该机构没有可机读的公开目录，无从核对',
  '未核':'号在机构目录中存在，但目录未标注它的模式地位',
  '疑似非机构编号':'形如「人名缩写+序号」的工作编号，疑似非机构保藏号',
  'NCBI 声明为模式号，机构无此记录':'GenBank 标为模式材料，但机构目录里查不到该号',
  '标本号':'这是标本馆号（如 CBS H-xxxx），不能拿培养物目录去核',
  '基名属':'记录用的是基名所属的属名（属名已变更），种加词未变',
  'CBS 记为别种':'CBS 目录把该号记在另一个种名下',
  'ICMP 记为别种':'ICMP 目录把该号记在另一个种名下',
  'CGMCC 记为非模式菌株':'CGMCC 目录明确标注该号不是模式菌株',
  '模式地位冲突':'GenBank 与机构目录对模式类型的说法不一致',
  '待核':'模式号归属待核'
};
function whyOf(r){
  const f=(r['发现']||'').trim(), c=(r['类别']||'').trim();
  if(f===c && FIND_WHY[c]) return FIND_WHY[c];
  return f;
}
// 号层（三列表逐号标记）在下方「逐号卡片」里已经一条一条显示，并且可以裁定。
// 本表默认只列【卡片里没有的】那部分：名字层 / 序列层 / 记录层 / 源库层。
// 勾上「含号层」才会列全（会与下方卡片重复）。
let SHOW_HAO=false;
function renderFindTable(){
  let rows = FIND_FILTER==='all' ? FIND_ROWS
             : FIND_ROWS.filter(function(r){ return (r['归类']||'')===FIND_FILTER; });
  const nBefore=rows.length;
  if(!SHOW_HAO) rows = rows.filter(function(r){ return (r['层级']||'')!=='号'; });
  const hid=nBefore-rows.length;
  const by=FIND_BY||{};
  const parts=Object.keys(by).sort().map(function(k){ return k+' '+by[k]; });
  document.getElementById('findsum').textContent =
    (FIND_ROWS.length? ('本表 '+rows.length+' 条'+
      (SHOW_HAO?'':'（另有 '+hid+' 条号层已在下方卡片，可逐条裁定）')+
      ' ｜ 全部 '+FIND_ROWS.length+' 条 ｜ '+parts.join(' ｜ '))
    : '还没有汇总结果，点「重新汇总」');
  if(!rows.length){ document.getElementById('findtable').innerHTML=''; return; }
  const cols=['序号','层级','类别','对象','发现','归类','来源'];
  let h='<table style="width:100%;border-collapse:collapse;font-size:11px">';
  h+='<tr style="background:#f3f4f6;position:sticky;top:0">'+
     cols.map(function(c){return '<th align=left>'+_escAuto(c)+'</th>';}).join('')+'</tr>';
  rows.forEach(function(r){
    const cls=r['归类']||'';
    const color = cls==='错值' ? '#b91c1c' : (cls==='矛盾' ? '#b91c1c' :
                  (cls==='不确定' ? '#92400e' : '#374151'));
    h+='<tr style="border-top:1px solid #e5e7eb">'+cols.map(function(c){
      const v=String(r[c]||'');
      const style=(c==='归类')?(' style="color:'+color+';font-weight:600"'):'';
      const raw=(c==='发现')?whyOf(r):v;
      const txt=(c==='发现'||c==='对象')?raw.slice(0,120):raw;
      const tip=(c==='发现'&&raw.length>40)?(' title="'+_escAuto(raw)+'"'):'';
      return '<td'+style+tip+'>'+_escAuto(txt)+'</td>';
    }).join('')+'</tr>';
  });
  document.getElementById('findtable').innerHTML=h+'</table>';
}
function _escAuto(s){
  return String(s==null?'':s).replace(/[&<>"]/g,function(c){
    return {'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;'}[c];});
}
async function runAuto(){
  const r=await fetch('/api/auto_verdict',{method:'POST',
    headers:{'Content-Type':'application/json'},
    body:JSON.stringify({pairs:document.getElementById('autopairs').value,
                         out:'tree_auto', min_cols:120, per_name:150})});
  const j=await r.json();
  document.getElementById('autolog').textContent = j.ok ? ('已提交 '+j.n+' 对…') : ('失败：'+(j.msg||''));
  if(j.ok) pollAuto();
}
async function pollAuto(){
  const r=await fetch('/api/autoresult',{cache:'no-store'});
  const j=await r.json();
  const lg=document.getElementById('autolog');
  lg.textContent=(j.log||[]).join('\n'); lg.scrollTop=1e9;
  renderAutoData(j);
  if(j.running) setTimeout(pollAuto,2000);
}
function renderAutoData(j){
  const rows=j.rows||[];
  document.getElementById('autosum').textContent =
    rows.length? (rows.length+' 对已出结果' + (j.running?'（仍在跑）':'')) : (j.running?'运行中…':'还没有结果');
  if(!rows.length){ document.getElementById('autotable').innerHTML=''; return; }
  let h='<table style="width:100%;border-collapse:collapse;font-size:12px">';
  h+='<tr style="background:#f3f4f6"><th align=left>争议对</th><th>层级</th><th>位点</th>'+
     '<th align=left>结论强度</th><th align=left>备注</th></tr>';
  rows.forEach(function(r){
    const k=Object.keys(r);
    h+='<tr style="border-top:1px solid #e5e7eb"><td>'+_escAuto(r[k[0]])+'</td>'+
       '<td align=center>'+_escAuto(r[k[1]])+'</td>'+
       '<td align=center>'+_escAuto(r[k[2]])+'</td>'+
       '<td>'+_escAuto(r[k[3]])+'</td>'+
       '<td class="tiny">'+_escAuto(r[k[4]])+'</td></tr>';
  });
  document.getElementById('autotable').innerHTML=h+'</table>';
}
async function runTree(){
  const body={a:$('ta').value,b:$('tb').value,workdir:$('twork').value,
    acc_range:$('trange').value,phylo_bin:$('tphylo').value,email:$('email').value};
  const r=await fetch('/api/tree',{method:'POST',headers:{'Content-Type':'application/json'},
    body:JSON.stringify(body)});
  const j=await r.json();
  if(!j.ok){alert(j.msg||'启动失败');return;}
  treePoll();
}
async function treePoll(){
  const r=await fetch('/api/treelog');const j=await r.json();
  $('tlog').textContent=j.log.join('\n');$('tlog').scrollTop=$('tlog').scrollHeight;
  if(j.running)setTimeout(treePoll,800);
  else{$('tlog').textContent+=j.rc===0?'\n\n✅ 结束':'\n\n（结束，rc='+j.rc+'）';}
}
async function loadErrors(){
  const r=await fetch('/api/errorlist?'+q());const j=await r.json();
  const rows=j.rows||[];
  $('errcnt').textContent='（'+rows.length+' 条）';
  if(!rows.length){$('errview').innerHTML='<div class="tiny">没有错值类候选（号不存在 / 未知机构代码 / CBS 记为别种 …）。</div>';return;}
  const src=`<div class="tiny" style="margin-bottom:6px;background:#f9fafb;border-left:3px solid #9ca3af;padding:4px 8px">
      数据来源：<code>${j.file||'—'}</code>${j.mtime?`　生成于 <b>${j.mtime}</b>`:''}
      　<span style="color:#92400e">（读磁盘上已有的表，不是本次运行算的）</span></div>`;
  $('errview').innerHTML=src+`<table class="et"><thead><tr>
      <th>种名</th><th>列</th><th>号</th><th>判定</th><th>机构记的名</th><th>接受名（MycoBank）</th><th>关系</th>
    </tr></thead><tbody>`+rows.map(x=>`<tr>
      <td>${x.species.replace('Colletotrichum ','C. ')}</td>
      <td>${x.column}</td><td><code>${x.item}</code></td>
      <td><span class="flag">${x.flags}</span></td>
      <td>${x.catalog_name||'—'}</td><td>${x.current_name||'—'}</td>
      <td>${x.verdict?`<span class="flag ${x.verdict.includes('非错')?'ok':'warn'}">${x.verdict}</span>`:'—'}</td>
    </tr>`).join('')+`</tbody></table>`;
}
window.addEventListener('load', fillRefs);
async function exportReviews(){
  const r=await fetch('/api/reviews?'+q());const j=await r.json();
  const tsv=['种名\t列\t号\t系统判定\t人工裁定\t备注\t时间']
    .concat((j.reviews||[]).map(x=>[x.species,x.column,x.item,x.system_flag,x.verdict,x.note,x.ts].join('\t'))).join('\n');
  const a=document.createElement('a');
  a.href=URL.createObjectURL(new Blob(['\ufeff'+tsv],{type:'text/tab-separated-values'}));
  a.download=CUR+'_裁定记录.tsv';a.click();
}

// ---- CBS 目录下载：点「⇩ 下载目录」→ 问范围 → 下载 → 自动把路径填回输入框 ----
async function dlCbs(){
  const g = prompt(
    '下载 CBS 目录\n\n要下载哪个属？\n（留空 = 下载【全量】，约 63,000 条 / 30 次请求 / 10–40 分钟）',
    $('genus').value);
  if (g === null) return;                       // 用户取消
  const genus = g.trim();
  let syn = false;
  if (genus) {
    syn = confirm('是否同时拉【同义属】？\n\n「确定」= 一起拉（慢，但能找回挂在有性型属名下的记录，' +
                  '如 CBS 113010 记作 Glomerella cingulata 而其实是炭疽菌）\n「取消」= 只按属名单词筛（快）');
  }
  const body = {genus: genus, email: $('email').value, out: $('out').value, no_syn: !syn};
  const r = await fetch('/api/cbs_fetch', {method:'POST',
      headers:{'Content-Type':'application/json'}, body: JSON.stringify(body)});
  const j = await r.json();
  if (!j.ok) { alert(j.msg || '启动失败'); return; }
  $('cbsstat').textContent = '⏳ 已开始下载' + (genus ? '（' + genus + '）' : '（全量）') + '……';
  dlPoll();
}
async function cbsStop(){ await fetch('/api/cbs_stop', {method:'POST'}); }
async function dlPoll(){
  const r = await fetch('/api/cbslog'); const j = await r.json();
  const last = ((j.log || []).slice(-1)[0] || '').replace(/</g, '&lt;').slice(0, 140);
  if (j.running) {
    $('cbsstat').innerHTML = '⏳ ' + last +
      '　<button class="gray" style="padding:1px 8px;font-size:11px" onclick="cbsStop()">■ 停止</button>';
    setTimeout(dlPoll, 1500);
    return;
  }
  if (j.out_tsv) {
    $('cbs').value = j.out_tsv;                 // ← 自动把路径填进上面的输入框
    $('cbsstat').innerHTML = '✅ 已下载并填入路径：<code>' + j.out_tsv + '</code>';
  } else {
    $('cbsstat').innerHTML = '下载已结束但没拿到产物路径 —— 看日志最后几行。';
  }
}

// ---------------- 详情弹窗 ----------------
// 为什么要这个：审稿人要能一眼看到「在哪里发现的、为什么判它是错」——
// 只给一个标记（如「CBS 记为别种」）等于把判断回给人瞎猜。
function showDetail(k){
  const c=CANDS.find(x=>encodeURIComponent(key(x))===k);
  if(!c){alert('找不到这条候选');return;}
  const rows=(c.detail||[]).map(([a,b])=>
    `<tr><td class="tiny" style="width:130px;color:#6b7280;vertical-align:top">${a}</td>`+
    `<td>${(b===null||b===undefined)?'—':String(b)}</td></tr>`).join('');
  $('dlgbody').innerHTML=
    `<div style="margin-bottom:6px"><span class="flag ${flagClass(c.system_flag)}">${shortBadge(c.badge||c.system_flag)}</span>`+
    `　<code>${c.item}</code></div>`+
    `<table class="et"><tbody>${rows}</tbody></table>`+
    `<div class="tiny" style="margin-top:10px;color:#6b7280">`+
    `裁定会写进 append-only 的 review_log.jsonl，永不静默删除。</div>`;
  $('dlg').style.display='block';
}
function closeDetail(){ $('dlg').style.display='none'; }

// ---------------- 受信菌株表 ----------------
let NA=null;
async function renderCov(){
  // 覆盖声明页：把「哪些种没扫全」摆出来。
  // 为什么必须单列：读者看到一张表，会默认它是"全查过了"。
  // 实际上记录数超限的种只扫了前 N 条 —— 不说，就是把没查完当成没问题。
  let rows=[];
  try{
    const r=await fetch('/api/coverage?'+q());const j=await r.json();
    rows=j.rows||[];
  }catch(e){ }
  let h=`<div class="tiny" style="margin-bottom:8px;background:#fffbeb;border-left:3px solid #f59e0b;padding:7px 9px">
    <b>覆盖声明</b>：下面这些种的 GenBank 记录条数超过本次扫描上限，
    <b>只扫了前面一部分</b>。因此它们的「号」可能不全 —— 这是<b>覆盖不足</b>，不是「查过且没问题」。
    <br>写报告时请把这类种单列，不要与完整扫描的种混在一起算错误率。
  </div>`;
  if(!rows.length){
    h+=`<div class="tiny">（本批没有记录数超限的种）</div>`;
  }else{
    h+=`<div class="tiny" style="margin-bottom:6px">共 <b>${rows.length}</b> 个种记录数超限</div>`;
    h+='<table class="et"><thead><tr><th>种名</th><th>覆盖情况</th></tr></thead><tbody>';
    for(const x of rows){
      h+=`<tr><td>${(x.species||'').replace('Colletotrichum ','')}</td><td class="tiny">${x.note||''}</td></tr>`;
    }
    h+='</tbody></table>';
  }
  $('covview').innerHTML=h;
}
async function renderNA(){
  if(NA===null){
    const r=await fetch('/api/undecidable?'+q());const j=await r.json();
    NA=j.rows||[];
  }
  const rows=NA;
  let h=`<div class="tiny" style="margin-bottom:8px;background:#fffbeb;border-left:3px solid #f59e0b;padding:7px 9px">
    <b>暂无法判定</b>：这些种<b>一个可核验的号都没有</b> —— 它们只有 GenBank 的 <code>/strain</code>（论文内部编号），
    或干脆是空列，本工具对它们<b>无从核验</b>。
    <br>⚠️ 把它们单列出来、而不是让它们消失，是因为「<b>无从检查</b>」与「<b>检查过、干净</b>」必须能分辨。
    <br>这不等于「没问题」—— 只是本工具答不了。
  </div>`;
  if(!rows.length){
    h+=`<div class="tiny">（这一类为空）</div>`;
  }else{
    h+=`<div class="tiny" style="margin-bottom:6px">共 <b>${rows.length}</b> 个种</div>`;
    h+='<table class="et"><thead><tr><th>种名</th><th>仅有的编号（未核验）</th><th>为什么无从判定</th></tr></thead><tbody>';
    rows.forEach(x=>{
      h+=`<tr><td>${x.species}</td><td><code>${x.strain_only||'—'}</code></td>`+
         `<td class="tiny">${(x.reason||'').slice(0,320)}</td></tr>`;
    });
    h+='</tbody></table>';
  }
  $('naview').innerHTML=h;
}

async function renderTrust(){  if(TRUSTED===null){
    const r=await fetch('/api/trusted?'+q());const j=await r.json();
    TRUSTED=j.rows||[];TSRC={file:j.file,mtime:j.mtime,sources:j.sources||[]};
  }
  const rows=TRUSTED;
  const note=`<div class="tiny" style="margin-bottom:8px;background:#f9fafb;border-left:3px solid #9ca3af;padding:5px 9px">
      <b>受信菌株表</b>：只放<b>没有任何问题标记</b>的号 —— 即经过检查、与外部目录/登记表一致的那些。
      <br><b>「受信」的含义仅限于「号与外部目录一致」</b>：不等于这株菌真是模式菌株，也不等于序列属于该种。
      <br>判定规则：① 同一个号在<b>任一列</b>只要出现「错值」或「不确定」类标记就不进本表（一票否决）；
      只有「非错」类标记（异名／基名属）或无标记的才算受信；
      ② 只收<b>真正的保藏号</b> —— GenBank <code>/strain</code> 里的实验室编号没有机构前缀，
      不在本工具核对范围内，一律不收（把「没检查」当「受信」是错的）。
      ${TSRC.mtime?`<br>来源：<code>${TSRC.file||'—'}</code>（生成于 ${TSRC.mtime}）`:''}
      ${(TSRC.sources&&TSRC.sources.length)?`　参与判定的源：${TSRC.sources.join('、')}`:''}
    </div>
    <div style="margin-bottom:8px"><button onclick="exportTrusted()">⤓ 导出受信菌株表（TSV）</button>
      <span class="tiny">　共 ${rows.length} 行</span></div>`;
  if(!rows.length){
    $('trustview').innerHTML=note+'<div class="tiny">没有受信菌株 —— 要么所有号都带标记，要么还没跑流水线。</div>';
    return;
  }
  $('trustview').innerHTML=note+`<table class="et"><thead><tr>
      <th>种名</th><th>模式菌株</th><th>菌株号</th><th>受信依据</th></tr></thead><tbody>`+
    rows.map(x=>`<tr><td>${x.species}</td><td><code>${x.type||'—'}</code></td>`+
      `<td><code>${x.strain||'—'}</code></td>`+
      `<td class="tiny">${(x.why_type||'')+(x.why_type&&x.why_strain?'　':'')+(x.why_strain||'')||
        '<span style="color:#6b7280">检查通过，无标记</span>'}</td></tr>`).join('')+
    `</tbody></table>`;
}
// 页面打开即拉一次「结论速览」——不联网、不依赖任何选择
try{ loadSummary(); }catch(e){}
function exportTrusted(){
  if(!TRUSTED)return;
  const tsv=['种名\t模式菌株\t菌株号']
    .concat(TRUSTED.map(x=>[x.species,x.type,x.strain].join('\t'))).join('\n');
  const a=document.createElement('a');
  a.href=URL.createObjectURL(new Blob(['\ufeff'+tsv],{type:'text/tab-separated-values'}));
  a.download=($('genus').value||'genus')+'_受信菌株表.tsv';a.click();
}
</script></body></html>
"""


class H(BaseHTTPRequestHandler):
    def _send(self, code, body, ctype="application/json; charset=utf-8"):
        b = body if isinstance(body, bytes) else str(body).encode("utf-8")
        self.send_response(code)
        # ⚠️ 禁止浏览器缓存本页。不加这几行时，改了前端后用户按 F5
        #    仍可能拿到旧页面 —— 本机实测症状：新加的「暂无法判定」选项卡
        #    在服务端已存在（页面 HTML 里查得到），但界面里没有这个按钮，
        #    于是那个种看起来还是「凭空消失」。
        self.send_header("Cache-Control", "no-store, no-cache, must-revalidate, max-age=0")
        self.send_header("Pragma", "no-cache")
        self.send_header("Expires", "0")
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(b)))
        self.end_headers()
        self.wfile.write(b)

    def _json(self, obj, code=200):
        self._send(code, json.dumps(obj, ensure_ascii=False))

    def _body(self):
        n = int(self.headers.get("Content-Length") or 0)
        try:
            return json.loads(self.rfile.read(n) or b"{}")
        except Exception:
            return {}

    def do_GET(self):
        u = urlparse(self.path)
        qs = {k: v[0] for k, v in parse_qs(u.query).items()}
        if u.path in ("/", "/index.html"):
            return self._send(200, PAGE, "text/html; charset=utf-8")
        if u.path == "/api/log":
            with LOCK:
                lines = list(LOG)
            st = "运行中…" if STATE["running"] else ("完成" if STATE["rc"] == 0 else
                                                 ("已停止" if STATE["rc"] == 130 else "空闲"))
            return self._json({"log": lines or ["（暂无日志）"], "running": STATE["running"],
                               "status": st, "rc": STATE["rc"], "deps": DEPS})
        if u.path == "/api/candidates":
            od = qs.get("out", "out_audit")
            p, cards = read_candidates(od, qs.get("genus", ""), sources_for(od))
            return self._json({"file": p, "cards": cards, "mtime": _mtime(p),
                               "classes": flag_classes(),
                               "progress": read_progress(od, qs.get("genus", ""))})
        if u.path == "/api/cbslog":
            with CLOCK:
                lines = list(CLOG)
            return self._json({"log": lines or ["（还没跑）"], "running": CBS_STATE["running"],
                               "rc": CBS_STATE["rc"], "mode": CBS_STATE["mode"],
                               "out_tsv": CBS_STATE.get("out_tsv")})
        if u.path == "/api/trusted":
            od = qs.get("out", "out_audit")
            p, rows = read_trusted(od, qs.get("genus", ""))
            return self._json({"file": p, "rows": rows, "mtime": _mtime(p),
                               "sources": sorted(sources_for(od).keys())})
        if u.path == "/api/undecidable":
            od = qs.get("out", "out_audit")
            p, rows = read_undecidable(od, qs.get("genus", ""))
            return self._json({"file": p, "rows": rows, "mtime": _mtime(p)})
        if u.path == "/api/summary":
            # 三属结论速览：直接读 results/summary.json。
            # 为什么要它：论文读者多半不写代码，打开网页只想看「结论是什么」，
            # 而审稿区是逐条卡片、看不出全局数字。
            sp = os.path.join(REPO, "results", "summary.json")
            try:
                d = json.load(io.open(sp, encoding="utf-8"))
            except Exception as e:
                d = {"error": "{}: {}".format(type(e).__name__, str(e)[:80])}
            return self._json(d)
        if u.path == "/api/coverage":
            # 覆盖声明：哪些种因为 GenBank 记录太多，本次只扫了一部分。
            # 为什么要单开一个接口：这件事决定「这张表可信到什么程度」，
            # 但三列表每格必须是「号＋[标记]」，塞不进去 —— 所以单列一处展示。
            od = qs.get("out", "out_audit")
            if not os.path.isabs(od):
                od = os.path.join(REPO, od)
            g = (qs.get("genus") or "").lower()
            cand = []
            if g:
                cand.append(os.path.join(od, "{}_覆盖声明.tsv".format(g)))
            cand.append(os.path.join(od, "覆盖声明.tsv"))
            p = next((x for x in cand if os.path.exists(x)), "")
            rows = []
            if p:
                try:
                    for l in io.open(p, encoding="utf-8-sig").read().splitlines()[1:]:
                        c = l.split("\t")
                        if c and c[0].strip():
                            rows.append({"species": c[0].strip(),
                                         "note": (c[1] if len(c) > 1 else "").strip()})
                except Exception:
                    rows = []
            return self._json({"file": p, "rows": rows, "mtime": _mtime(p)})
        if u.path == "/api/refs":
            return self._json(REFS)
        if u.path == "/api/enrichresult":
            return self._json({"running": ENRICH_STATE["running"], "rc": ENRICH_STATE["rc"],
                               "log": ENRICH_LOG[-300:]})
        if u.path == "/api/findings":
            od = qs.get("out", "out_audit")
            if not os.path.isabs(od):
                od = os.path.join(REPO, od)
            p, rows, cnt = read_findings(od, qs.get("genus", "colletotrichum"))
            return self._json({"file": p, "rows": rows, "by": cnt,
                               "running": FIND_STATE["running"], "log": FIND_LOG[-120:]})
        if u.path == "/api/autoresult":
            return self._json({"running": AUTO_STATE["running"], "rc": AUTO_STATE["rc"],
                               "log": AUTO_LOG[-300:], "rows": AUTO_ROWS})
        if u.path == "/api/speciesresult":
            return self._json({"running": SPEC_STATE["running"], "rc": SPEC_STATE["rc"],
                               "log": SPEC_LOG[-300:], "rows": SPEC_ROWS})
        if u.path == "/api/linkresult":
            return self._json({"running": LINK_STATE["running"], "rc": LINK_STATE["rc"],
                               "log": LINK_LOG[-200:],
                               "groups": LINK_RESULT.get("groups") or [],
                               "singletons": LINK_RESULT.get("singletons") or [],
                               "edges": LINK_RESULT.get("edges") or [],
                               "stats": LINK_RESULT.get("stats") or {}})
        if u.path == "/api/revlog":
            return self._json({"log": REV_LOG[-250:], "running": REV_STATE["running"],
                               "rc": REV_STATE["rc"]})
        if u.path == "/api/reverse":
            tsv = os.path.join(REV_OUT, "反向审计.tsv")
            rows = []
            if os.path.exists(tsv):
                for line in io.open(tsv, encoding="utf-8-sig").read().splitlines()[1:]:
                    c = line.split("\t")
                    if len(c) >= 7:
                        rows.append({"species": c[0], "number": c[1], "status": c[2],
                                     "hits": c[3], "verified": c[4],
                                     "verdict": c[5], "note": c[6]})
            from collections import Counter as _C
            return self._json({"ok": True, "rows": rows, "total": len(rows),
                               "summary": dict(_C(r["verdict"] for r in rows))})
        if u.path == "/api/treelog":
            with LOCK:
                lines = list(TLOG)
            return self._json({"log": lines or ["（还没跑）"], "running": TREE_STATE["running"],
                               "rc": TREE_STATE["rc"], "pair": TREE_STATE["pair"]})
        if u.path == "/api/errorlist":
            od = qs.get("out", "out_audit")
            p, rows = read_errors(od, qs.get("genus", ""), sources_for(od))
            return self._json({"file": p, "rows": rows, "mtime": _mtime(p)})
        if u.path == "/api/reviews":
            return self._json({"reviews": read_reviews(qs.get("out", "out_audit"))})
        return self._send(404, "not found", "text/plain; charset=utf-8")

    def do_POST(self):
        u = urlparse(self.path)
        d = self._body()
        if u.path == "/api/run":
            if STATE["running"]:
                return self._json({"ok": False, "msg": "已在运行中"})
            if not (d.get("genus") or "").strip():
                return self._json({"ok": False, "msg": "属名必填"})
            t = threading.Thread(target=worker, daemon=True, kwargs=dict(
                genus=d["genus"].strip(), out_dir=(d.get("out") or "out_audit").strip(),
                limit=(d.get("limit") or "").strip() or None,
                framework=(d.get("framework") or "").strip(),
                species_list=(d.get("species") or "").strip(),
                cbs=(d.get("cbs") or "").strip(), email=(d.get("email") or "").strip(),
                cfg=None, mycobank=(d.get("mycobank") or "").strip(),
                biocol=(d.get("biocol") or "").strip(),
                skip_done=bool(d.get("skip_done"))))
            t.start()
            return self._json({"ok": True})
        if u.path == "/api/tree":
            if TREE_STATE["running"]:
                return self._json({"ok": False, "msg": "建树已在运行"})
            a = (d.get("a") or "").strip()
            b = (d.get("b") or "").strip()
            if not a or not b:
                return self._json({"ok": False, "msg": "要填两个保藏号（如 ICMP:18578 / ICMP:18642）"})
            threading.Thread(target=tree_worker, daemon=True, kwargs=dict(
                workdir=(d.get("workdir") or "tree_work").strip(),
                pair=[a, b],
                acc_range=(d.get("acc_range") or "JX009400:JX010500").strip(),
                email=(d.get("email") or os.environ.get("NCBI_EMAIL") or "").strip(),
                phylo_bin=(d.get("phylo_bin") or "/home/yelu/miniforge3/envs/phylo/bin").strip(),
            )).start()
            return self._json({"ok": True})
        if u.path == "/api/enrich":
            if ENRICH_STATE["running"]:
                return self._json({"ok": False, "msg": "补全已在运行"})
            od = (d.get("out") or "out_audit").strip()
            if not os.path.isabs(od):
                od = os.path.join(REPO, od)
            threading.Thread(target=enrich_worker, daemon=True,
                             kwargs=dict(out_dir=od,
                                         genus=(d.get("genus") or "colletotrichum").strip(),
                                         link=bool(d.get("link")))).start()
            return self._json({"ok": True, "link": bool(d.get("link"))})

        if u.path == "/api/collect_findings":
            if FIND_STATE["running"]:
                return self._json({"ok": False, "msg": "汇总已在运行"})
            od = (d.get("out") or "out_audit").strip()
            if not os.path.isabs(od):
                od = os.path.join(REPO, od)
            threading.Thread(target=collect_worker, daemon=True,
                             kwargs=dict(out_dir=od,
                                         genus=(d.get("genus") or "colletotrichum").strip())
                             ).start()
            return self._json({"ok": True})

        if u.path == "/api/auto_verdict":
            if AUTO_STATE["running"]:
                return self._json({"ok": False, "msg": "自动裁决已在运行"})
            txt = (d.get("pairs") or "").strip()
            pairs = []
            for line in re.split(r"[\n;]+", txt):
                line = line.strip()
                if not line:
                    continue
                f = re.split(r"[\s,]+", line)
                if len(f) >= 2:
                    pairs.append((f[0], f[1]))
            if not pairs:
                return self._json({"ok": False, "msg": "请每行填一对种加词，如 celtidis karsti"})
            try:
                min_cols = int(d.get("min_cols") or 120)
                per_name = int(d.get("per_name") or 150)
            except Exception:
                min_cols, per_name = 120, 150
            threading.Thread(target=auto_worker, daemon=True,
                             kwargs=dict(pairs=pairs,
                                         out_root=os.path.join(REPO, d.get("out") or "tree_auto"),
                                         min_cols=min_cols, per_name=per_name)).start()
            return self._json({"ok": True, "n": len(pairs)})

        if u.path == "/api/species_status":
            if SPEC_STATE["running"]:
                return self._json({"ok": False, "msg": "四方核验已在运行"})
            names = [x.strip() for x in re.split(r"[\n;]+", d.get("names") or "") if x.strip()]
            names = [x if " " in x else "Colletotrichum " + x for x in names]
            if not names:
                return self._json({"ok": False, "msg": "请填种名（每行一个）"})
            out_dir = (d.get("out") or "out_rerun").strip()
            if not os.path.isabs(out_dir):
                out_dir = os.path.join(REPO, out_dir)
            threading.Thread(target=spec_worker, daemon=True,
                             kwargs=dict(names=names, out_dir=out_dir)).start()
            return self._json({"ok": True, "n": len(names)})

        if u.path == "/api/annotate_links":
            if LINK_STATE["running"]:
                return self._json({"ok": False, "msg": "已有任务在运行"})
            out_dir = (d.get("out") or "out_rerun").strip()
            if not os.path.isabs(out_dir):
                out_dir = os.path.join(REPO, out_dir)
            try:
                budget = int(d.get("icmp_budget") or 200)
            except Exception:
                budget = 200
            threading.Thread(target=annotate_worker, daemon=True,
                             kwargs=dict(out_dir=out_dir, icmp_budget=budget)).start()
            return self._json({"ok": True})
        if u.path == "/api/link":
            if LINK_STATE["running"]:
                return self._json({"ok": False, "msg": "号连接已在运行"})
            out_dir = (d.get("out") or "out_rerun").strip()
            if not os.path.isabs(out_dir):
                out_dir = os.path.join(REPO, out_dir)
            mode = (d.get("mode") or "manual").strip()
            if mode == "batch":
                nums = _numbers_from_table(out_dir)
            else:
                nums = [x.strip() for x in re.split(r"[\n,;]+", d.get("numbers") or "")
                        if x.strip()]
            if not nums:
                return self._json({"ok": False, "msg": "没有可处理的号（手动模式请贴号；批量模式需该目录已有三列表）"})
            try:
                budget = int(d.get("icmp_budget") or 200)
            except Exception:
                budget = 200
            threading.Thread(target=link_worker, daemon=True,
                             kwargs=dict(numbers=nums, out_dir=out_dir,
                                         icmp_budget=budget)).start()
            return self._json({"ok": True, "n": len(nums)})
        if u.path == "/api/reverse_run":
            if REV_STATE["running"]:
                return self._json({"ok": False, "msg": "反向审计已在运行"})
            threading.Thread(target=reverse_worker, daemon=True,
                             kwargs=dict(email=(d.get("email") or os.environ.get("NCBI_EMAIL") or "").strip())).start()
            return self._json({"ok": True})
        if u.path == "/api/treestop":
            TREE_STATE["stop"] = True
            tree_log("⛔ 已请求停止…")
            return self._json({"ok": True})
        if u.path == "/api/cbs_fetch":
            if CBS_STATE["running"]:
                return self._json({"ok": False, "msg": "下载已在运行中（先停止或等它跑完）"})
            threading.Thread(target=cbs_worker, daemon=True, kwargs=dict(
                genus=(d.get("genus") or "").strip(),
                out_dir=(d.get("out") or "out").strip(),
                email=(d.get("email") or os.environ.get("NCBI_EMAIL") or "").strip(),
                page=(d.get("page") or "").strip() or None,
                sleep=(d.get("sleep") or "").strip() or None,
                no_syn=bool(d.get("no_syn")),
                syn_min=(d.get("syn_min") or "").strip() or None)).start()
            return self._json({"ok": True})
        if u.path == "/api/cbs_stop":
            p = CBS_STATE.get("proc")
            if p:
                try:
                    p.terminate()
                except Exception:
                    pass
            cbs_log("⛔ 已请求停止。已拉的页有缓存，重跑会接着走。")
            return self._json({"ok": True})
        if u.path == "/api/stop":
            STATE["stop"] = True
            logln("⛔ 已请求停止（当前种跑完即停）…")
            return self._json({"ok": True})
        if u.path == "/api/review":
            try:
                append_review((d.get("out") or "out_audit").strip(), {
                    "genus": d.get("genus", ""), "species": d.get("species", ""),
                    "column": d.get("column", ""), "item": d.get("item", ""),
                    "system_flag": d.get("system_flag", ""), "verdict": d.get("verdict", ""),
                    "reviewer": d.get("reviewer", ""), "note": d.get("note", "")})
                return self._json({"ok": True})
            except Exception as e:
                return self._json({"ok": False, "msg": str(e)})
        return self._send(404, "not found", "text/plain; charset=utf-8")

    def log_message(self, *a):
        pass          # 静音访问日志


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--port", type=int, default=8765)
    ap.add_argument("--host", default="127.0.0.1")
    a = ap.parse_args()
    srv = ThreadingHTTPServer((a.host, a.port), H)
    url = "http://{}:{}/".format(a.host, a.port)
    print("=" * 60)
    print("类型材料审计 · 流水线 + 审稿系统")
    print("  打开：{}".format(url))
    print("  停止：Ctrl+C")
    print("  依赖：tkinter={}  biopython={}".format(
        DEPS.get("tkinter"), DEPS.get("biopython") or "缺失（页面上会提示修法）"))
    print("  解释器：{}".format(DEPS["python"]))
    print("=" * 60)
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        print("\n已停止服务。")


if __name__ == "__main__":
    main()


def run_py_script(script, args, log, label=None):
    """在【进程内】执行另一个脚本的 main()。

    ⚠️ 为什么不能 subprocess + sys.executable：
       程序被 PyInstaller 打包后，`sys.executable` 指向 **exe 本身**，不是 Python。
       于是 Popen([sys.executable, "xxx.py", ...]) 实际上是「又启动一次 exe」——
       它检测到服务已在运行、开个浏览器、退出码 0。
       结果：功能显示成功、实际什么都没做（实测「一键补全」就是这样假成功的）。
    """
    import contextlib, importlib.util, sys as _sys
    p = os.path.join(HERE, script)
    label = label or script
    log("{} …".format(label))
    t0 = time.time()

    class _W:
        def write(self, s):
            if s and s.strip():
                log("   " + s.rstrip()[:400])
            return len(s)
        def flush(self):
            pass

    old_argv = _sys.argv
    try:
        spec = importlib.util.spec_from_file_location(
            "_inproc_" + os.path.splitext(os.path.basename(script))[0], p)
        m = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(m)
        _sys.argv = [p] + [str(x) for x in (args or [])]
        # ⚠️ stderr 也要捕获：argparse 的错误、traceback 都走 stderr，
        #    只重定向 stdout 会让「因参数错而退出」看起来像「跑完但没输出」（实测 rc=2 却零日志）。
        with contextlib.redirect_stdout(_W()), contextlib.redirect_stderr(_W()):
            rc = m.main()
        log("  → 完成（{}，{:.0f}s）".format(rc, time.time() - t0))
        return rc
    except SystemExit as e:
        log("  → 退出（{}，{:.0f}s）".format(e.code, time.time() - t0))
        return e.code
    except Exception as e:
        log("  ❌ {}: {}".format(type(e).__name__, str(e)[:200]))
        return 1
    finally:
        _sys.argv = old_argv
