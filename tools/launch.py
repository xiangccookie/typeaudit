# -*- coding: utf-8 -*-
# --- console utf8 bootstrap（Windows 控制台默认 GBK，输出 ↔ ⚠️ 🔴 会崩；详见 src/typeaudit/_console.py）
import sys as _sys
for _s in (_sys.stdout, _sys.stderr):
    try:
        _s.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass
"""
launch.py —— 类型校对台的桌面启动器。

【它解决什么问题】
    原来要开网页版得先开终端、cd 到仓库、敲 `python tools/pipeline_web.py`。
    做成 exe 之后，双击就行：起服务 + 自动打开浏览器。

【为什么 exe 还要找仓库】
    这个工具不是独立程序 —— 它必须读仓库里的东西：
        src/typeaudit/   代码
        config/          数据源配置
        data/            保藏目录、机构登记表、语料
        out_*/           各轮输出
    所以 exe 里烧不进这些，只能「找到仓库再启动」。查找顺序：
        ① 环境变量 TYPEPAUDIT_REPO
        ② exe（或本脚本）同目录下的 repo.txt，里面写一行仓库路径
        ③ 从 exe 所在目录逐级往上找，看哪一级有 src/typeaudit/pipeline.py
        ④ 几个已知位置兜底

【用法】
    直接双击 exe（开了浏览器就说明成功）
    也可以带参数：
        typeaudit_web.exe --port 9000       换端口
        typeaudit_web.exe --no-browser      只起服务、不开浏览器（测试用）

【日志】
    所有启动信息写到 <仓库>/out_web/_launcher.log。
    why：exe 是「无控制台」模式，出错时不能弹黑窗，只能靠日志排查。
"""
import io
import os
import socket
import sys
import threading
import time
import traceback
import urllib.request


def _log_path(repo):
    d = os.path.join(repo, "out_web")
    try:
        os.makedirs(d, exist_ok=True)
    except Exception:
        d = repo
    return os.path.join(d, "_launcher.log")


def _write(logfile, msg):
    """写日志（同时尽量打到控制台，便于命令行调试）。"""
    line = "[{}] {}\n".format(time.strftime("%Y-%m-%d %H:%M:%S"), msg)
    try:
        with io.open(logfile, "a", encoding="utf-8", newline="\n") as f:
            f.write(line)
    except Exception:
        pass
    try:
        sys.stdout.write(line)
        sys.stdout.flush()
    except Exception:
        pass


def find_repo():
    """按上面写的四个顺序找仓库根目录，返回绝对路径或 None。"""
    marks = os.path.join("src", "typeaudit", "pipeline.py")

    # ① 环境变量
    p = os.environ.get("TYPEPAUDIT_REPO", "").strip()
    if p and os.path.exists(os.path.join(p, marks)):
        return os.path.abspath(p)

    # ② exe / 脚本同目录下的 repo.txt
    here = os.path.dirname(os.path.abspath(sys.executable if getattr(sys, "frozen", False)
                                          else __file__))
    for cand in (here, os.path.dirname(here)):
        rp = os.path.join(cand, "repo.txt")
        if os.path.exists(rp):
            p = io.open(rp, encoding="utf-8-sig").read().strip().splitlines()
            if p and os.path.exists(os.path.join(p[0].strip(), marks)):
                return os.path.abspath(p[0].strip())

    # ③ 从所在目录逐级往上找
    cur = here
    for _ in range(6):
        if os.path.exists(os.path.join(cur, marks)):
            return cur
        parent = os.path.dirname(cur)
        if parent == cur:
            break
        cur = parent

    # ④ 已知位置兜底
    for p in (r"C:\Users\Administrator\Desktop\sol\奇怪\colleto-type-audit",
              r"C:\Users\Administrator\Desktop\colleto-type-audit"):
        if os.path.exists(os.path.join(p, marks)):
            return os.path.abspath(p)
    return None


def already_running(url, timeout=2.0):
    """端口上已经有服务了吗？有就不重复起，直接开浏览器。"""
    try:
        with urllib.request.urlopen(url, timeout=timeout) as r:
            return r.status == 200
    except Exception:
        return False


def port_free(host, port):
    """端口空闲吗（用来提前给出友好提示，而不是让服务起不来才报错）。"""
    s = socket.socket()
    try:
        s.bind((host, port))
        return True
    except OSError:
        return False
    finally:
        s.close()


def main():
    # 先解析参数（不依赖 argparse，减少打包变量）
    port, host, open_browser = 8765, "127.0.0.1", True
    argv = sys.argv[1:]
    i = 0
    while i < len(argv):
        a = argv[i]
        if a in ("--port",) and i + 1 < len(argv):
            port = int(argv[i + 1]); i += 2; continue
        if a in ("--host",) and i + 1 < len(argv):
            host = argv[i + 1]; i += 2; continue
        if a in ("--no-browser",):
            open_browser = False; i += 1; continue
        i += 1

    repo = find_repo()
    if not repo:
        # 连日志都写不了，只能退到临时目录
        tmp = os.path.join(os.environ.get("TEMP", "."), "_launcher.log")
        _write(tmp, "❌ 找不到仓库（需要含 src/typeaudit/pipeline.py 的目录）。"
                    "请在 exe 旁边放一个 repo.txt，里面写仓库路径。")
        return 2

    log = _log_path(repo)
    _write(log, "启动器开始 | repo={} | port={} | frozen={}".format(
        repo, port, bool(getattr(sys, "frozen", False))))

    os.chdir(repo)
    for p in (repo, os.path.join(repo, "tools"), os.path.join(repo, "src")):
        if p not in sys.path:
            sys.path.insert(0, p)

    url = "http://{}:{}/".format(host, port)

    # 已经在跑 → 只开浏览器
    if already_running(url):
        _write(log, "服务已在运行，直接打开浏览器：{}".format(url))
        if open_browser:
            import webbrowser
            webbrowser.open(url)
        time.sleep(1.5)
        return 0

    if not port_free(host, port):
        _write(log, "❌ 端口 {} 被占用，但访问不通。换端口：--port 9000".format(port))
        return 3

    try:
        import pipeline_web as PW
    except Exception:
        _write(log, "❌ 导入 pipeline_web 失败:\n" + traceback.format_exc())
        return 4

    import http.server
    try:
        srv = PW.ThreadingHTTPServer((host, port), PW.H)
    except Exception:
        _write(log, "❌ 起服务失败:\n" + traceback.format_exc())
        return 5

    th = threading.Thread(target=srv.serve_forever, daemon=True)
    th.start()
    _write(log, "服务已起：{}".format(url))

    # 等端口真的能响应，再开浏览器（否则会看到「无法访问」）
    for _ in range(40):
        if already_running(url, timeout=0.5):
            break
        time.sleep(0.25)
    else:
        _write(log, "⚠️ 服务起来了但端口 10 秒内没响应，仍尝试打开浏览器")

    if open_browser:
        import webbrowser
        webbrowser.open(url)
        _write(log, "已请求打开浏览器")

    # 常驻：Ctrl+C 或关掉进程即停
    try:
        while True:
            time.sleep(1.0)
    except KeyboardInterrupt:
        _write(log, "收到中断，正在停止服务")
    finally:
        try:
            srv.shutdown()
        except Exception:
            pass
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except SystemExit:
        raise
    except Exception:
        try:
            rp = find_repo() or os.environ.get("TEMP", ".")
            _write(_log_path(rp) if rp else os.path.join(rp, "_launcher.log"),
                   "❌ 未捕获异常:\n" + traceback.format_exc())
        except Exception:
            pass
        sys.exit(9)
