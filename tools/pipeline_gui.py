#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""属级类型材料审计 · 图形小窗口（tkinter，零额外依赖）

用法
----
  双击本文件，或：  python tools/pipeline_gui.py

它干什么
--------
把 `typeaudit.pipeline`（给属名 → 出三列表）包成一个窗口：
  填属名 → 点「开始」→ 右侧实时看日志 → 完成后点「打开结果表」。

为什么要后台线程
----------------
一次审计要联网逐种取号，**几分钟到数小时**。若直接在主线程跑，窗口会假死、
点不动、也没法中止。所以：重活放 `threading.Thread`，
日志经 `queue.Queue` 回主线程刷新（tkinter 只能在主线程改控件）。

两个容易忘的坑（代码里都处理了）
--------------------------------
* `typeaudit` 未安装时，要能把仓库的 `src/` 塞进 sys.path；
* 用户点了「停止」后，正在跑的 `run()` 不会立刻中断——**只能标记，
  由流水线逐个种之间检查**（故停止是"当前种跑完后停"，不是硬杀）。
"""
import os
import queue
import subprocess
import sys
import threading
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

try:
    sys.stdout.reconfigure(encoding="utf-8")
    sys.stderr.reconfigure(encoding="utf-8")
except Exception:
    pass

# ---- 让脚本在"未安装"时也能跑 ----
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
SRC = os.path.join(REPO, "src")
if os.path.isdir(SRC) and SRC not in sys.path:
    sys.path.insert(0, SRC)

STOP = threading.Event()          # 用户点「停止」置位；流水线在种之间检查


def launch(app):
    """后台线程里跑流水线。"""
    try:
        from typeaudit import pipeline as P

        # 让流水线能感知「停止」：包一层 acquire_species，逐个种之间查 STOP
        orig_collect = P.collect

        def collect_guard(*a, **kw):
            if STOP.is_set():
                raise KeyboardInterrupt("用户停止")
            return orig_collect(*a, **kw)

        P.collect = collect_guard

        cfg = app.cfg_var.get().strip() or None
        args = dict(
            genus=app.genus_var.get().strip(),
            out_dir=app.out_var.get().strip() or "out",
            config=cfg,
            framework=app.fw_var.get().strip() or None,
            species_list=app.sp_var.get().strip() or None,
            limit=(None if app.nolimit_var.get() else (int(app.limit_var.get() or 0) or None)),
            email=app.email_var.get().strip() or None,
        )
        app.log("=" * 60)
        app.log("属名：{}".format(args["genus"] or "（未填！）"))
        app.log("输出：{}".format(args["out_dir"]))
        app.log("=" * 60)
        if not args["genus"]:
            app.log("❌ 请先填属名")
            return
        P.run(**args)
        app.log("\n✅ 完成。")
        app.done()
    except KeyboardInterrupt as e:
        app.log("\n⛔ 已停止：{}".format(e))
    except Exception as e:
        import traceback
        app.log("\n❌ 出错：{}".format(e))
        app.log(traceback.format_exc()[-1500:])
    finally:
        app.running = False
        app.sync_buttons()


class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("炭疽菌属 · 类型材料审计")
        self.geometry("980x680")
        self.q = queue.Queue()
        self.running = False
        self.last_out = None
        self._build()
        self.after(120, self._drain)
        self._selfcheck()

    def _selfcheck(self):
        """开窗就报依赖状态——别让用户点「开始」才撞 ImportError。
        实测：系统 Python 常缺 biopython（tkinter 一般自带）。"""
        self.log("依赖自检：")
        self.log("  python: {}".format(sys.executable))
        try:
            import tkinter as _tk
            self.log("  [OK] tkinter {}".format(_tk.TkVersion))
        except Exception as e:
            self.log("  [缺] tkinter -> {}".format(e))
        try:
            import Bio
            self.log("  [OK] biopython {}".format(Bio.__version__))
        except Exception:
            self.log("  [缺] biopython --- 点「开始」会 ImportError。装法与命令：")
            self.log('       "{}" -m pip install biopython'.format(sys.executable))
            self.log("       或用已带 biopython 的解释器打开本文件。")
        self.log("")

    # ---------------- 界面 ----------------
    def _row(self, parent, r, label, var, width=42, browse=None, tip=""):
        tk.Label(parent, text=label, anchor="e", width=12).grid(row=r, column=0, sticky="e", padx=(8, 4), pady=3)
        e = tk.Entry(parent, textvariable=var, width=width)
        e.grid(row=r, column=1, sticky="we", pady=3)
        if browse:
            tk.Button(parent, text="浏览", command=browse).grid(row=r, column=2, padx=4)
        else:
            tk.Label(parent, text=tip, fg="#777").grid(row=r, column=2, sticky="w", padx=4)
        return e

    def _build(self):
        top = ttk.LabelFrame(self, text="参数")
        top.pack(fill="x", padx=8, pady=(8, 4))
        top.columnconfigure(1, weight=1)

        self.genus_var = tk.StringVar(value="Colletotrichum")
        self.limit_var = tk.StringVar(value="10")
        self.nolimit_var = tk.BooleanVar(value=False)
        self.fw_var = tk.StringVar()
        self.sp_var = tk.StringVar()
        self.cbs_var = tk.StringVar()
        self.out_var = tk.StringVar(value=os.path.join(os.getcwd(), "out_audit"))
        self.email_var = tk.StringVar(value=os.environ.get("NCBI_EMAIL", ""))

        self._row(top, 0, "属名", self.genus_var, tip="必填，如 Colletotrichum")

        f = tk.Frame(top); f.grid(row=1, column=1, columnspan=2, sticky="we", pady=3)
        tk.Entry(f, textvariable=self.limit_var, width=8).pack(side="left")
        tk.Label(f, text="  只跑前 N 个种").pack(side="left")
        tk.Checkbutton(f, text="不限制（全属，很慢）", variable=self.nolimit_var).pack(side="left", padx=12)

        self._row(top, 2, "框架论文", self.fw_var, browse=lambda: self._pick(self.fw_var, "选框架论文（PDF 抽出的 txt）"))
        self._row(top, 3, "物种清单", self.sp_var, browse=lambda: self._pick(self.sp_var, "选物种清单（每行一个种名）"))
        self._row(top, 4, "CBS 目录", self.cbs_var, browse=lambda: self._pick(self.cbs_var, "选 cbs_catalogue.tsv"))
        self._row(top, 5, "输出目录", self.out_var, browse=lambda: self._pick(self.out_var, "选输出目录", dir=True))
        self._row(top, 6, "NCBI 邮箱", self.email_var, tip="NCBI 要求，随便填个格式正确的即可")

        btns = tk.Frame(self); btns.pack(fill="x", padx=8)
        self.btn_run = tk.Button(btns, text="▶ 开始", width=12, command=self.start)
        self.btn_run.pack(side="left")
        self.btn_stop = tk.Button(btns, text="■ 停止", width=12, command=self.stop, state="disabled")
        self.btn_stop.pack(side="left", padx=6)
        self.btn_open = tk.Button(btns, text="打开结果表", command=self.open_result, state="disabled")
        self.btn_open.pack(side="left", padx=6)
        tk.Label(btns, text="（停止 = 当前种跑完后停，不是硬杀）", fg="#777").pack(side="left", padx=10)

        mid = ttk.LabelFrame(self, text="日志")
        mid.pack(fill="both", expand=True, padx=8, pady=6)
        self.txt = tk.Text(mid, wrap="none", bg="#111", fg="#ddd", insertbackground="#ddd")
        self.txt.pack(side="left", fill="both", expand=True)
        sb = tk.Scrollbar(mid, command=self.txt.yview); sb.pack(side="right", fill="y")
        self.txt.config(yscrollcommand=sb.set)

    def _pick(self, var, title, dir=False):
        p = filedialog.askdirectory(title=title) if dir else filedialog.askopenfilename(title=title)
        if p:
            var.set(p)

    # ---------------- 日志 ----------------
    def log(self, s):
        self.q.put(s if s.endswith("\n") else s + "\n")

    def _drain(self):
        try:
            while True:
                s = self.q.get_nowait()
                self.txt.insert("end", s)
                self.txt.see("end")
        except queue.Empty:
            pass
        self.after(120, self._drain)

    # ---------------- 控制 ----------------
    def start(self):
        if self.running:
            return
        STOP.clear()
        self.running = True
        self.sync_buttons()
        self.txt.delete("1.0", "end")
        # 把 CBS 目录写进一份临时 config（省得用户手改 config/sources.json）
        cfg = None
        cbs = self.cbs_var.get().strip()
        if cbs:
            import json, tempfile
            cfg = os.path.join(tempfile.gettempdir(), "typeaudit_gui_sources.json")
            json.dump({"sources": {
                "cbs": {"adapter": "CBSCatalogue", "path": cbs},
                "genbank": {"adapter": "GenBankAPI", "path": os.path.join(self.out_var.get(), ".cache")},
                "taxonomy": {"adapter": "TaxonomyAPI", "path": os.path.join(self.out_var.get(), ".cache")},
            }}, open(cfg, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
        self.cfg_var = tk.StringVar(value=cfg or "")
        threading.Thread(target=launch, args=(self,), daemon=True).start()

    def stop(self):
        STOP.set()
        self.log("\n⛔ 已请求停止（当前种跑完即停）…")

    def sync_buttons(self):
        self.btn_run.config(state="disabled" if self.running else "normal")
        self.btn_stop.config(state="normal" if self.running else "disabled")

    def done(self):
        g = self.genus_var.get().strip().lower()
        p = os.path.join(self.out_var.get().strip(), "{}_三列表.tsv".format(g))
        if os.path.exists(p):
            self.last_out = p
            self.btn_open.config(state="normal")
            self.log("\n📄 结果表：{}".format(p))

    def open_result(self):
        p = self.last_out
        if not p or not os.path.exists(p):
            messagebox.showinfo("提示", "还没有生成结果表")
            return
        if sys.platform.startswith("win"):
            os.startfile(p)            # noqa: S606
        elif sys.platform == "darwin":
            subprocess.Popen(["open", p])
        else:
            subprocess.Popen(["xdg-open", p])


if __name__ == "__main__":
    App().mainloop()
