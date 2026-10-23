# -*- coding: utf-8 -*-
"""在 Origin 里直接生成图 2（机构构成堆叠柱）与图 3（不一致率 + 95% CI）。

用法：
    "C:\\Users\\Administrator\\AppData\\Local\\Programs\\Python\\Python312\\python.exe" tools\\make_origin_figs.py

前置条件（本机已满足）：
    HKCU\\Software\\Classes\\CLSID\\{B0F21977-...} 下必须有
      LocalServer32  = C:\\Program Files\\OriginLab\\Origin2024\\Origin64.exe
      InprocHandler32= ole32.dll
    （HKLM 里那条写的是坏掉的 8.3 短路径，所以必须在 HKCU 补一份）
"""

import os
import sys
import time

OUT_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "交付_三属对照_20261022", "07_图数据", "origin")

# ── 数据 ────────────────────────────────────────────────────
FIG2_HEADS = ["CBS", "ICMP", "CGMCC", "Other"]          # 列标题
FIG2_MAT = [                                             # 三行 = 三属
    [1085, 286, 182, 2080],
    [850,   21,  99, 3230],
    [542,   36,  91, 2601],
]
FIG3_MAT = [                                             # 属、率、CI 下限、CI 上限
    [1, 12.3, 10.8, 14.0],
    [2, 14.6, 12.6, 17.0],
    [3, 10.3,  8.2, 12.8],
]
FIG3_TOTAL = 12.6


def connect():
    """连 Origin；Dispatch 成功即为连上（本机无 Version 属性，不能用它判断）。"""
    import win32com.client as w
    app = w.Dispatch("Origin.Application")
    app.Visible = 1
    if not app.Execute('type -a "probe";'):
        raise RuntimeError("Origin 连上了但 Execute 调用失败")
    return app


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    app = connect()
    print("[OK] 已连接 Origin")

    # ══════ 图 2 ══════
    print("\n[图2] 建表并写数据 …")
    app.Execute('newbook name:="Fig2" sheet:=1 option:=lsname;')
    ok = app.PutWorksheet("Fig2", FIG2_MAT, 0, 0)        # 从第 0 行第 0 列开始写
    print("[图2] PutWorksheet:", ok)
    for j, h in enumerate(FIG2_HEADS):                   # 写列标题（长名）
        app.Execute('col({})[L]$="{}";'.format(j + 1, h))
    print("[图2] 画堆叠柱（213）…")
    print("[图2] plotxy:", app.Execute('win -a Fig2; plotxy iy:=(1,2:4) plot:=213;'))
    app.Execute('layer.x.from=0.5; layer.x.to=3.5; layer.x.inc=1;')
    f2 = os.path.join(OUT_DIR, "图2_机构构成.png")
    app.Execute('expGraph type:=png path:="{}" filename:="图2_机构构成" '
                'overwrite:=replace tr1.unit:=2 tr1.width:=2400;'.format(OUT_DIR.replace("\\", "/")))
    print("[图2] 导出:", f2, "存在:", os.path.exists(f2))

    # ══════ 图 3 ══════
    print("\n[图3] 建表并写数据 …")
    app.Execute('newbook name:="Fig3" sheet:=1 option:=lsname;')
    print("[图3] PutWorksheet:", app.PutWorksheet("Fig3", FIG3_MAT, 0, 0))
    for j, h in enumerate(["Index", "Rate", "CI_low", "CI_high"]):
        app.Execute('col({})[L]$="{}";'.format(j + 1, h))
    print("[图3] 画散点（201）…")
    print("[图3] plotxy:", app.Execute('win -a Fig3; plotxy iy:=(1,2) plot:=201;'))
    app.Execute('layer.x.from=0.5; layer.x.to=3.5; layer.x.inc=1;')
    app.Execute('layer.y.from=8; layer.y.to=18; layer.y.inc=2;')
    print("[图3] 画 CI 竖线与合计参考线 …")
    for idx, rate, lo, hi in FIG3_MAT:
        app.Execute('draw -l {' + str(idx) + ',' + str(lo) + ',' + str(idx) + ',' + str(hi) + '};')
    app.Execute('draw -l {0.5,' + str(FIG3_TOTAL) + ',3.5,' + str(FIG3_TOTAL) + '};')
    f3 = os.path.join(OUT_DIR, "图3_率与CI.png")
    app.Execute('expGraph type:=png path:="{}" filename:="图3_率与CI" '
                'overwrite:=replace tr1.unit:=2 tr1.width:=2400;'.format(OUT_DIR.replace("\\", "/")))
    print("[图3] 导出:", f3, "存在:", os.path.exists(f3))

    # ══════ 保存工程 ══════
    proj = os.path.join(OUT_DIR, "三属图件.opju").replace("\\", "/")
    app.Execute('save -i "{}";'.format(proj))
    print("\n[保存] OPJU:", proj, "存在:", os.path.exists(proj.replace("/", "\\")))

    print("\n=== 产物 ===")
    for f in sorted(os.listdir(OUT_DIR)):
        print("  {:<26} {:>9,} 字节".format(f, os.path.getsize(os.path.join(OUT_DIR, f))))


if __name__ == "__main__":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    main()
