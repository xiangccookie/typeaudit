# -*- coding: utf-8 -*-
"""临时脚本：按 rate_audit 的公开口径，重算三属的可核号与错值（新数据）。"""
import io
import os
import re
import sys
import importlib.util as ilu

R = r"C:\Users\Administrator\Desktop\sol\奇怪\colleto-type-audit"
spec = ilu.spec_from_file_location("_pw", os.path.join(R, "tools", "pipeline_web.py"))
pw = ilu.module_from_spec(spec)
spec.loader.exec_module(pw)
cls = pw.classify_flag

VER = {"CBS", "ICMP", "CGMCC"}
PRE = re.compile(r"^\((?:strain|isolate|culture|specimen|other)\)\s*", re.I)
TOK = re.compile(r"^([A-Za-z][\w\-\.]*)\s*[:#]?\s*([\d\.]+[A-Za-z]?)\s*(?:\[(.*?)\])?\s*$")

PAIRS = (("out_colleto_web", "colletotrichum"),
         ("out_tricho_web", "trichoderma"),
         ("out_dia_web", "diaporthe"))

tot_o = tot_e = 0
print("属                 可核号   错值")
for out, slug in PAIRS:
    p = os.path.join(R, out, slug + "_三列表.tsv")
    occ = err = 0
    for row in io.open(p, encoding="utf-8-sig").read().splitlines()[1:]:
        c = row.split("\t")
        if len(c) < 3:
            continue
        for cell in c[1:3]:
            for seg in cell.split("｜"):
                s = PRE.sub("", seg.strip())
                if not s or s == "—":
                    continue
                m = TOK.match(s)
                if not m or m.group(1).upper() not in VER:
                    continue
                occ += 1
                fl = (m.group(3) or "").strip()
                if fl and fl != "—" and cls(fl) == "错值":
                    err += 1
    print("{:<18} {:>6} {:>6}".format(out, occ, err))
    tot_o += occ
    tot_e += err
print("{:<18} {:>6} {:>6}".format("合计", tot_o, tot_e))
print()
print("可核率      : {}/11,103 = {:.1f}%".format(tot_o, 100.0 * tot_o / 11103))
print("未修正率    : {}/{} = {:.1f}%".format(tot_e, tot_o, 100.0 * tot_e / tot_o))
mv = tot_e - 156
print("移出 156 后 : {}/{} = {:.1f}%".format(mv, tot_o, 100.0 * mv / tot_o))
print("假阳性占比  : 156/{} = {:.1f}%".format(tot_e, 100.0 * 156 / tot_e))
print("高估倍数    : {:.2f}".format((100.0 * tot_e / tot_o) / (100.0 * mv / tot_o)))
