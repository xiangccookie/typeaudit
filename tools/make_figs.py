# -*- coding: utf-8 -*-
"""一次性生成论文三张图：PNG / PDF（矢量）/ TIFF（600 dpi）。

用法：
    C:\\Users\\Administrator\\AppData\\Local\\Programs\\Python\\Python312\\python.exe tools\\make_figs.py
（bio-genie 环境也自带 matplotlib，用哪个解释器都行）

产物目录：交付_三属对照_20261022/07_图数据/figures/
"""

import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ── 全局样式：期刊要求字体嵌入 + 字号下限 ──
plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei"]
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["pdf.fonttype"] = 42        # TrueType 嵌入（避免 Type 3 被拒）
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["axes.linewidth"] = 0.75
plt.rcParams["xtick.direction"] = "in"
plt.rcParams["ytick.direction"] = "in"
plt.rcParams["font.size"] = 7

NEW = r"C:\Users\Administrator\Desktop\sol\奇怪\colleto-type-audit"
OUT = os.path.join(NEW, "交付_三属对照_20261022", "07_图数据", "figures")
os.makedirs(OUT, exist_ok=True)

BLUE, ORANGE, GREEN, GREY = "#0072B2", "#E69F00", "#009E73", "#999999"

# ── 数据 ───────────────────────────────────────────────
TOTAL, VERIF_T, UNVER_T = 11103, 3192, 7911
SAME, DIFF, REAL, FALSEP = 2634, 558, 402, 156
GEN = ["Colletotrichum", "Trichoderma", "Diaporthe"]
VERIF = [1553, 970, 669]
UNVER = [2080, 3230, 2601]
CBS, ICMP, CGMCC = [1085, 850, 542], [286, 21, 36], [182, 99, 91]
BEFORE = [(19.4, 17.6, 21.5), (17.9, 15.7, 20.5), (12.3, 10.0, 15.0)]
AFTER = [(12.3, 10.8, 14.0), (14.6, 12.6, 17.0), (10.3, 8.2, 12.8)]


def save(fig, stem):
    """同图存三种格式：PNG（预览）、PDF（矢量投稿）、TIFF（600 dpi 位图）。"""
    paths = []
    p = os.path.join(OUT, stem + ".png"); fig.savefig(p, dpi=600); paths.append(p)
    p = os.path.join(OUT, stem + ".pdf"); fig.savefig(p); paths.append(p)
    p = os.path.join(OUT, stem + ".tif")
    fig.savefig(p, dpi=600, pil_kwargs={"compression": "tiff_lzw"})
    # 期刊只收 RGB/灰度，不收回 alpha 的 RGBA
    from PIL import Image as _Img
    _im = _Img.open(p).convert("RGB")
    _im.save(p, dpi=(600, 600), compression="tiff_lzw")
    paths.append(p)
    plt.close(fig)
    for q in paths:
        print("  {:<34} {:>10,} 字节".format(os.path.basename(q), os.path.getsize(q)))


# ══════════════════ 图 1：流程框图 ══════════════════
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
fig, ax = plt.subplots(figsize=(6.85, 4.5), dpi=600)
ax.set_xlim(-0.35, 10.7); ax.set_ylim(1.65, 9.45); ax.axis("off")
BOX_W, BOX_H, XC = 6.6, 1.02, 3.70
YS = [8.55, 6.55, 4.55, 2.55]


def box(y, title, n, sub, color):
    ax.add_patch(FancyBboxPatch((XC - BOX_W / 2, y - BOX_H / 2), BOX_W, BOX_H,
                 boxstyle="round,pad=0.04,rounding_size=0.08",
                 linewidth=1.0, edgecolor=color, facecolor="white"))
    ax.text(XC - BOX_W / 2 + 0.22, y + 0.23, title, fontsize=7.2, va="center", ha="left")
    ax.text(XC + BOX_W / 2 - 0.22, y + 0.23, "n = {:,}".format(n), fontsize=8.0,
            va="center", ha="right", fontweight="bold", color=color)
    ax.text(XC - BOX_W / 2 + 0.22, y - 0.25, sub, fontsize=6.0, va="center",
            ha="left", color="#555555")


def arrow(y1, y2):
    ax.add_patch(FancyArrowPatch((XC, y1), (XC, y2), arrowstyle="-|>",
                 mutation_scale=9, linewidth=0.9, color="#444444"))


def excl(y, text, num):
    x_from = XC + 0.22                       # 从流程竖箭头右侧引出
    x_to = XC + BOX_W / 2 + 0.55             # 越过框右缘后停
    ax.add_patch(FancyArrowPatch((x_from, y), (x_to, y), arrowstyle="-|>",
                 mutation_scale=8, linewidth=0.8, color="#888888", linestyle=(0, (2, 1.6))))
    ax.text(x_to + 0.18, y + 0.17, text, fontsize=6.3, va="center", ha="left", color="#444444")
    ax.text(x_to + 0.18, y - 0.26, num, fontsize=6.3, va="center", ha="left", color="#888888")


box(YS[0], "三属全部号级单元格", TOTAL,
    "Colletotrichum 3,633 ｜ Trichoderma 4,200 ｜ Diaporthe 3,270", BLUE)
arrow(YS[0] - BOX_H / 2 - 0.02, YS[1] + BOX_H / 2 + 0.02)
excl((YS[0] + YS[1]) / 2, "无可机读目录，或仅为研究者工作编号", "n = 7,911（71.3%）")
box(YS[1], "可核号（CBS / ICMP / CGMCC）", VERIF_T,
    "占全部号 28.7%；仅此部分进入率的分母", BLUE)
arrow(YS[1] - BOX_H / 2 - 0.02, YS[2] + BOX_H / 2 + 0.02)
excl((YS[1] + YS[2]) / 2, "记录一致（GenBank 与机构目录无冲突）", "n = 2,634（82.5%）")
box(YS[2], "号级不一致（修正前）", DIFF, "占可核号 17.5%", BLUE)
arrow(YS[2] - BOX_H / 2 - 0.02, YS[3] + BOX_H / 2 + 0.02)
excl((YS[2] + YS[3]) / 2, "反向审计：判为机构未收录（实为假阳性）", "n = 156（28.0%）")
box(YS[3], "号级不一致（修正后）", REAL,
    "占可核号 12.6%；保守口径（95% CI 11.5–13.8%）", ORANGE)
print("[图 1]")
save(fig, "Fig1_audit_flow")

# ══════════════════ 图 2：两面板 ══════════════════
fig, axes = plt.subplots(1, 2, figsize=(6.8, 3.05), dpi=600,
                         gridspec_kw={"width_ratios": [1, 1]})
x = np.arange(3)
ax = axes[0]
ax.bar(x, VERIF, 0.55, color=BLUE, label="可核")
ax.bar(x, UNVER, 0.55, bottom=VERIF, color=GREY, label="不可核")
for i, (v, u) in enumerate(zip(VERIF, UNVER)):
    ax.text(i, v / 2, "{:,}".format(v), ha="center", va="center", color="white", fontsize=6.5)
    ax.text(i, v + u / 2, "{:,}".format(u), ha="center", va="center", color="#333333", fontsize=6.5)
    ax.text(i, v + u + 110, "可核 {:.1f}%".format(100 * v / (v + u)),
            ha="center", fontsize=6.3, color=BLUE)
ax.set_xticks(x); ax.set_xticklabels(GEN, fontsize=6.5)
ax.set_ylabel("号级单元格数", fontsize=7); ax.tick_params(labelsize=6.5)
ax.set_ylim(0, 5600)
ax.legend(fontsize=6.5, frameon=False, loc="upper right", ncol=2,
          handlelength=1.2, columnspacing=1.0)
ax.set_title("A  可核与不可核（绝对数）", fontsize=7, pad=4)
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

ax = axes[1]
tot = np.array(VERIF, float)
ax.bar(x, np.array(CBS) / tot * 100, 0.55, color=BLUE, label="CBS")
ax.bar(x, np.array(ICMP) / tot * 100, 0.55, bottom=np.array(CBS) / tot * 100,
       color=ORANGE, label="ICMP")
ax.bar(x, np.array(CGMCC) / tot * 100, 0.55,
       bottom=(np.array(CBS) + np.array(ICMP)) / tot * 100, color=GREEN, label="CGMCC")
for i in range(3):
    ax.text(i, CBS[i] / tot[i] * 100 / 2, "{:.0f}%".format(100 * CBS[i] / tot[i]),
            ha="center", va="center", color="white", fontsize=6.5)
    y_ic = (CBS[i] + ICMP[i] / 2) / tot[i] * 100
    txt = "ICMP {:.1f}%".format(100 * ICMP[i] / tot[i])
    if ICMP[i] / tot[i] * 100 > 8:
        ax.text(i, y_ic, txt, ha="center", va="center", color="#4a2f00", fontsize=5.8)
    else:
        ax.text(i, y_ic + 3.2, txt, ha="center", va="bottom", color="#8a5a00", fontsize=5.8)
ax.set_xticks(x); ax.set_xticklabels(GEN, fontsize=6.5)
ax.set_ylabel("占该属可核号的比例（%）", fontsize=7); ax.tick_params(labelsize=6.5)
ax.set_ylim(0, 124)
ax.legend(fontsize=6.5, frameon=False, loc="upper center", ncol=3,
          handlelength=1.2, columnspacing=1.0)
ax.set_title("B  可核号内的机构构成", fontsize=7, pad=4)
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
fig.tight_layout()
print("[图 2]")
save(fig, "Fig2_verifiability_composition")

# ══════════════════ 图 3：率 + CI ══════════════════
from matplotlib.lines import Line2D
fig, ax = plt.subplots(figsize=(4.72, 3.4), dpi=600)
off = 0.13
for i, (m, lo, hi) in enumerate(BEFORE):
    ax.plot([i - off] * 2, [lo, hi], color=GREY, lw=1.0)
    ax.plot(i - off, m, "o", ms=4.2, mfc="white", mec=GREY, mew=1.0, zorder=3)
for i, (m, lo, hi) in enumerate(AFTER):
    ax.plot([i + off] * 2, [lo, hi], color=BLUE, lw=1.0)
    ax.plot(i + off, m, "o", ms=4.2, mfc=BLUE, mec=BLUE, zorder=3)
ax.axhline(17.5, color=GREY, ls=(0, (4, 2)), lw=0.8)
ax.axhline(12.6, color=BLUE, ls=(0, (4, 2)), lw=0.8)
ax.text(-0.45, 17.85, "合计 17.5%", fontsize=6, color=GREY, ha="left")
ax.text(-0.45, 12.95, "合计 12.6%", fontsize=6, color=BLUE, ha="left")
ax.set_xticks(x); ax.set_xticklabels(GEN, fontsize=6.8)
ax.set_ylabel("号级不一致率（%）", fontsize=7.2)
ax.set_ylim(7, 23); ax.set_xlim(-0.5, 2.5); ax.tick_params(labelsize=6.5)
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
lg = [Line2D([], [], color=GREY, marker="o", mfc="white", mec=GREY, mew=1.0, ms=4,
             lw=1.0, label="修正前"),
      Line2D([], [], color=BLUE, marker="o", mfc=BLUE, mec=BLUE, ms=4, lw=1.0,
             label="修正后（保守）")]
ax.legend(handles=lg, fontsize=6.5, frameon=False, loc="lower right")
fig.tight_layout()
print("[图 3]")
save(fig, "Fig3_rate_CI")

print("\n全部完成，目录：", OUT)
