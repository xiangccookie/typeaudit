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
OUT = os.path.join(NEW, "交付_三属对照_20261022", "07_图数据", "figures_en")
os.makedirs(OUT, exist_ok=True)

BLUE, ORANGE, GREEN, GREY = "#0072B2", "#E69F00", "#009E73", "#999999"

# ── 数据 ───────────────────────────────────────────────
TOTAL, VERIF_T, UNVER_T = 11116, 3040, 8076
SAME, DIFF, REAL, FALSEP = 2481, 559, 403, 156
GEN = ["Colletotrichum", "Trichoderma", "Diaporthe"]
VERIF = [1487, 946, 607]
UNVER = [2152, 3259, 2665]
# 机构构成：来自 tools/rate_audit.py（证据链口径，与论文分母 3,040 同源）
CBS, ICMP, CGMCC = [1032, 837, 494], [280, 21, 36], [175, 88, 77]
BEFORE = [(20.1, 18.1, 22.2), (18.3, 16.0, 20.9), (14.3, 11.8, 17.3)]
AFTER = [(12.6, 11.0, 14.4), (14.9, 12.8, 17.3), (12.2, 9.8, 15.0)]


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
ax.set_xlim(-0.35, 12.4); ax.set_ylim(1.65, 9.45); ax.axis("off")
BOX_W, BOX_H, XC = 6.6, 1.24, 3.70
YS = [8.55, 6.55, 4.55, 2.55]


def box(y, title, n, sub, color):
    ax.add_patch(FancyBboxPatch((XC - BOX_W / 2, y - BOX_H / 2), BOX_W, BOX_H,
                 boxstyle="round,pad=0.04,rounding_size=0.08",
                 linewidth=1.0, edgecolor=color, facecolor="white"))
    ax.text(XC - BOX_W / 2 + 0.22, y + 0.23, title, fontsize=6.4, va="center", ha="left")
    ax.text(XC + BOX_W / 2 - 0.16, y + 0.23, "n = {:,}".format(n), fontsize=7.2,
            va="center", ha="right", fontweight="bold", color=color)
    ax.text(XC - BOX_W / 2 + 0.22, y - 0.30, sub, fontsize=5.8, va="center",
            ha="left", color="#555555", linespacing=1.25)


def arrow(y1, y2):
    ax.add_patch(FancyArrowPatch((XC, y1), (XC, y2), arrowstyle="-|>",
                 mutation_scale=9, linewidth=0.9, color="#444444"))


def excl(y, text, num):
    x_from = XC + 0.22                       # 从流程竖箭头右侧引出
    x_to = XC + BOX_W / 2 + 0.55             # 越过框右缘后停
    ax.add_patch(FancyArrowPatch((x_from, y), (x_to, y), arrowstyle="-|>",
                 mutation_scale=8, linewidth=0.8, color="#888888", linestyle=(0, (2, 1.6))))
    ax.text(x_to + 0.18, y + 0.17, text, fontsize=5.8, va="center", ha="left", color="#444444")
    ax.text(x_to + 0.18, y - 0.26, num, fontsize=5.8, va="center", ha="left", color="#888888")


box(YS[0], "All number-level cells (three genera)", TOTAL,
    "Colletotrichum 3,639 ｜ Trichoderma 4,205 ｜ Diaporthe 3,272", BLUE)
arrow(YS[0] - BOX_H / 2 - 0.02, YS[1] + BOX_H / 2 + 0.02)
excl((YS[0] + YS[1]) / 2, "no machine-readable catalogue, or private working numbers", "n = 8,076（72.7%）")
box(YS[1], "Verifiable numbers (CBS / ICMP / CGMCC)", VERIF_T,
    "27.4% of all numbers; only these enter the denominator\n"
    "Colletotrichum 1,487 ｜ Trichoderma 946 ｜ Diaporthe 607", BLUE)
arrow(YS[1] - BOX_H / 2 - 0.02, YS[2] + BOX_H / 2 + 0.02)
excl((YS[1] + YS[2]) / 2, "Concordant records (no conflict)", "n = 2,481（81.6%）")
box(YS[2], "Number-level discrepancies (before)", DIFF,
    "18.4% of verifiable numbers\n"
    "Colletotrichum 299 ｜ Trichoderma 173 ｜ Diaporthe 87", BLUE)
arrow(YS[2] - BOX_H / 2 - 0.02, YS[3] + BOX_H / 2 + 0.02)
excl((YS[2] + YS[3]) / 2, "reverse audit: judged as not held (false positives)", "n = 156（27.9%）")
box(YS[3], "Number-level discrepancies (after)", REAL,
    "13.3% of verifiable numbers; conservative (95% CI 12.1–14.5%)\n"
    "Colletotrichum 188 ｜ Trichoderma 141 ｜ Diaporthe 74", ORANGE)
print("[图 1]")
save(fig, "Fig1_audit_flow")

# ══════════════════ 图 2：两面板 ══════════════════
fig, axes = plt.subplots(1, 2, figsize=(6.8, 3.05), dpi=600,
                         gridspec_kw={"width_ratios": [1, 1]})
x = np.arange(3)
ax = axes[0]
ax.bar(x, VERIF, 0.55, color=BLUE, label="Verifiable")
ax.bar(x, UNVER, 0.55, bottom=VERIF, color=GREY, label="Non-verifiable")
for i, (v, u) in enumerate(zip(VERIF, UNVER)):
    ax.text(i, v / 2, "{:,}".format(v), ha="center", va="center", color="white", fontsize=6.5)
    ax.text(i, v + u / 2, "{:,}".format(u), ha="center", va="center", color="#333333", fontsize=6.5)
    ax.text(i, v + u + 110, "Verifiable {:.1f}%".format(100 * v / (v + u)),
            ha="center", fontsize=6.3, color=BLUE)
ax.set_xticks(x); ax.set_xticklabels(GEN, fontsize=6.5)
ax.set_ylabel("Number-level cells", fontsize=7); ax.tick_params(labelsize=6.5)
ax.set_ylim(0, 5600)
ax.legend(fontsize=6.5, frameon=False, loc="upper right", ncol=2,
          handlelength=1.2, columnspacing=1.0)
ax.set_title("A  Verifiable vs non-verifiable (counts)", fontsize=7, pad=4)
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
ax.set_ylabel("Share of verifiable numbers per genus (%)", fontsize=7); ax.tick_params(labelsize=6.5)
ax.set_ylim(0, 124)
ax.legend(fontsize=6.5, frameon=False, loc="upper center", ncol=3,
          handlelength=1.2, columnspacing=1.0)
ax.set_title("B  Composition by holding collection", fontsize=7, pad=4)
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
ax.axhline(18.4, color=GREY, ls=(0, (4, 2)), lw=0.8)
ax.axhline(13.3, color=BLUE, ls=(0, (4, 2)), lw=0.8)
ax.text(-0.45, 17.85, "Pooled 18.4%", fontsize=6, color=GREY, ha="left")
ax.text(-0.45, 12.95, "Pooled 13.3%", fontsize=6, color=BLUE, ha="left")
ax.set_xticks(x); ax.set_xticklabels(GEN, fontsize=6.8)
ax.set_ylabel("Number-level discrepancy rate (%)", fontsize=7.2)
ax.set_ylim(7, 23); ax.set_xlim(-0.5, 2.5); ax.tick_params(labelsize=6.5)
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
lg = [Line2D([], [], color=GREY, marker="o", mfc="white", mec=GREY, mew=1.0, ms=4,
             lw=1.0, label="Before reverse audit"),
      Line2D([], [], color=BLUE, marker="o", mfc=BLUE, mec=BLUE, ms=4, lw=1.0,
             label="After (conservative)")]
ax.legend(handles=lg, fontsize=6.5, frameon=False, loc="lower right")
fig.tight_layout()
print("[图 3]")
save(fig, "Fig3_rate_CI")


# ══════════════════ 图 · 发现的四类构成（100% 堆叠条）══════════════════
# 独立成图，不占图 1–3 的编号；文件名不带编号，插进论文时再定图号。
FIND = {
    "gen": ["Colletotrichum", "Trichoderma", "Diaporthe", "Total"],
    "Inconsistent":   [604, 604, 297, 1505],
    "Uncertain":   [1110, 1626, 1131, 3867],
    "Not an error":     [85, 106, 46, 237],
    "Species-level conflict": [33, 7, 120, 160],
}
FCOL = {"Inconsistent": BLUE, "Uncertain": GREY, "Not an error": "#88CCEE", "Species-level conflict": ORANGE}
FORDER = ["Uncertain", "Inconsistent", "Not an error", "Species-level conflict"]

fig, ax = plt.subplots(figsize=(6.8, 3.45), dpi=600)
ys = np.arange(len(FIND["gen"]))[::-1]
tot = np.array([sum(FIND[k][n] for k in FORDER) for n in range(len(ys))], float)
segs = []                       # (row, cat, x0, width%)
for k in FORDER:
    left = np.zeros(len(ys))
    for kk in FORDER:
        if kk == k:
            break
        left += 100.0 * np.array(FIND[kk], float) / tot
    vals = np.array(FIND[k], float)
    w = 100.0 * vals / tot
    ax.barh(ys, w, left=left, height=0.46, color=FCOL[k],
            edgecolor="white", linewidth=0.6, label=k)
    for n in range(len(ys)):
        segs.append((n, k, left[n], w[n]))
    left += w

# 标签：段宽 ≥8% 写在段内；否则引到段外，同一行内按顺序逐级抬高，避免相邻重叠
for n, y in enumerate(ys):
    row = [t for t in segs if t[0] == n]
    row.sort(key=lambda t: t[2])          # 按横坐标从左到右
    lvl = 0
    for _, k, x0, w in row:
        if w >= 8:
            ax.text(x0 + w / 2, y, "{:.1f}%".format(w), ha="center", va="center",
                    fontsize=6.2, color=("white" if k in ("Inconsistent", "Uncertain") else "black"))
        else:
            dy = 0.26 + 0.26 * lvl
            lvl += 1
            ax.annotate("{:.1f}%".format(w), xy=(x0 + w / 2, y + 0.24),
                        xytext=(x0 + w / 2, y + dy), ha="center", va="bottom",
                        fontsize=6.0, color="#222222",
                        arrowprops=dict(arrowstyle="-", lw=0.5, color="#777777",
                                        shrinkA=1, shrinkB=0))

for n, y in enumerate(ys):
    ax.text(103, y, "{:,}".format(int(tot[n])), va="center", ha="left",
            fontsize=6.4, color="#333333")
ax.set_yticks(ys); ax.set_yticklabels(FIND["gen"], fontsize=6.8)
ax.set_xlim(0, 112); ax.set_xticks([0, 25, 50, 75, 100])
ax.set_xticklabels(["0", "25", "50", "75", "100%"], fontsize=6.2)
ax.set_xlabel("Share of all findings in the row", fontsize=7.0)
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
ax.tick_params(axis="y", length=0)
ax.legend(ncol=4, fontsize=6.2, frameon=False, loc="upper center",
          bbox_to_anchor=(0.42, -0.32), handlelength=1.1, columnspacing=1.2)
fig.tight_layout()
print("[图 · 四类构成]")
save(fig, "Fig_findings_composition")


# ══════════════════ 图 · 号级不一致的内部构成（分组条形）══════════════════
# 对应原「表 1」。用分组条形而不是堆叠条：与图 3（100% 堆叠）区分开，
# 且能同时看出「Trichoderma 的疑似非机构编号 396 是唯一的大块」这一事实。
SUBCAT = ["Not an institutional number", "Declared type, no record in collection", "Number absent", "Other minor classes"]
SUBD = {
    "Colletotrichum": [134, 86, 131, 253],
    "Trichoderma":    [396, 73, 43, 92],
    "Diaporthe":      [138, 35, 14, 110],
}
SUBCOL = {"Colletotrichum": BLUE, "Trichoderma": "#56B4E9", "Diaporthe": ORANGE}

fig, ax = plt.subplots(figsize=(6.8, 2.9), dpi=600)
yy = np.arange(len(SUBCAT))[::-1]
h = 0.24
for i, (g, vals) in enumerate(SUBD.items()):
    pos = yy + (1 - i) * h                     # 三个属上下错开
    ax.barh(pos, vals, height=h * 0.88, color=SUBCOL[g], label=g,
            edgecolor="white", linewidth=0.4)
    for p, v in zip(pos, vals):
        ax.text(v + 6, p, str(v), va="center", ha="left", fontsize=5.8, color="#333333")

tot = [sum(SUBD[g][k] for g in SUBD) for k in range(len(SUBCAT))]
for k, y in enumerate(yy):
    ax.text(556, y, "{:,}".format(tot[k]), va="center", ha="right",
            fontsize=6.4, color="#111111", fontweight="bold")

ax.set_yticks(yy); ax.set_yticklabels(SUBCAT, fontsize=6.6)
ax.set_xlim(0, 560); ax.set_ylim(-0.62, len(SUBCAT) - 0.38)
ax.set_xticks([0, 100, 200, 300, 400]); ax.tick_params(labelsize=6.2)
ax.set_xlabel("Number of findings (bold: three-collection total)", fontsize=6.9)
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
ax.tick_params(axis="y", length=0)
ax.legend(ncol=3, fontsize=6.2, frameon=False, loc="upper center",
          handlelength=1.0, columnspacing=1.2, bbox_to_anchor=(0.45, -0.26))
fig.tight_layout()
print("[图 · 不一致内部构成]")
save(fig, "Fig_discrepancy_breakdown")


# ══════════════════ 图 · 判定流程（两个层级：种 / 记录）══════════════════
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
fig, ax = plt.subplots(figsize=(6.6, 4.6), dpi=600)
ax.set_xlim(0, 9.4); ax.set_ylim(1.55, 8.85); ax.axis("off")


def nd(x, y, w, h, txt, ec, fc="white", fs=6.4, bold=False):
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h,
                 boxstyle="round,pad=0.03,rounding_size=0.07",
                 linewidth=0.95, edgecolor=ec, facecolor=fc))
    ax.text(x, y, txt, ha="center", va="center", fontsize=fs,
            fontweight=("bold" if bold else "normal"), linespacing=1.3)


def ar(x1, y1, x2, y2, dash=False):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>",
                 mutation_scale=8, linewidth=0.85, color="#555555",
                 linestyle=((0, (2.4, 1.8)) if dash else "solid")))


# ── 层级一：一个种（名字层，不经三步）──
nd(2.15, 8.35, 2.30, 0.55, "One species", "#B85C00", fc="#FFF1E0", fs=6.9, bold=True)
ar(3.32, 8.35, 4.72, 8.35)
ax.text(4.02, 8.47, "databases disagree", fontsize=5.5, color="#B85C00", ha="center")
nd(6.75, 8.35, 3.90, 0.60, "Species-level conflict  160\nNCBI keeps the species, MycoBank merges it", "#B85C00", fc="#FFF1E0", fs=5.9)

# ── 层级二：一条记录（号层，三步）──
XM, XR = 2.75, 6.95
nd(XM, 7.20, 3.10, 0.55, "One record (species × number)", "#333333", fs=6.9, bold=True)
ar(XM, 6.92, XM, 6.52)
ax.text(XM + 0.09, 6.74, "yes", fontsize=5.6, color="#666666", ha="left")

nd(XM, 6.12, 3.50, 0.74, "① Is the number checkable?\npublic catalogue available, number retrievable", BLUE)
ar(XM, 5.75, XM, 5.35)
ax.text(XM + 0.09, 5.57, "yes", fontsize=5.6, color="#666666", ha="left")

nd(XM, 4.94, 3.50, 0.74, "② One of six known relations?\nnew combination / gender ending / spelling / infraspecific / reassigned / unchanged", BLUE)
ar(XM, 4.57, XM, 4.17)
ax.text(XM + 0.09, 4.39, "yes", fontsize=5.6, color="#666666", ha="left")

nd(XM, 3.76, 3.50, 0.74, "③ Do both sides assert positively?\n(not just one side 'not found')", BLUE)
ar(XM, 3.39, XM, 3.00)
ax.text(XM + 0.09, 3.21, "yes", fontsize=5.6, color="#666666", ha="left")

nd(XM, 2.68, 2.55, 0.56, "Inconsistent  1,505", ORANGE, fs=7.0, bold=True)

for y, txt in ((6.12, "Uncertain  3,867"), (4.94, "Not an error  237"), (3.76, "Uncertain (as above)")):
    ar(XM + 1.78, y, XR - 1.25, y, dash=True)
    ax.text((XM + 1.78 + XR - 1.25) / 2, y + 0.10, "no", fontsize=5.6, color="#666666", ha="center")
    nd(XR, y, 2.45, 0.54, txt, GREY, fs=6.6)

ax.text(XM - 1.28, 2.05,
        "Three classes need the protologue: number absent | name mismatch | type status mismatch\n→ 30-case validation sample (10 per class)",
        fontsize=5.9, color="#7a4a00", ha="left", va="top", linespacing=1.35)
print("[图 · 判定流程]")
save(fig, "Fig_decision_flow")

print("\n全部完成，目录：", OUT)
