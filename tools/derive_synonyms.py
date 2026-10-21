# -*- coding: utf-8 -*-
"""
derive_synonyms.py —— 从 MycoBank 全表【反推】某个属的同义属清单。

【为什么要反推，而不是手写】
    手写同义属 = 把记忆当证据。本工具的作者就干过一次：手写
    "Glomerella,Gloeosporium,Vermicularia" 三个，反推后发现 ≥2 次的有 20 个，
    漏了 17 个（其中还有名字里就带 Colletotrichum 的 Colletotrichopsis）。
    所以清单必须可复现、可溯源，换任何属都能自动得到。

【规则（一行）】
    凡「当前名」首词是属 G 的记录，它们的「Taxon name」首词里出现的其他属名，
    就是 G 的同义属候选。出现次数越多，越是真正的同义关系。

【输出是候选，不是权威 —— 必须过筛】
    反推会带进噪声，两种典型：
      ① 形态差太远的属（如 Penicillium 下出现 Mucor）
      ② 双向互挂（Aspergillus 挂在 Penicillium 下，Penicillium 又挂在 Aspergillus 下）
         —— 逻辑上不可能同时成立，本工具会给这类打 mutual 标记
    所以输出 tsv 里带 n 和 mutual 两列，供你按需筛。
    但注意 mutual 只是「提示去看」，不是丢弃依据 —— 实测 Gloeosporium(×37)、
    Vermicularia(×37) 这两个最重要的真同义属也会被标 mutual（确有物种被移动过）。

用法
----
    python tools/derive_synonyms.py Colletotrichum
    python tools/derive_synonyms.py Colletotrichum Fusarium Penicillium Aspergillus
    python tools/derive_synonyms.py Colletotrichum --min-count 1        # 连只出现 1 次的也收
    python tools/derive_synonyms.py --list                              # 看已经派生过哪些属

产物
----
    data/synonyms_<属>.tsv        genus / synonym_genus / n / n_synonym / statuses / mutual
    data/synonyms_<属>.meta.json  溯源（源文件/大小/修改时间/扫描行数/规则/参数/时间）
"""
import sys
try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

import argparse, io, json, os, re, time
from collections import Counter, defaultdict

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# 本机已知的 MycoBank 全表位置（按优先级）
MYCOBANK_CANDIDATES = [
    os.path.join(REPO, "data", "mycobank_full.tsv"),
    r"C:\Users\Administrator\Desktop\Tool\脚本\调取mycobank数据，分析sp_curret,编码，是否为复合体\mycobank_full.tsv",
]


def resolve_mycobank(explicit=None):
    for p in ([explicit] if explicit else []) + MYCOBANK_CANDIDATES:
        if p and os.path.exists(p):
            return p
    return None


def first_word(name):
    """取名字首词并清洗；不像样的属名返回空（过滤作者缩写、空值、杂字符）。"""
    if not name:
        return ""
    s = name.strip()
    if not s:
        return ""
    tok = re.sub(r"[^A-Za-z]", "", s.split()[0])
    if len(tok) < 4 or not tok[0].isupper():
        return ""
    return tok


def scan(path, verbose=True):
    """
    扫描一遍 MycoBank 全表，建立两个倒排表：
        inv[当前名属][Taxon name 属] = 出现次数
        st [当前名属][Taxon name 属][Name status] = 出现次数
    一次扫描即可服务任意多个目标属，还能免费拿到「双向互挂」标记。
    """
    inv = defaultdict(Counter)
    st = defaultdict(lambda: defaultdict(Counter))
    n = 0
    with io.open(path, encoding="utf-8", errors="replace") as f:
        head = f.readline().rstrip("\n").rstrip("\r").split("\t")
        ci = {h: i for i, h in enumerate(head)}
        i_cur = ci.get("Current name (bare)", ci.get("Current name"))
        i_tax = ci.get("Taxon name")
        i_st = ci.get("Name status")
        if i_tax is None or i_cur is None:
            raise SystemExit("MycoBank 表缺少 'Taxon name' 或 'Current name' 列，现有：{}".format(head))
        for line in f:
            if not line.strip():
                continue
            p = line.rstrip("\n").rstrip("\r").split("\t")
            if len(p) <= max(i_cur, i_tax):
                continue
            n += 1
            cur = first_word(p[i_cur])
            tax = first_word(p[i_tax])
            if not cur or not tax or cur.lower() == tax.lower():
                continue
            inv[cur][tax] += 1
            if i_st is not None and len(p) > i_st:
                st[cur][tax][p[i_st].strip()] += 1
    if verbose:
        print("  扫描完成：{:,} 行，反推出 {:,} 个「当前名属」，共 {:,} 组属对".format(
            n, len(inv), sum(len(v) for v in inv.values())))
    return inv, st, n, head


def derive(genus, inv, st, min_count=2, status_filter=None):
    """取出目标属的同义属候选。"""
    cand = inv.get(genus, Counter())
    rows = []
    for tax, n in cand.most_common():
        if n < min_count:
            continue
        statuses = st.get(genus, {}).get(tax, Counter())
        if status_filter and statuses.get(status_filter, 0) == 0:
            continue
        # 双向互挂：对方也把本属列为同义属 → 逻辑上不可能同时成立，标出来
        mutual = inv.get(tax, Counter()).get(genus, 0) > 0
        rows.append({
            "genus": genus, "synonym_genus": tax, "n": n,
            "n_synonym": statuses.get("Synonym", 0),
            "statuses": ";".join("{}:{}".format(k, v) for k, v in statuses.most_common()),
            "mutual": 1 if mutual else 0,
        })
    return rows


def main():
    ap = argparse.ArgumentParser(description="从 MycoBank 反推某属的同义属候选")
    ap.add_argument("genera", nargs="*", help="一个或多个属名（首字母大写）")
    ap.add_argument("--mycobank", default=None, help="MycoBank 全表路径")
    ap.add_argument("--min-count", type=int, default=2,
                    help="最少出现几次才算候选（默认 2；1 会把偶然/笔误也收进来）")
    ap.add_argument("--status", default=None,
                    help="只保留含该 Name status 的候选。【实测基本没用】："
                         "Name status 记的是名字本身的状态（Legitimate/Invalid/…），"
                         "不是同义关系标记，Synonym 计数恒为 0")
    ap.add_argument("--out-dir", default=os.path.join(REPO, "data"))
    ap.add_argument("--list", action="store_true", help="列出已派生过的属，然后退出")
    a = ap.parse_args()

    if a.list:
        d = a.out_dir
        got = sorted(f for f in os.listdir(d) if f.startswith("synonyms_") and f.endswith(".tsv")) \
            if os.path.isdir(d) else []
        print("已派生的同义属文件（{}）：".format(d))
        for f in got:
            print("   {}".format(f))
        if not got:
            print("   （还没有）")
        return

    if not a.genera:
        raise SystemExit("请给出至少一个属名，例如：python tools/derive_synonyms.py Colletotrichum")

    my = resolve_mycobank(a.mycobank)
    if not my:
        raise SystemExit("没找到 MycoBank 全表。用 --mycobank 指定，或跑 tools/fetch_refdata.py 探测。")
    size = os.path.getsize(my)
    print("=" * 96)
    print("从 MycoBank 反推同义属")
    print("=" * 96)
    print("  源表：{}  （{:,} 字节，修改时间 {}）".format(
        my, size, time.strftime("%Y-%m-%d %H:%M", time.localtime(os.path.getmtime(my)))))
    inv, st, nrows, head = scan(my)

    if not os.path.isdir(a.out_dir):
        os.makedirs(a.out_dir)

    for genus in a.genera:
        rows = derive(genus, inv, st, min_count=a.min_count, status_filter=a.status)
        tsv = os.path.join(a.out_dir, "synonyms_{}.tsv".format(genus))
        with io.open(tsv, "w", encoding="utf-8", newline="") as f:
            f.write("genus\tsynonym_genus\tn\tn_synonym\tstatuses\tmutual\n")
            for r in rows:
                f.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(
                    r["genus"], r["synonym_genus"], r["n"], r["n_synonym"],
                    r["statuses"], r["mutual"]))
        meta = {
            "genus": genus,
            "derived_at": time.strftime("%Y-%m-%d %H:%M:%S"),
            "source_file": my, "source_size": size, "source_mtime": os.path.getmtime(my),
            "source_columns": head, "rows_scanned": nrows,
            "rule": "凡 Current name 首词为属 G 的记录，其 Taxon name 首词（≠G）即为 G 的同义属候选",
            "min_count": a.min_count, "status_filter": a.status,
            "n_candidates": len(rows), "output": tsv,
        }
        with io.open(os.path.join(a.out_dir, "synonyms_{}.meta.json".format(genus)),
                     "w", encoding="utf-8") as f:
            json.dump(meta, f, ensure_ascii=False, indent=1)
        print()
        print("=" * 96)
        print("【{}】候选 {} 个（最少出现 {} 次）".format(genus, len(rows), a.min_count))
        print("=" * 96)
        mut = [r for r in rows if r["mutual"]]
        for r in rows:
            print("   {:<26} ×{:<5} Synonym×{:<5} {}".format(
                r["synonym_genus"], r["n"], r["n_synonym"],
                "  <-- 双向互挂，需人工核" if r["mutual"] else ""))
        if mut:
            print()
            print("   ! 双向互挂 {} 个（对方也把 {} 列为同义属）：{}".format(
                len(mut), genus, ", ".join(r["synonym_genus"] for r in mut)))
            print("     这不可能同时成立，属噪声信号，建议逐个人工核过再用。")
        print()
        print("   已写 {} 和同名 .meta.json".format(tsv))


if __name__ == "__main__":
    main()
