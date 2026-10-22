#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""把所有来源的发现汇总成一张表 —— 「全部发现」。

为什么要它
----------
原来每类发现散在不同产物里：
  · 号层标记        → <out>/<genus>_三列表.tsv
  · 名字层分歧      → <out>/四方核验.tsv
  · 记录层声明缺失  → out_reverse/反向审计.tsv
  · 目录层聚合问题  → out/quality_report.json
用户在网页上只看得到号层那一层，于是产生「怎么没抓到并名/异名」的错觉。
本脚本把它们合成一张表，每条都带【层级 / 归类 / 来源 / 证据】。

归类口径（与网页一致）
--------------------
  错值   = 有外部依据能判对错的
  不确定 = 能指出问题但判不出谁错（分歧、待核、声明缺失）
  非错   = 常见的"看起来像错但其实不是"（异名、分类学滞后）
"""
import argparse, io, json, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(REPO, "src"))

COLS = ["序号", "层级", "类别", "对象", "发现", "归类", "来源", "证据"]


def load_flag_classes():
    """复用网页的归类表，保证同一标记在两边归类一致。"""
    import importlib.util
    p = os.path.join(HERE, "pipeline_web.py")
    spec = importlib.util.spec_from_file_location("_pw", p)
    m = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(m)
        return m.classify_flag
    except Exception:
        def _f(x):
            return "不确定"
        return _f


def from_table(out_dir, genus, classify):
    """① 号层：三列表里每个带标记的单元格 → 每条标记一行。"""
    rows = []
    p = os.path.join(out_dir, "{}_三列表.tsv".format(genus.lower()))
    if not os.path.exists(p):
        return rows
    for line in io.open(p, encoding="utf-8-sig").read().splitlines()[1:]:
        f = line.split("\t")
        if len(f) < 3:
            continue
        sp = f[0].strip()
        for col, cname in ((1, "菌株号"), (2, "模式菌")):
            for cell in f[col].split(" ｜ "):
                cell = cell.strip()
                if not cell or cell == "—":
                    continue
                m = re.match(r"^(.*?)\s*\[(.*)\]\s*$", cell)
                if not m:
                    continue
                num = m.group(1).strip()
                for fl in [x.strip() for x in m.group(2).split(";")
                           if x.strip() and x.strip() not in ("—", "-")]:
                    # 「同名异物：…」「已核：…」是种级裁决（我的名字层给的判决），
                    # 不是号层标记；同一现象已由「异名/独立种分歧」记过一条，
                    # 在这里再按号记一次会让 4 个种各多 1 条（实测 1,596 → 1,600）。
                    if fl.startswith(("同名异物", "已核")):
                        continue
                    rows.append({"层级": "号", "类别": fl.split("（")[0],
                                 "对象": "{} / {}".format(sp, num),
                                 "发现": fl, "归类": classify(fl),
                                 "来源": "三列表·{}列".format(cname),
                                 "证据": "out={} ｜ 列={} ｜ 号={}".format(
                                     os.path.basename(out_dir), cname, num)})
    return rows


def from_species(out_dir, classify):
    """② 名字层：四方核验里的矛盾 / 待核 / 名称分歧。"""
    rows = []
    p = os.path.join(out_dir, "四方核验.tsv")
    if not os.path.exists(p):
        return rows
    lines = io.open(p, encoding="utf-8").read().splitlines()
    hdr = lines[0].split("\t")
    for l in lines[1:]:
        r = dict(zip(hdr, l.split("\t")))
        sp = r.get("种", "")
        for field, tier, cls in (("矛盾明细", "名字", "错值"),
                                 ("待核明细", "名字", "不确定")):
            v = (r.get(field) or "").strip()
            for item in [x for x in v.split(" ｜ ") if x.strip()]:
                # 去重：名字层分歧已由「名称分歧」单独出一条，
                # 这里再记一遍会让同一件事出现两次（实测 4 个分歧各出现两遍）。
                if field == "矛盾明细" and "当独立种" in item:
                    continue
                rows.append({"层级": tier, "类别": field.replace("明细", ""),
                             "对象": sp, "发现": item, "归类": cls,
                             "来源": "四方核验", "证据": "out={} ｜ {}".format(
                                 os.path.basename(out_dir), field)})
        nj = (r.get("名称分歧") or "").strip()
        if nj:
            try:
                d = json.loads(nj)
            except Exception:
                d = {}
            if d.get("verdict") == "拼写变体":
                # 拼写变体：同一个名字的两种写法（MycoBank 标 Orthographic variant）。
                # 它不该进「异名/独立种分歧」，也不该算不一致 —— 归「非错」。
                rows.append({"层级": "名字", "类别": "拼写变体（非错）",
                             "对象": sp,
                             "发现": "MycoBank 标为 Orthographic variant（{}）；NCBI 另立 taxid".format(
                                 d.get("mycobank", "")),
                             "归类": "非错",
                             "来源": "四方核验·名称分歧",
                             "证据": "status={} ｜ taxid={}".format(
                                 d.get("status", ""), r.get("NCBI_taxid", ""))})
            elif d.get("verdict") == "分歧":
                rows.append({"层级": "名字", "类别": "异名/独立种分歧",
                             "对象": sp,
                             "发现": "NCBI：{} ｜ MycoBank：{}".format(
                                 d.get("ncbi", ""), d.get("mycobank", "")),
                             "归类": d.get("level") or "不确定",
                             "来源": "四方核验·名称分歧",
                             "证据": "taxid={} ｜ 接受名={}".format(
                                 r.get("NCBI_taxid", ""), r.get("MycoBank接受名", "")[:60])})
    return rows



def _load_pw():
    """加载 tools/pipeline_web.py 模块（只为拿 name_relation 这一个函数）。"""
    import importlib.util
    p = os.path.join(HERE, "pipeline_web.py")
    spec = importlib.util.spec_from_file_location("_pw_names", p)
    m = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(m)
        return m
    except Exception:
        return None


def from_names(out_dir, genus, repo, classify):
    """⑤ 名字层（逐号）：目录记的名 vs MycoBank 接受名 —— 异名/并名的直接证据。

    判定统一走网页的 name_relation()，**不在这里另写一套**：
    两套判定漂移过一次（网页说「一致」、汇总说「异名」），修起来极难发现。

    归类口径：
      一致                        → 不算发现，跳过
      基名属 / 异名 / 种下单元      → 非错（看起来像错，实际不是）
      两皆接受名（须查原始出处）    → 不确定（两个都现行有效，须查原始发表）
    """
    rows = []
    p = os.path.join(out_dir, "{}_三列表.tsv".format(genus.lower()))
    if not os.path.exists(p):
        return rows
    try:
        from typeaudit import sources as S
    except Exception:
        return rows
    srcs = S.load(os.path.join(repo, "config", "sources.json"))
    mod = _load_pw()
    if mod is None or not hasattr(mod, "name_relation"):
        return rows
    seen = set()
    for line in io.open(p, encoding="utf-8-sig").read().splitlines()[1:]:
        f = line.split("\t")
        if len(f) < 3:
            continue
        sp = f[0].strip()
        for col in (1, 2):
            for cell in f[col].split(" ｜ "):
                m = re.match(r"^(CBS\s*[:#]?\s*[\w\.\-]+)\s*\[", cell.strip(), re.I)
                if not m:
                    continue
                item = m.group(1).strip()
                k = (item.upper().replace(" ", ""), sp)
                if k in seen:
                    continue
                seen.add(k)
                try:
                    r = mod.name_relation(item, sp, srcs) or {}
                except Exception:
                    continue
                v = (r.get("verdict") or "").strip()
                if not v or v.startswith("一致"):
                    continue
                # 同名异物 / 已核 这两类裁决是「种级」判断，同一现象已由
                # from_fourway 的「异名/独立种分歧」记过一条；此处再按号记一次
                # 会让 4 个种各多出 1 条重复（实测 1,596 → 1,600）。
                if v.startswith("同名异物") or v.startswith("已核"):
                    continue
                cls = "非错" if "非错" in v else "不确定"
                rows.append({"层级": "名字", "类别": v,
                             "对象": "{} / {}".format(sp, item),
                             "发现": "声明为 {}；该机构记名 {}；MycoBank 接受名 {}".format(
                                 sp.split()[-1],
                                 (r.get("catalog_name") or "—"),
                                 (r.get("current_name") or "—")),
                             "归类": cls, "来源": "名字层·逐号",
                             "证据": "该机构 status/taxon_name + MycoBank.current_name"})
    return rows



def from_reverse(classify):
    """③ 记录层：反向审计里「GenBank 无声明」与「以他号声明」。

    真实列名（实测）：种名 / 目录号 / 目录说的类型地位 / 检索命中 / 实证核验 / 结论 / 说明
    """
    rows = []
    p = os.path.join(REPO, "out_reverse", "反向审计.tsv")
    if not os.path.exists(p):
        return rows
    lines = io.open(p, encoding="utf-8-sig").read().splitlines()
    hdr = lines[0].split("\t")
    for l in lines[1:]:
        f = l.split("\t")
        if len(f) < len(hdr):
            continue
        r = dict(zip(hdr, f))
        verdict = (r.get("结论") or "").strip()
        if verdict not in ("GenBank 无声明", "该种类型材料以他号声明"):
            continue
        cls = "不确定"          # 未提交 ≠ 错；以他号声明也可能是同一菌的另一个号
        rows.append({"层级": "记录", "类别": verdict,
                     "对象": "{} / {}".format(r.get("种名", ""), r.get("目录号", "")),
                     "发现": ("该机构称该号是本种模式材料，但 GenBank 里查不到任何记录声明它"
                              if verdict == "GenBank 无声明"
                              else "该机构记的类型号在 GenBank 无声明；GenBank 里这个种用的是另一个号"),
                     "归类": cls, "来源": "反向审计",
                     # ⚠️ 必须把「说明」带进来：原来说明里写着"该种的类型材料记录声明了
                     #    CPC:28612（CBS 号缺失，可能是同菌异号）"—— 这是本条最有用的信息，
                     #    只抄「检索命中 0」等于把结论丢了（用户实测看不出所以然）。
                     "证据": "{} ｜ {} ｜ 检索命中 {} ｜ 实证 {}".format(
                         (r.get("该机构说的类型地位") or "")[:40],
                         (r.get("说明") or "")[:150],
                         r.get("检索命中", ""), r.get("实证核验", ""))})
    return rows



def from_seq(out_dir, classify):
    """⑤ 序列层：同一保藏号在 GenBank 里被标成几个不同的种。

    这是【直接证据】而非统计推断：同一管菌的若干条记录，organism 字段分别写着
    两个不同的种名（实测 CBS 127535：MH877934 标 boninense，JQ005199… 标 karstii）。
    同时它也是 FLAG_MULTINAME「一株多名」的真实实现 ——
    那个标记过去只定义、从未产生过。
    """
    rows = []
    # 分类库碰撞（NCBI Taxonomy 的 Type Material 字段）：同一号被 ≥2 个种登记为正模
    pc = os.path.join(out_dir, "分类库_一号多模式.tsv")
    if os.path.exists(pc):
        Lc = io.open(pc, encoding="utf-8").read().splitlines()
        if len(Lc) > 1:
            hc = Lc[0].split("\t")
            for l in Lc[1:]:
                f = l.split("\t")
                if len(f) < len(hc):
                    continue
                dc = dict(zip(hc, f))
                rows.append({"层级": "序列", "类别": "一株多名（NCBI 分类库把同号登记为两个种的正模）",
                             "对象": dc.get("号") or "",
                             "发现": "该号被 {} 个种同时登记为模式材料".format(dc.get("被几个种登记为正模")),
                             "归类": "错值", "来源": "分类库·一号多模式",
                             "证据": (dc.get("种名") or "")[:140]})
    p = os.path.join(out_dir, "序列轴_同号多名.tsv")
    if not os.path.exists(p):
        return rows
    L = io.open(p, encoding="utf-8").read().splitlines()
    if len(L) < 2:
        return rows
    h = L[0].split("\t")
    for l in L[1:]:
        f = l.split("\t")
        if len(f) < len(h):
            continue
        d = dict(zip(h, f))
        kind = (d.get("判定") or "").strip()
        cat = ("一株多名（GenBank 里被标成两个不同的种）" if kind == "两个不同的种"
               else "同种 + 种下单元（非错）")
        rows.append({"层级": "序列", "类别": cat,
                     "对象": "{} / {}".format(d.get("种") or "", d.get("号") or ""),
                     "发现": "该号在 GenBank 的记录里出现 {} 个不同的种名".format(
                         len((d.get("GenBank里出现的种名") or "").split("｜"))),
                     "归类": "错值" if kind == "两个不同的种" else "非错",
                     "来源": "序列轴·同号多名",
                     "证据": "{} ｜ 登录号 {}".format(
                         (d.get("GenBank里出现的种名") or "")[:90],
                         (d.get("登录号") or "")[:80])})
    return rows



def reclassify_missing(out_dir):
    """把「号不存在」按 GenBank 反查结果改写。

    【为什么要改】全属实测：77 条「号不存在」里，**77 条**在 GenBank 都有记录
    引用该号、且种名相符（如 CBS 126514 有 2 条、CBS 14231 有 224 条）。
    也就是说号是真实存在的，只是该机构的【公开目录没收录它】。
    把这种情况叫「号不存在」是错的 —— 它属于数据可得性问题（不确定），不是记录错。

    判据：核对1_号不存在_GenBank反查.csv 里 GenBank命中数 > 0 且种名相符。
    """
    import csv as _csv
    p = os.path.join(out_dir, "核对1_号不存在_GenBank反查.csv")
    if not os.path.exists(p):
        return None
    try:
        rows = list(_csv.DictReader(io.open(p, encoding="utf-8-sig")))
    except Exception:
        return None
    fix = {}
    for r in rows:
        try:
            n = int(r.get("GenBank命中数") or 0)
        except Exception:
            n = 0
        if n > 0 and "种名相符" in (r.get("判定") or ""):
            fix[(r.get("声明种", "").strip(), r.get("号", "").strip().upper().replace(" ", ""))] = n
    return fix


def from_quality(classify):
    """④ 源库层：quality_report.json 里 CBS 记录自身的聚合发现。"""
    rows = []
    p = os.path.join(REPO, "out", "quality_report.json")
    if not os.path.exists(p):
        return rows
    d = json.load(io.open(p, encoding="utf-8"))
    META = {
        "F14": ("CBS 记作 Glomerella 属名", "非错", "有性型属名，历史写法"),
        "F15": ("受控词表拼写错误（Oomyceet）", "错值", "CBS 自身拼写错误"),
        "F16": ("CBS 名的属 ≠ MycoBank 接受名的属", "错值", "两家对属的归属不一致"),
        "F17": ("CBS 大量条目未标模式地位", "非错", "标注缺失，非错"),
        "F18": ("正名占位符被当数据（空数组）", "错值", "数据格式缺陷"),
        "F19": ("编号写法撞车（100.06 与 10006）", "错值", "归一化冲突"),
    }
    for k, v in d.items():
        if not isinstance(v, dict) or k not in META:
            continue
        n = v.get("n")
        desc, cls, note = META[k]
        rows.append({"层级": "源库", "类别": desc, "对象": "CBS 记录",
                     "发现": "{}（命中 {} 处）".format(note, n), "归类": cls,
                     "来源": "quality_checks·{}".format(k),
                     "证据": "out/quality_report.json ｜ {} ｜ n={}".format(
                         (v.get("check") or "")[:60], n)})
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="out_audit")
    ap.add_argument("--genus", default="colletotrichum")
    ap.add_argument("--dst", default=None, help="输出路径（缺省 <out>/<genus>_发现汇总.tsv）")
    a = ap.parse_args()

    out_dir = a.out if os.path.isabs(a.out) else os.path.join(REPO, a.out)
    classify = load_flag_classes()
    rows = []
    rows += from_table(out_dir, a.genus, classify)
    rows += from_species(out_dir, classify)
    rows += from_names(out_dir, a.genus, REPO, classify)

    # 【号不存在 → 机构目录查无】用 GenBank 反查结果改写（见 reclassify_missing 说明）
    try:
        _fix = reclassify_missing(out_dir)
        if _fix:
            _n = 0
            for r in rows:
                if "号不存在" not in (r.get("类别") or ""):
                    continue
                obj = (r.get("对象") or "")
                sp = obj.split(" / ")[0].strip() if " / " in obj else ""
                num = obj.split(" / ")[-1].strip().upper().replace(" ", "") if " / " in obj else ""
                hit = _fix.get((sp, num))
                if hit:
                    r["类别"] = "机构目录查无（GenBank 有 {} 条记录引用）".format(hit)
                    r["归类"] = "不确定"
                    r["发现"] = "机构公开目录查无该号，但 GenBank 有 {} 条记录引用它且种名相符 → 号真实存在，属目录未收录".format(hit)
                    _n += 1
            if _n:
                print("  号不存在改写：{} 条 → 机构目录查无（GenBank 有引用）".format(_n))
    except Exception as _e:
        print("  ⚠️ 号不存在改写失败：{}".format(_e))
    rows += from_seq(out_dir, classify)
    rows += from_reverse(classify)
    rows += from_quality(classify)

    # 去重：同一【对象 + 类别主干】只留一条。
    # 为什么需要：管线在号层已对某些号给出「两皆接受名」标记，名字层会再算一次，
    # 于是同一条事实出现两遍（实测 6 个号各出现两次，12 条）。类别主干 = 去掉括号里的参数。
    seen = set()
    dedup = []
    for r in rows:
        base = (r.get("类别") or "").split("（")[0].strip()
        key = ((r.get("对象") or "").upper().replace(" ", ""), base)
        if key in seen:
            continue
        seen.add(key)
        dedup.append(r)
    if len(dedup) != len(rows):
        print("  去重：{} → {} 条（同一对象+同类别的重复合并）".format(len(rows), len(dedup)))
    rows = dedup

    # 【归类统一】与网页同一口径（2026-10-20 用户裁定）：
    # 除「未核 / 无目录」外一律判错值；「非错 / 矛盾」保持不变。
    # 放在这里而不是各层里：各层赋值点有 5 处，逐个改必然漏（实测踩过同类问题）。
    try:
        import importlib.util as _ilu3
        _s3 = _ilu3.spec_from_file_location("_pwcls", os.path.join(HERE, "pipeline_web.py"))
        _m3 = _ilu3.module_from_spec(_s3); _s3.loader.exec_module(_m3)
        _cls3 = _m3.classify_flag
    except Exception:
        _cls3 = None
    if _cls3 is not None:
        n_chg = 0
        for r in rows:
            old = (r.get("归类") or "").strip()
            if old in ("非错", "矛盾"):
                # 非错：各层已判「我们确定它没错」，不要回退成错值
                # 矛盾：两层名录相反，保留"方向待裁决"的性质（同样进确定错）
                continue
            new = _cls3(r.get("类别") or "")
            if new != old:
                n_chg += 1
                r["归类"] = new
        if n_chg:
            print("  归类统一：{} 条按「除未核/无公开数据库外一律错值」重定".format(n_chg))

    for i, r in enumerate(rows, 1):
        r["序号"] = i
    dst = a.dst or os.path.join(out_dir, "{}_发现汇总.tsv".format(a.genus.lower()))
    with io.open(dst, "w", encoding="utf-8", newline="") as f:
        f.write("\t".join(COLS) + "\n")
        for r in rows:
            f.write("\t".join(str(r.get(c, "")).replace("\t", " ").replace("\n", " ")
                              for c in COLS) + "\n")
    from collections import Counter
    print("汇总 {} 条 → {}".format(len(rows), dst))
    print("  按层级:", dict(Counter(r["层级"] for r in rows)))
    print("  按归类:", dict(Counter(r["归类"] for r in rows)))
    print("\n  按类别:")
    for k, v in Counter(r["类别"] for r in rows).most_common(14):
        print("    {:<40} {}".format(k[:40], v))
    return 0


if __name__ == "__main__":
    sys.exit(main())
