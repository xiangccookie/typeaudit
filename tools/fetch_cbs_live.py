# -*- coding: utf-8 -*-
"""
fetch_cbs_live.py —— 从 CBS 官方后端实时拉取菌株目录，并与本地快照对比。

【为什么需要它】
    本地 data/cbs_catalogue.tsv 是 2026-10-11 的历史快照。
    CBS 网站改版成 Angular 单页应用（SPA）后网页抓不到数据，
    但它的后端 BioloMICS REST API 仍然开放匿名访问 —— 本工具走这条路。

【破解要点：这两个坑踩了很久，务必别改回去】
    坑一：SortColumn 必须非空！
          传 "" 时服务端返回 0 条记录，不报错也不提示，极难发现。
          必须传 SortColumn="name"。
    坑二：Fields 里必须写「长数字 FieldKey」（如 14682616000000219）。
          写网页源码里的短名（e13 / name / trlink1634）服务端不认：
          记录照样返回，但 ListFields 是空的，字段值一个都不给。

【用法】
    python tools/fetch_cbs_live.py --check              # 只问「现在有多少条」
    python tools/fetch_cbs_live.py --limit 3000         # 小样本试跑
    python tools/fetch_cbs_live.py                      # 拉全量 + 与本地快照对比
    python tools/fetch_cbs_live.py --no-diff            # 只拉，不对比

【产物】
    cache/cbs_live/page_XXXXX.json   每页原始响应（断线重跑自动续传）
    out/cbs_live_<日期>.tsv          拉下来的目录
    out/cbs_live_diff_<日期>.tsv     与本地快照的差异（新增号 / 消失号 / 状态变化）
"""
import sys
try:
    sys.stdout.reconfigure(encoding="utf-8")   # Windows 控制台默认 GBK，中文会炸
except Exception:
    pass

import argparse, ast, gzip, io, json, os, re, ssl, time
import hashlib
from collections import Counter, defaultdict
import urllib.request, urllib.error, http.client

# ==================== 接口常量 ====================
API_BASE = "https://webservices.bio-aware.com/cbsdatabase"
WEBSITE_ID = "88"            # CBS 那个网站在 BioAware 上的租户号
TABLE_KEY = "14682616000000011"   # "CBS strains" 表的键
ANY_TEXT_FIELD = -102             # 系统内置的「任意文本字段」——服务端全文搜索用它

# 字段表：(输出列名, FieldKey)。FieldKey=None 表示值在返回的 Name / Id 里。
FIELDS = [
    ("cbs_number",   None),                    # 值在 Name，如 "CBS 100.06"
    ("status",       "14682616000000219"),     # Status of the strain —— 模式菌株状态
    ("taxon_full",   "14682616000000069"),     # Taxon name new —— 含作者 + [MB#xxxxx]
    ("organism",     "14682616000001122"),     # Type of organism —— 酵母/丝状真菌
    ("other_coll",   "14682616000000218"),     # Other collections —— 其他保藏机构号
    ("name_changes", "14682616000000228"),     # Name changes —— 名称变更史
]
FIELD_KEYS = [fk for _, fk in FIELDS if fk]

PAGE_SIZE = 2000      # 实测 2000 可行（11,573,012 字节 / 19.42 秒）；
                      # 5000 会返回 36 MB 被掐断（IncompleteRead）
REQ_SLEEP = 1.3       # 限流：每次请求之间歇一下，别把人家服务器打疼

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# ==================== 取数纪律（请遵守，这是对 CBS 服务器的基本尊重）====================
# 1) 诚实标明身份和联系方式 —— 万一给人家添了麻烦，人家能找到你。这既是礼貌也是自保。
# 2) 全量只拉一次。手上已有 2026-10-11 的快照，加这一次就是两个时间点，够做完审稿。
# 3) cache/ 里的每页响应别删。重跑直接读缓存，不重复请求人家服务器。
# 4) 想更温和：--sleep 3（每次请求间隔 3 秒）或 --page 500（每页小一点）。
# 5) 如果哪天接口关了或要认证 —— 那就是「不」，不要绕。不要试图绕过认证。
#
# 实测负担参考（本机实测，可复核）：
#   每页 2000 条 → 全量 32 次请求、约 346 MB、服务端渲染约 621 秒
#   对照人工在网页上翻完同一个表：每页 100 条需点 628 次「下一页」
CONTACT_EMAIL = "your.name@example.com"   # ← 改成你自己的邮箱再跑（不改也能跑，但请改）

HEADERS = {
    "User-Agent": "Mozilla/5.0 (academic metadata audit)",
    "WebsiteId": WEBSITE_ID,          # ← 少了这个头，所有接口都 401
    "Accept": "application/json",
    "Content-Type": "application/json",
    "Referer": "https://wi.knaw.nl/",
}

_CTX = ssl.create_default_context()
_CTX.check_hostname = False
_CTX.verify_mode = ssl.CERT_NONE

# 本次进程实际向服务器发出的请求统计（读缓存不算）—— 用来如实交代「打了几次、拿了多少」
_STATS = {"requests": 0, "bytes": 0, "seconds": 0.0}


# ==================== 请求层 ====================
def _decode(resp):
    """读响应体；如果服务端做了 gzip 压缩就先解开。"""
    raw = resp.read()
    if resp.headers.get("Content-Encoding") == "gzip":
        raw = gzip.decompress(raw)
    return raw.decode("utf-8", "replace")


def query(fields, start=0, length=PAGE_SIZE, sort_col="name", sort_dir="asc",
          retries=4, verbose=True, conditions=None, complex_query=None):
    """
    发一次 SearchForSummaryGrid 请求。

    fields   : FieldKey 列表（长数字！）
    start    : 从第几条开始（分页用）
    length   : 要几条
    sort_col : 【必须非空】否则服务端返回 0 条
    返回     : (解析后的 JSON, 收到的字节数, 耗时秒)
    """
    body = {
        "TableKey": TABLE_KEY,
        "Fields": list(fields),
        "iDisplayLength": length,
        "iDisplayStart": start,
        "ListOfConditionEntity": list(conditions or []),   # 过滤条件，空 = 不过滤
        "ListOfOrderEntity": [],
        # ComplexQuery 是「条件编号表达式」：Q0 指第 0 个条件，多条件写成 "Q0 Or Q1"。
        # 有过滤条件时必须给非空且【列全】的表达式，否则条件会被静默忽略（返回全量）。
        "ComplexQuery": (complex_query if complex_query is not None
                         else (" Or ".join("Q{}".format(i)
                                           for i in range(len(conditions))) if conditions else "")),
        "SortColumn": sort_col,          # ← 坑一：留空就返回 0 条
        "SortDirection": sort_dir,
        "LoadOwnerRecord": False,
    }
    last_err = None
    ask = length          # 本次实际要多少条；失败后减半，别让服务器重做同样大的活
    MIN_ASK = 200
    for attempt in range(retries):
        try:
            body["iDisplayLength"] = ask
            req = urllib.request.Request(
                API_BASE + "/api/Search/SearchForSummaryGrid",
                headers=HEADERS, data=json.dumps(body).encode(), method="POST")
            t0 = time.time()
            with urllib.request.urlopen(req, timeout=240, context=_CTX) as r:
                txt = _decode(r)
            dt = time.time() - t0
            _STATS["requests"] += 1
            _STATS["bytes"] += len(txt)
            _STATS["seconds"] += dt
            return json.loads(txt), len(txt), dt
        except (http.client.IncompleteRead, urllib.error.URLError,
                TimeoutError, OSError, json.JSONDecodeError) as e:
            last_err = e
            if verbose:
                nxt = max(MIN_ASK, ask // 2)
                print("     [重试 {}/{}] {} —— 页大小 {} → {} 条，退避 {} 秒".format(
                    attempt + 1, retries, type(e).__name__, ask, nxt, 3 + attempt * 3))
            ask = max(MIN_ASK, ask // 2)   # 减半重试：既更可能成功，也让服务器少做无用功
            time.sleep(3 + attempt * 3)
    raise last_err


def dig(v):
    """
    拆嵌套值。
    有些字段（正名、国家）返回的是一段 Python 字面量文本，形如：
        [{'TemplateId': 0, 'Name': {'Value': 'Absidia glauca'}}]
    我们要的是里面所有 Value 拼起来。
    """
    if v is None:
        return ""
    s = str(v).strip()
    if s[:2] in ("[{", "{'"):
        try:
            obj = ast.literal_eval(s)
            acc = []

            def walk(o):
                if isinstance(o, dict):
                    for k, vv in o.items():
                        if k == "Value" and not isinstance(vv, (dict, list)):
                            if str(vv).strip():
                                acc.append(str(vv))
                        else:
                            walk(vv)
                elif isinstance(o, (list, tuple)):
                    for vv in o:
                        walk(vv)

            walk(obj)
            return " | ".join(acc) if acc else s
        except Exception:
            return s
    return s


def to_row(rec):
    """把一条 API 记录拍平成一行字典。"""
    row = {"id": rec.get("Id"), "cbs_number": (rec.get("Name") or "").strip()}
    for lf in (rec.get("ListFields") or []):
        fv = lf.get("FieldValue") or {}
        raw = fv.get("DataToView")
        if raw in (None, ""):
            raw = fv.get("Value")
        row[(lf.get("Name") or "").strip()] = dig(raw)
    return row


MB_RE = re.compile(r"\[MB#(\d+)\]")


def parse_taxon(taxon_full):
    """
    把 "Absidia glauca | Absidia glauca Hagem, ... (1908) [MB#221208]"
    拆成 (纯名, MycoBank号)。MycoBank 号是审稿时做名称对照的钥匙。
    """
    s = taxon_full or ""
    first = s.split("|")[0].strip()
    m = MB_RE.search(s)
    return first, (m.group(1) if m else "")


def build_conditions(genus, synonyms=None):
    """
    构造查询条件列表 + ComplexQuery 表达式。属名与同义属放在一起做【或】查询。

    【为什么需要同义属 —— 实测出来的召回缺口】
        只按 "Colletotrichum" 一个词筛，会漏掉 CBS 挂在【有性型属名】下的记录：
        CBS 113010 / 114054 / 116868 / 116869 / 118192 / 122487 六条记作
        "Glomerella cingulata"，而 MycoBank 当前名是 "Colletotrichum gloeosporioides"。
        它们通篇文本里不出现 "Colletotrichum"，所以按属名单词筛一条都捞不到
        （实测漏 6 / 856 = 0.70%）。

    【三条必须同时满足，少一条服务端就拒或不认】
      1) 每个条件都要有 Index，从 0 递增。
      2) 第一个条件的 OperatorOfCondition 是空串 ""，后面的填 "Or"（或 "And"）。
      3) ComplexQuery 必须把【所有】条件都列出来："Q0 Or Q1 Or Q2"。

    实测数据（Colletotrichum 850 条、Glomerella 1333 条，两者差 498 条）：
        两个条件带 Index + "Q0 Or Q1"  → 1348  ← 正好是并集（850 + 498）
        两个条件带 Index + "Q0 And Q1" → 835   ← 正好是交集
        条件缺少 Index                 → Success:false，整条查询被拒
        ComplexQuery 给空串            → 62751，条件被【静默忽略】返回全量
        两个条件但 ComplexQuery 只写 Q0 → HTTP 503

    条件里的操作符名 "Contains" 大小写敏感，写小写会被静默忽略（见下）。
    """
    names = [genus] + [x for x in (synonyms or []) if x]
    conds = []
    for i, nm in enumerate(names):
        conds.append({
            "Index": i,
            "FieldKey": ANY_TEXT_FIELD,        # -102 = 任意文本字段（全文搜索）
            "Value": nm,
            "OperatorOfField": "Contains",     # 大小写敏感！小写会被静默忽略
            "OperatorOfCondition": "" if i == 0 else "Or",
            "NotOfCondition": False,
        })
    cq = " Or ".join("Q{}".format(i) for i in range(len(conds)))
    return conds, cq


def is_target_genus(taxon_name, genus):
    """客户端精确判属：正名的第一个词是否就是该属。"""
    nm = (taxon_name or "").strip()
    if not nm:
        return False
    return nm.split()[0].lower() == (genus or "").strip().lower()


def norm_number(code):
    """
    归一化 CBS 保藏号 —— 【必须保留小数点】。判重、比对一律用它。

    【为什么绝不能把小数点删掉】
        CBS 有两套并行编号体系：
            CBS 100.06  = 100 号菌株的第 06 支（老体系，带分号）
            CBS 10006   = 10006 号菌株（新体系，整号）
        这两条是【完全不同的菌株】。本机实测 5,000 条里有 426 组「删点后撞车」，
        逐组核对正名，没有任何一组正名相同（真重复 = 0 组）。
        也就是说：删点归一化会凭空制造 426 个假重复。
        另有 CBS 101.77A / 101.77B 这种带字母后缀的写法，同理。
    """
    return re.sub(r"\s+", "", (code or "").upper())


def norm_loose(code):
    """
    宽松归一化：连小数点一起删掉。
    【只用来做「编号写法撞车」观察，绝不可用来判重】—— 它会把 CBS 100.06 和
    CBS 10006 合成一个键，那是两株不同的菌。
    """
    return re.sub(r"[^A-Z0-9]", "", (code or "").upper())


# ==================== 主流程 ====================
def selfcheck(rows, total, verbose=True):
    """
    分页拉取自检：漏条 / 重条，另附「编号写法撞车」观察。

    【为什么必须查漏条重条】
        分页靠 SortColumn="name" 排序 + iDisplayStart 偏移量定位下一批。
        如果两次请求之间某条记录的 CBS 号被改动（改号或被并进另一条），
        它在排序里的位置就变了 —— 结果可能漏掉一条，也可能同一条被抓两次。
        所以「拉到多少行」本身不能证明拉全了，必须对账。

    【为什么另外报「撞车」】
        CBS 同时存在 CBS 100.06（100 号第 06 支）与 CBS 10006（10006 号）两种写法。
        删掉小数点后两者会合并成同一个键。这不是数据库错误，是编号体系本身的特点，
        但任何做保藏号匹配的人都会踩到 —— 值得作为【观察】记进审稿日志。
    """
    keys = [k for k in (norm_number(r.get("cbs_number", "")) for r in rows) if k]
    cnt = Counter(keys)
    dup = {k: v for k, v in cnt.items() if v > 1}

    # 撞车观察：宽松键相同、但原始写法不同
    loose = defaultdict(set)
    for r in rows:
        loose[norm_loose(r.get("cbs_number", ""))].add(r.get("cbs_number", ""))
    collide = {k: sorted(v) for k, v in loose.items() if len(v) > 1}

    missing = max(0, total - len(rows))
    res = {
        "rows": len(rows),
        "total_reported": total,
        "unique_keys": len(cnt),
        "missing": missing,
        "n_duplicate_keys": len(dup),
        "duplicate_sample": sorted(dup.items())[:20],
        "n_number_collisions": len(collide),
        "collision_sample": sorted(collide.items())[:8],
        "clean": (len(rows) == total) and (not dup),
    }
    if verbose:
        print()
        print("=" * 90)
        print("分页拉取自检（对账：漏条 / 重条）")
        print("=" * 90)
        print("  拉到行数              {:,}".format(res["rows"]))
        print("  服务端自报总数        {:,}".format(total))
        print("  去重后的保藏号数      {:,}".format(res["unique_keys"]))
        print("  漏条                  {:,}".format(missing))
        print("  重复的保藏号          {:,}".format(res["n_duplicate_keys"]))
        if missing:
            print("    ! 少了 {:,} 条 —— 拉取窗口期内数据库可能被改动过".format(missing))
        if dup:
            print("    ! 下列号出现多次（前 20 个）：")
            for k, v in res["duplicate_sample"]:
                print("        {:<20} 出现 {} 次".format(k, v))
        print("  结论                  {}".format(
            "干净 —— 行数、去重数、服务端总数三者一致"
            if res["clean"] else
            "不干净 —— 见上面 ! 行。这不是工具 bug，是数据在动的证据：先记入审稿日志，再重拉那一段"))

        print()
        print("  附：编号写法撞车观察（不是错误，是 CBS 编号体系的特点）")
        print("      把小数点删掉后会同键的写法组数：{:,}".format(res["n_number_collisions"]))
        for k, v in res["collision_sample"]:
            print("        宽松键 {:<16} ← {}".format(k, v))
        if res["n_number_collisions"]:
            print("      含义：按保藏号做任何匹配时【必须保留小数点】，否则这批数据里会凭空")
            print("            多出 {:,.0f} 组假重复。可写进论文的「数据使用注意」一节。".format(
                res["n_number_collisions"]))
    return res


def genus_from_status(status):
    """
    从 status 字符串里推出属名。
        "Type of Colletotrichum hymenocallidis ..."  →  "Colletotrichum"

    【为什么需要它】
        本地快照 cbs_catalogue_raw.jsonl 的 taxon 字段是【空的】
        （实测 62,751 行全为空值），能用的属名线索只剩 status。
        代价：status 只覆盖模式菌株（实测 15,109 / 62,751 = 24.1%），
        所以这样得到的属级清单是【部分覆盖】，报告里必须说明。
    """
    st = status or ""
    low = st.lower()
    for kw in ("type of ", "holotype of ", "neotype of ", "isotype of ",
               "ex-type of ", "syntype of ", "lectotype of ", "paratype of ",
               "authentic strain of ", "epitype of "):
        i = low.find(kw)
        if i >= 0:
            rest = st[i + len(kw):].strip()
            return rest.split()[0] if rest else ""
    return ""


def load_local_snapshot(path, genus=None):
    """
    读本地快照 → {严格键: 行字典}。
    支持 .jsonl（列名 name/id/status/taxon/organism）与 .tsv（cbs_number/...）。
    genus 非空时只保留该属的行 —— 比对前必须把两边收窄到同一个属，
    否则线上只有该属、本地是全量，会凭空报出几万个「消失」。

    返回 (字典, 覆盖说明)。覆盖说明非 None 时，表示本地这一侧是【部分覆盖】。
    """
    loc = {}
    note = None
    if path.lower().endswith(".jsonl"):
        with io.open(path, encoding="utf-8", errors="replace") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                try:
                    d = json.loads(line)
                except Exception:
                    continue
                num = d.get("name") or d.get("cbs_number") or ""
                tax = d.get("taxon") or ""
                if genus:
                    if tax.strip():
                        if genus.lower() not in tax.lower():
                            continue
                    else:
                        # 本地 taxon 为空 → 退回从 status 推属名（只覆盖模式菌株）
                        if genus_from_status(d.get("status")).lower() != genus.lower():
                            continue
                        note = ("本地快照没有 taxon 字段，属名是从 status 推出来的，"
                                "因而只覆盖模式菌株（本地 status 有值 15,109 / 62,751 = 24.1%）")
                k = norm_number(num)
                if k:
                    loc[k] = {"cbs_number": num, "status": d.get("status") or "", "taxon": tax}
        return loc, note

    with io.open(path, encoding="utf-8", errors="replace") as f:
        head = f.readline().rstrip("\n").split("\t")
        if "cbs_number" not in head:
            print("  本地快照缺少 cbs_number 列（现有列：{}），无法用同一把尺子比对".format(head))
            return None, None
        if genus and "taxon" not in head:
            print("  本地快照没有 taxon 列，无法把范围收窄到 {}，跳过对比".format(genus))
            print("  建议用带 taxon 列的原始快照：cbs_catalogue_raw.jsonl")
            return None, None
        for line in f:
            pr = line.rstrip("\n").split("\t")
            if len(pr) < 2:
                continue
            d = dict(zip(head, pr))
            if genus and genus.lower() not in (d.get("taxon") or "").lower():
                continue
            # 本地快照里的 norm_key 是旧版「删点」键，会撞车，不能用；
            # 从原始 cbs_number 现算严格键，两边才是同一把尺子。
            k = norm_number(d.get("cbs_number", ""))
            if k:
                loc[k] = d
    return loc, note


def coverage_report(genus):
    """
    跑一个属【之前】先看数据够不够。

    【为什么需要它】
        名称轴的判定靠 CBS 目录的 status 字段（里面写着模式地位 + 名字）。
        但实测 CBS 目录 62,740 行里只有 15,098 行（24.1%）有标注，
        能推出属名的 14,763 行只涉及 2,535 个属 —— 而且极不均匀：
        ≥100 条标注的只有 15 个属，≥10 条的 282 个。
        冷门属（地衣、蘑菇、锈菌、大量小属）**一条标注都没有**，
        跑下去只会得到满屏「未核」。这个自检就是提前把这件事说清楚。

    判据（N = CBS 目录里 status 含该属名的行数）：
        N = 0     名称轴不可用 —— 该属的号只能得到「未核/号不存在」
        N < 10    名称轴可用但命中有限
        N < 100   名称轴可用
        N >= 100  名称轴强（全库只有 15 个属到这个量级）
    """
    sys.path.insert(0, os.path.join(REPO, "src"))
    print("=" * 88)
    print("覆盖率自检：{}".format(genus))
    print("=" * 88)
    try:
        from typeaudit import sources as SN
        cfg = SN.load(os.path.join(REPO, "config", "sources.json"))
    except Exception as e:
        print("  载入数据源失败：{}".format(str(e)[:80]))
        return 1
    cbs, myb = cfg.get("cbs"), cfg.get("mycobank")
    if cbs is None:
        print("  ❌ CBS 目录没装载 —— 名称轴不可用。先跑：python tools/setup_data.py")
        return 1
    if myb is None:
        print("  ⚠️ MycoBank 没装载 —— 判「异名/基名属」与接受名链会失效。"
              "先跑：python tools/setup_data.py")
    g = genus.lower()
    # CBS：status 里含属名 且 出现在模式地位关键词之后
    KWS = ("type of ", "holotype of ", "neotype of ", "isotype of ", "ex-type of ",
           "syntype of ", "lectotype of ", "paratype of ", "authentic strain of ", "epitype of ")
    n_annot = 0
    for r in (getattr(cbs, "rows", []) or []):
        st = (r.get("status") or "").lower()
        if not st:
            continue
        for kw in KWS:
            i = st.find(kw)
            if i >= 0:
                if st[i + len(kw):].strip().startswith(g + " "):
                    n_annot += 1
                break
    n_rows = len(getattr(cbs, "rows", []) or [])
    print("  CBS 目录：{:,} 行，本属有模式地位标注 {:,} 条".format(n_rows, n_annot))
    n_my = 0
    if myb is not None:
        byname = getattr(myb, "byname", {}) or {}
        n_my = sum(1 for k in byname if k.startswith(g + " "))
        print("  MycoBank：本属作为【异名】出现 {:,} 条（即有多少老名会归到别的接受名）".format(n_my))
    print()
    if n_annot == 0:
        print("  ❌ 判读：名称轴【不可用】—— 本属在 CBS 目录里没有任何模式地位标注。")
        print("     跑下去你只会看到「未核 / 号不存在」。出路：")
        print("       · 自备该属所属机构的目录：放成 data/catalogue_<机构代码>.tsv")
        print("       · 或换轴（序列 / 文献），本工具答不了这个属")
    elif n_annot < 10:
        print("  ⚠️ 判读：名称轴【可用但命中有限】—— 只有 {} 条标注。".format(n_annot))
    elif n_annot < 100:
        print("  ✅ 判读：名称轴【可用】。")
    else:
        print("  ✅ 判读：名称轴【强】（全库只有 15 个属 ≥100 条标注）。")
    print()
    print("  记住边界：本工具只核【两列保藏号本身】有没有问题，")
    print("            不核序列是否属于该种、不核菌株是否真是模式、不核论文声明。")
    return 0


def load_derived_synonyms(genus, min_n=15):
    """
    读 tools/derive_synonyms.py 派生的同义属清单（data/synonyms_<属>.tsv）。

    为什么要有这一步：同义属清单绝不能手写（本工具作者手写过一次，漏了 17 个）。
    清单必须来自可复现的反推，并且带得出「什么时候、用哪张表反推的」。
    返回 (清单, 说明文字)；文件不存在时清单为空并给出去哪生成。
    """
    # 文件名大小写不敏感地找 —— 否则在区分大小写的系统上会静默找不到清单
    data_dir = os.path.join(REPO, "data")
    path = None
    want = "synonyms_{}.tsv".format(genus).lower()
    for fn in (os.listdir(data_dir) if os.path.isdir(data_dir) else []):
        if fn.lower() == want:
            path = os.path.join(data_dir, fn)
            break
    if not path:
        return [], ("没找到派生清单 data/synonyms_{}.tsv —— 先跑：\n"
                    "      python tools/derive_synonyms.py {}".format(genus, genus))
    names, skipped, n_mutual = [], 0, 0
    with io.open(path, encoding="utf-8") as fh:
        head = fh.readline().rstrip("\n").split("\t")
        i_g = head.index("synonym_genus") if "synonym_genus" in head else 1
        i_n = head.index("n") if "n" in head else 2
        i_m = head.index("mutual") if "mutual" in head else None
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) <= i_n:
                continue
            try:
                n = int(p[i_n])
            except Exception:
                continue
            if n < min_n:
                skipped += 1
                continue
            # 注意：这里【不】因 mutual=1 而丢弃。
            # 实测 Colletotrichum 的 Gloeosporium(×37)、Vermicularia(×37) 都被标成
            # 双向互挂，可它们是最重要的真同义属 —— 互相出现只是因为确有物种被移动过。
            # 所以 mutual 只是「提示你人工看一眼」，不是丢弃依据。
            if i_m is not None and len(p) > i_m and p[i_m].strip() == "1":
                n_mutual += 1
            names.append(p[i_g].strip())
    note = ("已从 {} 自动载入 {} 个同义属（出现 ≥{} 次；跳过 {} 个低频候选）".format(
        os.path.relpath(path, REPO), len(names), min_n, skipped))
    if n_mutual:
        note += ("\n      其中 {} 个被标为「双向互挂」（{}）——".format(
            n_mutual, ", ".join(names[:4]) + ("…" if len(names) > 4 else "")))
        note += "\n      这只表示「确有物种被移动过」，不是噪声证据，已全部保留，建议人工过一遍。"
    meta_p = path[:-4] + ".meta.json"
    if os.path.exists(meta_p):
        try:
            m = json.loads(io.open(meta_p, encoding="utf-8").read())
            note += "\n      溯源：{} 反推于 {}（扫 {} 行）".format(
                os.path.basename(m.get("source_file", "?")), m.get("derived_at", "?"),
                m.get("rows_scanned", "?"))
        except Exception:
            pass
    return names, note


def resolve_local(path_opt, genus=None):
    """找本地快照。命令行给了就用给的，否则按几个已知位置找。"""
    cands = [path_opt] if path_opt else []
    if genus:
        # 属模式优先用带 taxon 列的原始快照，否则无法把本地收窄到同一个属
        cands += [
            os.path.join(REPO, "data", "cbs_catalogue_raw.jsonl"),
            r"C:\Users\Administrator\Desktop\DeepSeek\bio原理\out\cbs_catalogue_raw.jsonl",
        ]
    cands += [
        os.path.join(REPO, "data", "cbs_catalogue.tsv"),
        r"C:\Users\Administrator\Desktop\DeepSeek\bio原理\out\cbs_catalogue.tsv",
    ]
    for p in cands:
        if p and os.path.exists(p):
            return p
    return None


def pull(limit=None, verbose=True, page=PAGE_SIZE, genus=None, synonyms=None):
    """
    按页拉取，返回 (行列表, 服务端自报总数)。

    genus 非空 → 只拉该属（服务端过滤），并把缓存隔离到独立子目录。
    缓存必须隔离：否则 genus_Colletotrichum 的第 0 页会覆盖全量拉取的第 0 页，
    两者内容完全不同，会造成静默串数据。
    """
    if genus:
        tag = re.sub(r"[^A-Za-z0-9]", "_", genus)
        if synonyms:
            # 同义属清单用【短哈希】而不是把名字全拼进去：
            # 实测把 20 个属名拼进目录名后，整条路径超过 Windows MAX_PATH(260)，
            # makedirs 直接抛 FileNotFoundError（WinError 3），且报错中文还是乱码，
            # 很难看出是路径长度问题。
            sig = hashlib.sha1("|".join(sorted(synonyms)).encode("utf-8")).hexdigest()[:8]
            tag += "_syn{}_{}".format(len(synonyms), sig)
        cache_dir = os.path.join(REPO, "cache", "cbs_live", "genus_" + tag)
    else:
        cache_dir = os.path.join(REPO, "cache", "cbs_live")
    if not os.path.isdir(cache_dir):
        os.makedirs(cache_dir)
    # 属名 + 同义属一起做【或】查询，一次请求拿全（服务端支持多条件或运算）
    conds, cq = build_conditions(genus, synonyms) if genus else ([], "")

    # ---- 先问总数。首页若已在缓存里，就只要 1 条，别白拉一整页 ----
    cache0 = os.path.join(cache_dir, "page_{:06d}.json".format(0))
    have0 = os.path.exists(cache0)
    if verbose:
        if have0:
            print("  首页已在缓存中，只需取 1 条来问最新总数……", flush=True)
        else:
            # 提示必须打在请求【之前】并 flush，否则用户看到的是一片空白，
            # 分不清「程序在跑」还是「程序卡死了」。
            print("  正在请求第 1 页（{} 条）…… 首次响应通常要 20-35 秒，请耐心等".format(
                page), flush=True)
    j, nb, dt = query(FIELD_KEYS, 0, 1 if have0 else page, conditions=conds, complex_query=cq)
    total = int(((j.get("Data") or {}).get("TotalCount")) or 0)
    if not total:
        raise SystemExit("服务端报总数为 0，接口形状可能又变了（检查 SortColumn 是否为空）")
    if not have0:
        # 首页响应立刻落缓存 —— 否则循环回到 start=0 时会向服务器再要一次同一页
        io.open(cache0, "w", encoding="utf-8").write(json.dumps(j))
    if verbose:
        print("  服务端自报总数 TotalCount = {:,} 条".format(total))
        if not have0:
            print("  首页 {} 条  {:,} 字节  {:.1f} 秒".format(
                len((j.get("Data") or {}).get("RecordEntityList") or []), nb, dt))

    if limit:
        total = min(total, int(limit))
        if verbose:
            print("  --limit 生效，本次只拉前 {:,} 条".format(total))

    rows, start = [], 0
    while start < total:
        cache = os.path.join(cache_dir, "page_{:06d}.json".format(start))
        if os.path.exists(cache):                     # 断线重跑：直接读缓存
            j = json.loads(io.open(cache, encoding="utf-8").read())
            src = "缓存"
        else:
            want = min(page, total - start)   # 最后一页不足一整页时别多要
            time.sleep(REQ_SLEEP)
            if verbose:
                print("      → 正在请求第 {:,} 条起（{} 条）……".format(start + 1, want),
                      flush=True)
            j, nb, dt = query(FIELD_KEYS, start, want, conditions=conds, complex_query=cq)
            io.open(cache, "w", encoding="utf-8").write(json.dumps(j))
            src = "{:,}B {:.1f}s".format(nb, dt)
        recs = (j.get("Data") or {}).get("RecordEntityList") or []
        if not recs:
            if verbose:
                print("  第 {:,} 条起返回空，提前结束".format(start))
            break
        rows += [to_row(r) for r in recs]
        start += len(recs)
        if verbose:
            print("  已拉 {:,} / {:,} 条   （{}）".format(len(rows), total, src))
    rows = rows[:total]                       # 最后一页可能超出，截掉
    return rows, total


def main():
    global REQ_SLEEP        # 必须在函数最开头声明，否则 --sleep 默认值那行会先用到它
    ap = argparse.ArgumentParser(description="实时拉取 CBS 菌株目录并与本地快照对比")
    ap.add_argument("--check", action="store_true", help="只问总数，不拉数据")
    ap.add_argument("--limit", type=int, default=None, help="只拉前 N 条（试跑用）")
    ap.add_argument("--no-diff", action="store_true", help="只拉取，不与本地对比")
    ap.add_argument("--local", default=None, help="本地快照 tsv 路径")
    ap.add_argument("--out-dir", default=os.path.join(REPO, "out"))
    ap.add_argument("--page", type=int, default=PAGE_SIZE,
                    help="每页拉多少条（默认 {}；越大约束服务器请求次数越少）".format(PAGE_SIZE))
    ap.add_argument("--sleep", type=float, default=REQ_SLEEP,
                    help="每次请求之间的间隔秒数（默认 {}，想更温和就调大）".format(REQ_SLEEP))
    ap.add_argument("--syn-min", type=int, default=15,
                    help="自动载入同义属的最低出现次数（默认 15）。这是【省流量的成本旋钮】，"
                         "不是正确性筛子：实测 n>=15 恰好得到 Gloeosporium/Vermicularia/Glomerella 三个，"
                         "粗筛 1,365 条；降到 2 会变成 20 个、7,345 条，慢 5 倍。")
    ap.add_argument("--check-coverage", action="store_true",
                    help="只做覆盖率自检：这个属的数据够不够（不联网、不拉取）")
    ap.add_argument("--no-synonyms", action="store_true",
                    help="完全不自动载入同义属（最快；只按属名单词筛）")
    ap.add_argument("--synonyms", default=None,
                    help="同义属/有性型属清单，逗号分隔（例 Glomerella,Gloeosporium）——"
                         "不加会漏掉挂在有性型属名下的记录")
    ap.add_argument("--genus", default=None,
                    help="只拉指定属（服务端过滤，例 --genus Colletotrichum）—— 比全量小一个数量级")
    ap.add_argument("--email", default=None,
                    help="写进 User-Agent 的联系邮箱（诚实标明身份用）")
    a = ap.parse_args()

    # 覆盖率自检是【只读】操作：先返回，免得打印取数身份/页大小那套毫无关系的信息
    if a.check_coverage:
        if not a.genus:
            print("请同时给 --genus，例如：--genus Colletotrichum --check-coverage")
            return 1
        return coverage_report(a.genus)


    # ---- 把「礼貌」落实到本次运行 ----
    REQ_SLEEP = a.sleep
    mail = a.email or CONTACT_EMAIL
    HEADERS["User-Agent"] = ("colleto-type-audit/1.0 "
                             "(academic metadata review; contact: {})".format(mail))
    print("  本次取数身份：{}".format(HEADERS["User-Agent"]))
    print("  每页 {} 条 / 请求间隔 {} 秒".format(a.page, a.sleep))
    if mail.startswith("your.name"):
        print("  ⚠ 提醒：还没设置联系邮箱，建议用 --email 你的邮箱@单位 跑一次（对人家礼貌）")

    print("=" * 90)
    print("CBS 目录实时拉取（BioAware / BioloMICS REST API）")
    print("=" * 90)

    if a.check:
        print("  正在查询……", flush=True)
        j, nb, dt = query(FIELD_KEYS, 0, 10)
        d = j.get("Data") or {}
        print("  服务端报总数：{:,} 条".format(int(d.get("TotalCount") or 0)))
        print("  首页样例：")
        for r in (d.get("RecordEntityList") or [])[:3]:
            row = to_row(r)
            nm, mb = parse_taxon(row.get("Taxon name new", ""))
            print("     {:<14} {:<28} {}  MB#{}".format(
                row.get("cbs_number", ""), nm[:28], row.get("Status of the strain", "")[:34], mb))
        return

    syns = [x.strip() for x in (a.synonyms or "").split(",") if x.strip()]
    if a.genus and not syns and not a.no_synonyms:
        # 没手给同义属 → 自动读派生清单（手写清单是不可复现的，必须走这条路）
        auto, note = load_derived_synonyms(a.genus, min_n=a.syn_min)
        print("  " + note)
        if auto:
            syns = auto
    elif syns:
        print("  同义属（手工指定）：{}".format(", ".join(syns)))
    rows, total = pull(limit=a.limit, page=a.page, genus=a.genus, synonyms=syns)
    print("\n  共拉到 {:,} 行".format(len(rows)))

    # 对账自检：漏条 / 重条
    check = selfcheck(rows, total)

    # ---- 属别判定：服务端是全文粗筛，必须按正名分开两类 ----
    if a.genus:
        n_target = n_syn = n_other = 0
        syn_low = set(x.lower() for x in syns)
        counts = {}
        for r in rows:
            nm = (r.get("Taxon name new") or "").split("|")[0].strip()
            first = (nm.split()[0] if nm else "").lower()
            if first == a.genus.lower():
                n_target += 1
            elif first in syn_low:
                n_syn += 1
                counts[first] = counts.get(first, 0) + 1
            else:
                n_other += 1
        print()
        print("=" * 90)
        print("属别判定（服务端是全文粗筛，这里按正名分开）")
        print("=" * 90)
        print("  正名就是 {:<26} {:>6} 条".format(a.genus, n_target))
        if syns:
            for sy in syns:
                print("  正名是同义属名 {:<18} {:>6} 条".format(sy, counts.get(sy.lower(), 0)))
            print("      ↑ 这一批就是【只按单词筛会漏掉】的记录")
        print("  正名是别的属（文本里提到而已）      {:>6} 条".format(n_other))
        print("      ↑ 粗筛的假阳性，做属级清单时必须丢弃")
        print()
        print("  正名确属该属（可直接用）   ：{:,} 条".format(n_target))
        if n_syn:
            print("  同义属名下、待人工定属     ：{:,} 条".format(n_syn))
            print("      工具【不】替你判定这 {} 条属不属于该属 —— 因为它确实不能自动判。".format(n_syn))
            print("      实测理由：Colletotrichum 的同义属名那 22 条里，MycoBank 只认 6 条")
            print("      （Glomerella cingulata → Colletotrichum gloeosporioides），")
            print("      另 16 条 MycoBank 自己也仍用原属名（Glomerella tucumanensis × 15、")
            print("      Glomerella lagenaria、Gloeosporium solani）—— 属【权威分歧】，")
            print("      得人工按你采用的那个权威定，不能当成已确属自动算进去。")
        print("  粗筛上限（全算入，仅供参照） ：{:,} 条".format(n_target + n_syn))

    # 如实交代本次打了几次、拿了多少（读缓存的不算在内）
    print()
    print("  本次实际向服务器发出 {} 次请求，共收到 {:,} 字节，服务端渲染累计 {:.0f} 秒".format(
        _STATS["requests"], _STATS["bytes"], _STATS["seconds"]))

    # ---- 落盘 ----
    if not os.path.isdir(a.out_dir):
        os.makedirs(a.out_dir)
    stamp = time.strftime("%Y%m%d")
    if a.genus:
        out_tsv = os.path.join(a.out_dir, "cbs_{}_{}.tsv".format(
            re.sub(r"[^A-Za-z0-9]", "", a.genus), stamp))
    else:
        out_tsv = os.path.join(a.out_dir, "cbs_live_{}.tsv".format(stamp))
    # ⚠️ 打印绝对路径。相对路径会误导两拨人：
    #    ① 人：复制粘贴到别处就找不到；
    #    ② 程序：sources.load() 解析相对路径是按【配置目录】找的，
    #       找不到就静默回退到自动探测的那份 —— 用户以为在用刚下的，其实用的是别的。
    out_tsv = os.path.abspath(out_tsv)
    cols = ["cbs_number", "id", "Status of the strain", "Taxon name new",
            "Type of organism", "Other collections", "Name changes"]
    with io.open(out_tsv, "w", encoding="utf-8", newline="") as f:
        f.write("cbs_number\tid\tstatus\ttaxon_full\ttaxon_name\tmb_id\torganism\tother_coll\tname_changes\tnumber_key\n")
        for r in rows:
            nm, mb = parse_taxon(r.get("Taxon name new", ""))
            f.write("\t".join(str(x).replace("\t", " ").replace("\n", " ") for x in [
                r.get("cbs_number", ""), r.get("id", ""),
                r.get("Status of the strain", ""), r.get("Taxon name new", ""),
                nm, mb, r.get("Type of organism", ""),
                r.get("Other collections", ""), r.get("Name changes", ""),
                norm_number(r.get("cbs_number", "")),
            ]) + "\n")
    print("  已写 {}".format(out_tsv))

    if a.no_diff:
        return

    # ---- 与本地快照对比 ----
    local = resolve_local(a.local, a.genus)
    if not local:
        print("\n  没找到本地快照（--local 指定或放进 data/cbs_catalogue.tsv），跳过对比")
        return
    print("\n  本地快照：{}".format(local))
    loc, cov_note = load_local_snapshot(local, a.genus)
    if loc is None:
        return
    if cov_note:
        print("  ! 覆盖说明：{}".format(cov_note))
    live = {}
    for r in rows:
        k = norm_number(r.get("cbs_number", ""))
        if k:
            live[k] = r

    # 防线：本地侧一条都没有、线上却有一大堆 —— 这绝不可能是真实差异，
    # 而是本地快照没能按这个属过滤成功。宁可不出表，也不产出满屏假「新增」。
    if not loc and live:
        print("  本地侧一个号都没有，线上却有 {:,} 个 —— 这不是真实差异。".format(len(live)))
        print("  说明本地快照未能按属 {} 过滤成功。为免产出全是假「新增」的对比表，已跳过对比。".format(
            a.genus or "(全量)"))
        print("  你可以改用 --local 指定一个带正名的本地文件。")
        return

    only_live = sorted(set(live) - set(loc))
    only_local = sorted(set(loc) - set(live))
    both = set(live) & set(loc)
    st_chg = []
    for k in sorted(both):
        a1 = (loc[k].get("status") or "").strip()
        b1 = (live[k].get("Status of the strain") or "").strip()
        if a1 != b1:
            st_chg.append((live[k].get("cbs_number", k), a1, b1))

    print("\n" + "=" * 90)
    print("与本地快照的差异")
    print("=" * 90)
    print("  本地 {:,} 个号   线上 {:,} 个号   共有 {:,} 个".format(
        len(loc), len(live), len(both)))
    print("  线上有、本地没有（新增）：{:,} 个".format(len(only_live)))
    print("  本地有、线上没有（消失）：{:,} 个".format(len(only_local)))
    print("  模式状态有变化：{:,} 个".format(len(st_chg)))
    for k in only_live[:8]:
        print("     [新增] {}".format(live[k].get("cbs_number")))
    for k in only_local[:8]:
        print("     [消失] {}".format(loc[k].get("cbs_number")))
    for nm, o, n in st_chg[:8]:
        print("     [状态] {:<14} {!r} → {!r}".format(nm, o[:30], n[:30]))

    if a.genus:
        diff_tsv = os.path.join(a.out_dir, "cbs_{}_diff_{}.tsv".format(
            re.sub(r"[^A-Za-z0-9]", "", a.genus), stamp))
    else:
        diff_tsv = os.path.join(a.out_dir, "cbs_live_diff_{}.tsv".format(stamp))
    diff_tsv = os.path.abspath(diff_tsv)
    with io.open(diff_tsv, "w", encoding="utf-8", newline="") as f:
        f.write("kind\tcbs_number\told_value\tnew_value\n")
        for k in only_live:
            f.write("新增\t{}\t\t\n".format(live[k].get("cbs_number")))
        for k in only_local:
            f.write("消失\t{}\t\t\n".format(loc[k].get("cbs_number")))
        for nm, o, n in st_chg:
            f.write("状态变化\t{}\t{}\t{}\n".format(
                nm, o.replace("\t", " "), n.replace("\t", " ")))
    print("  已写 {}".format(diff_tsv))


if __name__ == "__main__":
    main()
