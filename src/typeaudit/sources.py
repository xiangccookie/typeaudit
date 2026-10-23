# -*- coding: utf-8 -*-
"""外部源 adapter 层：把「保藏目录 / 登记表 / 异名库 / NCBI」统一成同一接口。

为什么要这一层：
  J 规则库里 13 条需要外部源。若每条规则各自去读文件、各自解析格式，
  格式一变就要改 N 处。故把「取数」与「判据」分开——
  **源只管把它的原生格式变成统一字段；规则只面对统一字段。**

统一接口（最小集，够 J-04/J-05/J-09/J-16/J-19 用）：

    status_of(number)        -> {'rank': 'Neotype', 'name': 'Colletotrichum falcatum'} 或 None
    type_numbers_of(name)    -> ['CBS 127945', ...]      按学名反查模式号
    institution_type(code)   -> 'culture' | 'herbarium' | None
    current_name(name)       -> 接受名（异名关系）

源的三种形态与本项目实测：
  offline  已有本地文件，直接读（CBS 目录 / Biocollections / MycoBank）
  api      有公开接口（CGMCC `https://www.cgmcc.net/cgmccapi`）
  scrape   只能抓网页（BPI / MyCoPortal、ICMP）

⚠️ 数据来源与许可：本模块**不含任何数据**，只读用户自备的文件。
   使用者须自行确认各来源的使用条款（CBS 目录为 Westerdijk Institute 公开查询页）。
"""
import hashlib
import io
import time
import json
import os
import re

# ---------- 统一的状态串解析 ----------
RANK_PAT = r"(Ex-holotype|Ex-type|Ex-neotype|Holotype|Neotype|Isotype|Epitype|Paratype|Syntype|Type)"
RE_STATUS = re.compile(r"^\s*" + RANK_PAT + r"\s+of\s+(.+?)\s*$", re.I)


def parse_status(s):
    """把目录的 status 原文解析成 (rank, name)。解析不了返回 None。

    实测原文样例（本地 CBS 数据）：
      'Neotype of Colletotrichum falcatum'
      'Holotype of Chrysosporium magnasporum'
      'holotype of Colletotrichum hymenocallidis'      （小写也行）
      'Ex-holotype of Colletotrichum incanum'
    """
    if not s:
        return None
    m = RE_STATUS.match(s.strip())
    if not m:
        return None
    rank = m.group(1)
    rank = {"ex-holotype": "Ex-holotype", "ex-type": "Ex-type", "ex-neotype": "Ex-neotype"}.get(
        rank.lower(), rank.capitalize())
    return {"rank": rank, "name": m.group(2).strip()}


def norm_key(code, number):
    """
    统一键：'CBS' + '147945' -> 'CBS147945'。**保留小数点。**

    【为什么必须保留点 —— 这是对上一版的一次纠正】
        旧版连点一起删（'CBS 100.06' → 'CBS10006'），理由是「CBS 目录的 norm_key 列
        本来就是去点的」。但那是**照着有损的存量列对齐**，代价是：
            'CBS 100.06'（100 号菌株的第 06 支）与 'CBS 10006'（10006 号菌株）
            压成同一个键。
        本机实测：目录 62,761 行只有 60,052 个不同 norm_key ⇒ 2,709 行在撞车里被覆盖。
        而 CBSCatalogue 用 bykey[norm_key] 建字典、撞车时后写覆盖先写 ⇒
        **查 A 会静默返回 B 的目录记录**。这两者是不同的菌：
            CBS 100.06 = Absidia glauca ／ CBS 10006 = Pseudozyma thailandica。
        正确做法不是继续迁就有损的存量列，而是**从原始 cbs_number 重算键**
        （见 CBSCatalogue.__init__），这样两边才用同一把无损的尺子。
    """
    # 剥掉整个国家标签 <THA>（只删尖括号会留下 "THA"，导致同一号在两条调用路径下
    # 得到不同的键：norm_bare('MFLU<THA>','090228')='MFLU090228' 而 norm_key 得到
    # 'MFLUTHA090228' —— 本项目实测抓到的不一致）
    # 去空白、去国家标签、**去短横**，但【保留点】：
    #   点有语义 —— CBS 100.06（100 号第 06 支）≠ CBS 10006（10006 号），删点会撞车；
    #   短横没有语义 —— `LLTA- 01` 与 `LLTA 01` 是同一个号，保留短横会去重失效
    #   （本机实测：保留短横后三列表里同一个号占了两格）。
    t = re.sub(r"<[^>]*>", "", str(code) + str(number))
    return re.sub(r"[\s\-]", "", t).upper()


class Source(object):
    """所有源的基类。子类只需实现自己有的方法。"""
    id = "base"
    kind = "offline"
    license_note = ""

    def status_of(self, code, number):
        return None

    def code_map(self, code, number):
        """把【别家机构的号】映射回本机构的号。**没有这个能力就返回 None。**

        返回值三态，必须分清（本项目踩过坑）：
            None              = 我这个源【根本没有】跨库映射能力
            {"matches": []}   = 我有能力，但【查了没有】这条映射
            {"matches": [...]} = 查到了，附 via（映射依据的字段名）

        ⚠️ 为什么 None 与 [] 必须分开：基类 `type_numbers_of` 曾经 `return []`，
           于是「这个源没实现」在调用方看来与「查过确实没有」一模一样。
           本项目因此误判过两次（先说 ICMP 不记跨机构号，又说 MycoBank 没有）。
           **没有能力 ≠ 查无。** 新接口一律用 None 表达前者。
        """
        return None

    def type_numbers_of(self, name):
        return []

    def institution_type(self, code):
        return None

    def current_name(self, name):
        return None


# ---------------- 1. CBS / Westerdijk 保藏目录 ----------------
class CBSCatalogue(Source):
    """CBS 保藏目录（Westerdijk Institute）。

    kind = offline：读一份已抓取的 TSV（列：cbs_number / norm_key / status / has_status）。
    抓取 URL 模式：https://wi.knaw.nl/page/fungal_display/name/CBS%20<NUM>
    ⚠️ 本模块不代抓——抓取频率与条款须由使用者自行遵守。
    """
    id = "cbs"
    kind = "offline"
    license_note = "Westerdijk Institute 公开查询页；使用前请确认其条款"

    # 号码形态：'CBS 100.06' / 'ICMP:18578' / 'MFLU<THA>:090228'
    RE_NUM = re.compile(r"^\s*([A-Za-z][A-Za-z\.\-]{0,12}(?:<[^>]*>)?)\s*[:#]?\s*([\w\.\-]+)\s*$")

    def __init__(self, path):
        self.rows, self.bykey, self.collisions = [], {}, {}
        self.n_parsed = self.n_skipped = 0
        if path and os.path.exists(path):
            with io.open(path, encoding="utf-8-sig", newline="", errors="replace") as f:
                import csv as _csv
                self.rows = list(_csv.DictReader(f, delimiter="\t"))
            for r in self.rows:
                # ⚠️ 键从【原始 cbs_number】重算，不用存量的 norm_key 列 ——
                #    那一列是去点的（有损），正是撞车的来源。这一步同时丢掉非号码行
                #    （实测目录里混着作者名/注释片段，如 ' knoxdaviesiana'）。
                m = self.RE_NUM.match((r.get("cbs_number") or "").strip())
                if not m:
                    self.n_skipped += 1
                    continue
                self.n_parsed += 1
                k = norm_key(m.group(1), m.group(2))
                if k in self.bykey:
                    # 撞车：**绝不静默覆盖**（旧版就是覆盖，于是查 A 返回 B 的记录）。
                    # 两边都挪进 collisions，查询时判为「无从核验」。
                    self.collisions.setdefault(k, [self.bykey.pop(k)])
                    self.collisions[k].append(r)
                else:
                    self.bykey[k] = r

    def status_of(self, code, number):
        k = norm_key(code, number)
        if k in self.collisions:
            # 键有歧义 ⇒ 无从核验。宁可标「未核」，也不给出可能张冠李戴的记录。
            return {"rank": None, "name": None, "raw": "", "unannotated": True,
                    "ambiguous": True,
                    "collision": [(x.get("cbs_number") or "") for x in self.collisions[k]]}
        r = self.bykey.get(k)
        if not r:
            return None
        parsed = parse_status(r.get("status"))
        if parsed:
            return parsed
        # 没有模式地位标注时，**名字仍然要从这一行取**。
        # ⚠️ 原来这里直接把 name 置 None，于是 76% 没标注的记录「连名字都没有」，
        #    异名/并名判定根本跑不起来（用户实测发现：CBS 149616 目录明明记着
        #    Colletotrichum acericola，详情里却显示「目录记的名 —」）。
        #    名字与模式地位是两件事：没有地位标注 ≠ 没有名字。
        nm = (r.get("taxon_name") or "").strip()
        if not nm:
            tf = (r.get("taxon_full") or "").strip()
            nm = tf.split("|")[0].strip() if tf else ""
        return {"rank": None, "name": nm or None,
                "raw": r.get("status") or "", "unannotated": True}

    def type_numbers_of(self, name):
        out = []
        pat = re.compile(r"\b" + re.escape(name) + r"\b", re.I)
        for r in self.rows:
            if r.get("status") and pat.search(r["status"]):
                out.append(r.get("cbs_number") or "")
        return out

    def n_rows(self):
        # 报「可查询的条目数」= 无歧义键 + 撞车键下的所有记录（已跳过非号码行）
        return len(self.bykey) + sum(len(v) for v in self.collisions.values())


# ---------------- 1b. 通用机构目录（任意机构，放一份表就能核） ----------------
class NumberCatalogue(Source):
    """
    通用机构目录适配器 —— 让「无目录」变成「可核验」。

    【为什么需要它】
        本工具原来只有 CBS 一份机读目录，于是 ICMP / CMM / MFLU / IMI / BRIP 这些机构的号
        一律标「无目录」。而「无目录」的含义是「**无从核验**」——既不是「号不存在」，
        也不是「号有问题」。这是**适配器缺口**，不是数据缺陷。
        本机实测：ICMP 18578 / ICMP 19118 都是真号（机构在登记表里），却只能标无目录。

    【为什么不内置爬虫】
        各机构查询库的页面结构互不相同，硬写爬虫既脆弱又牵扯对方条款。
        所以这里只做**装载层**：你把导出的表放进 data/，工具立刻就能核。
        目录一旦就位，该机构还能多判一档 —— **号不存在**（原来判不了，只能说无目录）。

    文件约定（放 data/ 下）
        data/catalogue_<机构代码大写>.tsv，两列：
            number <TAB> status
        number 不带机构前缀（写 '18578'，不写 'ICMP 18578'）
        status 是模式地位原文，可为空
    """
    id = "numcat"
    kind = "offline"
    license_note = "由使用者自行导出并遵守来源机构条款"

    def __init__(self, dirpath):
        self.bycode = {}
        self.n_files = 0
        if dirpath and os.path.isdir(dirpath):
            for fn in sorted(os.listdir(dirpath)):
                m = re.match(r"^catalogue_([A-Za-z][A-Za-z0-9]{0,11})\.tsv$", fn)
                if not m:
                    continue
                code = m.group(1).upper()
                d = {}
                try:
                    with io.open(os.path.join(dirpath, fn), encoding="utf-8-sig",
                                 newline="", errors="replace") as f:
                        head = f.readline().rstrip("\n").rstrip("\r").split("\t")
                        ci = {h.strip().lower(): i for i, h in enumerate(head)}
                        i_n = ci.get("number", 0)
                        i_s = ci.get("status", 1 if len(head) > 1 else 0)
                        for line in f:
                            p = line.rstrip("\n").rstrip("\r").split("\t")
                            if len(p) > max(i_n, i_s) and p[i_n].strip():
                                d[norm_key("", p[i_n])] = p[i_s]
                except Exception:
                    continue
                if d:
                    self.bycode[code] = d
                    self.n_files += 1

    def has_code(self, code):
        return (code or "").strip().upper() in self.bycode

    def status_of(self, code, number):
        d = self.bycode.get((code or "").strip().upper())
        if not d:
            return None
        k = norm_key("", number)
        if k not in d:
            return None
        return parse_status(d[k]) or {"rank": None, "name": None, "raw": d[k] or "",
                                      "unannotated": True}

    def n_rows(self):
        return sum(len(v) for v in self.bycode.values())


# ---------------- 2. NCBI Biocollections 登记表 ----------------
def _lev(a, b, cap):
    """带上限的编辑距离（超过 cap 就提前返回，省时间）。"""
    if a == b:
        return 0
    if abs(len(a) - len(b)) > cap:
        return cap + 1
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i]
        row_min = i
        for j, cb in enumerate(b, 1):
            v = min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (ca != cb))
            cur.append(v)
            if v < row_min:
                row_min = v
        if row_min > cap:
            return cap + 1
        prev = cur
    return prev[-1]


class Biocollections(Source):
    """NCBI Biocollections（原 GrBio）：机构代码 -> 机构类型（培养物保藏 / 标本馆）。

    kind = offline：读一份 CSV（列含 code / name / type 之类；列名按实际表头自适配）。
    用途：判定规则 ③（标本/活培养物配对）依赖「机构类型」，该输入错了整批号会走错分支。
    """
    id = "biocollections"
    kind = "offline"
    license_note = "NCBI Biocollections 登记表"

    def __init__(self, path):
        self.bycode = {}
        if path and os.path.exists(path):
            import csv as _csv
            with io.open(path, encoding="utf-8-sig", newline="", errors="replace") as f:
                # ⚠️ 必须自适应分隔符：本适配器原来写死逗号，
                #    而使用者很可能给制表符分隔的表（实测把 TSV 喂进来时，
                #    整份文件被当成一行、code 列等于整行，于是**所有机构都查不到**）。
                head = f.readline()
                f.seek(0)
                delim = "\t" if head.count("\t") > head.count(",") else ","
                rd = _csv.DictReader(f, delimiter=delim)
                cols = {c.lower().strip(): c for c in (rd.fieldnames or [])}
                c_code = next((cols[k] for k in cols if "code" in k or k in ("abbreviation", "acronym")), None)
                c_kind = next((cols[k] for k in cols if "type" in k or "kind" in k or "category" in k), None)
                if c_code:
                    for r in rd:
                        code = (r.get(c_code) or "").strip().upper()
                        if code:
                            self.bycode[code] = (r.get(c_kind) or "").strip().lower() if c_kind else ""

    def institution_type(self, code):
        v = self.bycode.get((code or "").strip().upper())
        if v is None:
            return None
        if any(k in v for k in ("culture", "collection of microorganisms", "yeast")):
            return "culture"
        if any(k in v for k in ("herbarium", "fungarium", "museum", "botanic")):
            return "herbarium"
        return v or None

    def nearest_codes(self, code, max_dist=2, min_len=4, limit=5, only=None):
        """
        在登记表里找「差一两个字母」的候选代码 —— 用来区分「代码写错」和「根本不是代码」。

        【动机（本机实测）】
            `未知机构代码` 这个标记把两种完全不同的东西混在一起：
              ① 代码写错了 —— 如 CGMG，登记表里有 CGMCC（编辑距离 2）⇒ 真错值候选
              ② 根本不是机构代码 —— 如 NK 24 / BPDI I2（实验室内部编号）、
                 Paphi 2-1（疑似兰花属名 Paphiopedilum 被当成代码）⇒ 纯假阳性
            两组都查不到，光凭「不在登记表」分不开，必须再看「有没有近似代码」。

        min_len=4：代码太短时近似匹配毫无意义 —— 实测 NK（2 个字母）在编辑距离 1 上
        有 AK/ANK/BK/GK 一大堆邻居，全是噪声。

        only=：只在给定代码集合里找（如「有可查目录的机构」）。
        **必须过滤后再排序**，否则想要的候选会被 limit 截掉。
        """
        c = (code or "").strip().upper()
        if not c:
            return []
        out = []
        for k in self.bycode:
            # ⚠️ only= 必须在【排序之前】过滤。否则真正想要的那个代码
            #    （如 CGMCC，编辑距离 2）会按字母序排在别的距离 2 代码之后，
            #    被 limit 截掉 —— 本机实测 CGMG→CGMCC 就是这样丢的。
            if only is not None and k not in only:
                continue
            if len(k) < min_len or abs(len(k) - len(c)) > max_dist:
                continue
            d = _lev(c, k, max_dist)
            if d <= max_dist:
                out.append((d, k))
        out.sort()
        return [k for _, k in out[:limit]]

    def n_rows(self):
        return len(self.bycode)


# ---------------- 3. MycoBank 异名库 ----------------
class MycoBank(Source):
    """MycoBank：学名 -> 接受名（异名关系）。

    kind = offline：读一份 TSV（列含 name / current name 之类，按表头自适配）。
    ⚠️ 实测：**MycoBank 本身也会并名**（`C. hymenocallidis` 的 Current name 是
       `C. siamense`），所以它给出的「接受名」同样属于转引层判断，不能当终审。
    """
    id = "mycobank"
    kind = "offline"
    license_note = "MycoBank（Westerdijk Institute + CBS-KNAW）"

    def __init__(self, path):
        self.byname = {}
        self.odd_status = {}      # 非 Legitimate 的状态（拼写变体 / 无效 / 不合法）
        self.species_cur = {}     # 两词归并：种名 → {接受名…}（供并名自洽性检查）
        self.dup_conflict = {}    # 同名多行且接受名互相矛盾 → {名: [(状态, 接受名), …]}
        self._rows = {}           # 临时：名 → [(状态, 接受名), …]
        if path and os.path.exists(path):
            import csv as _csv
            with io.open(path, encoding="utf-8-sig", newline="", errors="replace") as f:
                rd = _csv.DictReader(f, delimiter="\t")
                cols = {c.lower().strip(): c for c in (rd.fieldnames or [])}
                c_n = next((cols[k] for k in cols if k in ("name", "taxon name", "epithet") or "name" == k), None)
                c_c = next((cols[k] for k in cols if "current" in k or "accepted" in k or "synonym" in k), None)
                c_s = next((cols[k] for k in cols if "status" in k), None)
                if c_n and c_c:
                    for r in rd:
                        n = (r.get(c_n) or "").strip()
                        if c_s:
                            st = (r.get(c_s) or "").strip()
                            if st and st.lower() != "legitimate" and n:
                                self.odd_status[n.lower()] = st
                        c = (r.get(c_c) or "").strip()
                        if n:
                            self._rows.setdefault(n.lower().strip(), []).append(
                                ((r.get(c_s) or "").strip() if c_s else "", c))
                        if n and c and n.lower() != c.lower():
                            self.byname[n.lower()] = c
                            # 两词（属+种）归并后的「名 → 接受名集合」。
                            # ⚠️ 必须在适配器里统一做：让各个检查项自己去截断，
                            #    会因遍历顺序/丢失规则不同而算出不同的数
                            #    （本机实测同一份数据得出 13 与 54 两个互相并名计数）。
                            a2 = " ".join(n.lower().split()[:2])
                            b2 = " ".join(c.lower().split()[:2])
                            if a2 and b2 and a2 != b2:
                                self.species_cur.setdefault(a2, set()).add(b2)
        # ── 同名重复条目消歧（必须在 __init__ 里做）──────────────────
        # ⚠️ 实测 MycoBank 同一个名字可能有多行，且**接受名互相矛盾**：
        #    `Diaporthe biguttulata` 一行 Illegitimate→pseudobiguttulata、
        #    另一行 Legitimate→eres；`chinensis`/`longispora`/`salinicola` 同理。
        #    旧逻辑 byname 取【最后一行】、odd_status 取【第一行非 Legitimate】，
        #    名字与状态来自不同的行 —— 两者不可能同时成立。
        #    现在：优先取 Legitimate 行的接受名；矛盾对子记入 dup_conflict 供上层标"待消歧"。
        for _k, _vs in self._rows.items():
            _curs = {v[1] for v in _vs if v[1]}
            if len(_curs) > 1:
                self.dup_conflict[_k] = _vs
                # ⚠️ 不能"优先取 Legitimate 行"：同名的两个条目可能是【同名异物】
                #    （不同作者、不同年份的同一个双名），Legitimate 那行很可能是
                #    另一个作者描述的另一物种。
                #    实测 Diaporthe biguttulata：NCBI 的 taxon 是 F. Huang 2015
                #    （模式 ZJUD 47），而 MycoBank 里 Legitimate 的那行是 Y.H. Gao &
                #    L. Cai 2014（后被并入 eres）—— 取 Legitimate 行会张冠李戴。
                #    正确做法：条目间作者无法区分时不猜，删掉映射并留给上层标"待消歧"。
                self.byname.pop(_k, None)
        self._rows = {}

    def current_name(self, name):
        return self.byname.get((name or "").strip().lower())

    def status_of(self, name):
        """名字在 MycoBank 的 Name status。

        只保存【非 Legitimate】的（约几百条），不存 42 万行全部 —— 省内存。
        为什么要它：`Colletotrichum karstii` 的状态是 Orthographic variant，
        它的接受名就是 `karsti`。这种情形属于【拼写变体】，不是"两个种都有效"。
        """
        return self.odd_status.get((name or "").strip().lower(), "")

    def n_rows(self):
        return len(self.byname)


# ---------------- 注册与装载 ----------------


# 各适配器「本该在哪」的已知候选位置 —— 配置里路径不存在时按序回退
# （用户实测：config/sources.json 写的是 data/cbs_catalogue.tsv，而真实文件在别处，
#   于是所有 CBS 号都变成「无目录」/「未核」，真错值一条都出不来）
_FALLBACK_PATHS = {
    "CBSCatalogue": ["data/cbs_catalogue.tsv",
                     r"C:\Users\Administrator\Desktop\DeepSeek\bio原理\out\cbs_catalogue.tsv"],
    "MycoBank": ["data/mycobank_full.tsv",
# ⚠️ 2026-10 目录搬迁：Tool 从 Desktop\Tool 搬到 Desktop\sol\Tool，
#    旧路径失效会让 MycoBank 静默空载 ⇒ 异名判定退化成「目录记为别种」的假阳性。
                 r"C:\Users\Administrator\Desktop\sol\Tool\脚本\调取mycobank数据，分析sp_curret,编码，是否为复合体\mycobank_full.tsv",
                 r"C:\Users\Administrator\Desktop\Tool\脚本\调取mycobank数据，分析sp_curret,编码，是否为复合体\mycobank_full.tsv"],
    "Biocollections": ["data/biocollections_codes.csv"],
}


# ---------------------------------------------------------------- 缓存根
# 仓库根：本文件在 <REPO>/src/typeaudit/sources.py，往上三层
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
CACHE_ROOT = os.path.join(REPO_ROOT, ".cache")

# 这四个源的 path 语义是【缓存目录】，不是数据文件。
# ⚠️ 它们必须有一个与 CWD 无关的确定位置。原来这四个在 _FALLBACK_PATHS 里
#    **根本没有条目**，path 为空时缓存目录退化成 CWD 相对路径 —— 本机实测
#    把 270.88 MiB 的 GenBank 缓存写进了 config/ 配置目录，另有
#    out_audit/out_audit/.cache 这种双层目录（相对路径被拼了两次）。
# CBSMap 的联网查询缓存目录（原实现查完 CBS 写在这里）
_CBSMAP_CACHE = os.path.join(REPO_ROOT, ".cache", "cbsmap")

_CACHE_ADAPTERS = {"GenBankAPI", "TaxonomyAPI", "CGMCCAPI", "ICMPAPI"}

_ORIGIN_RE = None


def _strip_origin(txt):
    """剥掉 GenBank 记录里的 ORIGIN 碱基序列块，只留 LOCUS/FEATURES。

    为什么（本机实测）：本项目只读 FEATURES 里的限定符
    （/culture_collection、/strain、/type_material、/note），**碱基序列一个字都不用**。
    但带基因组记录的种，ORIGIN 占文件 74.6% ——
    cc_Colletotrichum_abscissum_0.gb 90,859,971 B 里 67,822,896 B 是碱基。
    按 418 个种外推，全文落盘要 35.37 GiB，剥掉后 8.97 GiB。

    解析口径不变：所有 regex 只匹配限定符与 XML 标签，ORIGIN 里没有这些。
    末尾的 `//` 记录终止符**保留**（只删 ORIGIN 到 `//` 之间的内容），
    因为可能有代码按 `//` 切记录。
    """
    global _ORIGIN_RE
    if _ORIGIN_RE is None:
        import re as _r
        _ORIGIN_RE = _r.compile(r"^ORIGIN\b.*?(?=^//)", _r.S | _r.M)
    return _ORIGIN_RE.sub("", txt)


def _find_fallback(adapter, base):
    """配置路径不存在时，按候选表 + find_refdata 的探测结果找一份能用的。"""
    for rel in _FALLBACK_PATHS.get(adapter, []):
        p = rel if os.path.isabs(rel) else os.path.join(base, rel)
        if os.path.exists(p):
            return p
    # 再问 find_refdata 写过没有
    try:
        probe = os.path.expanduser(r"~\.dsh\typeaudit_refs.json")
        if os.path.exists(probe):
            j = json.load(io.open(probe, encoding="utf-8"))
            for k, v in (j.get("paths") or j or {}).items():
                # 探测结果可能是两种形状：{key: "路径"} 或 {key: {"path": "路径", ...}}
                # —— 只认字符串会把嵌套那种整条跳过（本机实测：临时配置下
                #    biocollections 因此空载，行数=0，该源不参与判定）。
                if isinstance(v, dict):
                    v = v.get("path") or ""
                if isinstance(v, str) and os.path.exists(v):
                    if adapter == "MycoBank" and "mycobank" in k.lower():
                        return v
                    if adapter == "CBSCatalogue" and ("cbs" in k.lower()):
                        return v
                    if adapter == "Biocollections" and ("biocol" in k.lower() or "grbio" in k.lower()):
                        return v
    except Exception:
        pass
    return None


def attach_memo(sources, table, only=None):
    """把 NCBI 结果表挂到「联网型」源上（只挂缓存型适配器）。

    为什么不挂在 CBS/MycoBank 这些本地表上：它们读的是本地 tsv，
    本来就不联网，记忆层只会白占内存。
    """
    from . import ncbitable as _nt
    n = 0
    for name, obj in (sources or {}).items():
        if obj is None:
            continue
        if only is not None and name not in only:
            continue
        if type(obj).__name__ not in _CACHE_ADAPTERS:
            continue
        try:
            n += _nt.attach(obj, table)
        except Exception:
            pass
    return n


def load(config_path):
    """按 config/sources.json 装载所有可用源。缺文件即跳过（不报错、不臆造）。"""
    cfg = json.load(io.open(config_path, encoding="utf-8")) if os.path.exists(config_path) else {}
    base = os.path.dirname(os.path.abspath(config_path))
    out = {}
    for name, spec in (cfg.get("sources") or {}).items():
        cls = KINDS.get(spec.get("adapter"))
        if not cls:
            continue
        p = spec.get("path") or ""
        # 缓存型适配器（联网源）：path 是【缓存目录】。给它一个与 CWD 无关的
        # 确定位置，然后直接构造，不走下面那套「数据文件」的两级相对路径解析
        # —— 那套解析会把相对路径按配置所在目录再拼一次，网页下就是双层目录。
        if spec.get("adapter") in _CACHE_ADAPTERS:
            if not p:
                p = CACHE_ROOT
            elif not os.path.isabs(p):
                p = os.path.normpath(os.path.join(REPO_ROOT, p))
            try:
                out[name] = cls(p)
            except Exception:
                pass
            continue
        # ⚠️ 空路径也必须走回退。原来只在「给了路径但文件不存在」时回退，
        #    于是配置里写 {"adapter":"MycoBank","path":""} 会走到 cls("")，
        #    **静默装载成空表** —— 而空表的后果不是「缺数据」，是「所有名字都查不到」，
        #    比报错危险得多（本机实测：网页版留空 MycoBank 时，三列表里的
        #    「目录记为别种」无法被 MycoBank 解释，于是本来非错的号被留在错值栏）。
        need = (not p) or (not os.path.exists(p if os.path.isabs(p) else os.path.join(base, p))
                           and not os.path.exists(os.path.join(os.path.dirname(base), p)))
        if need:
            alt = _find_fallback(spec.get("adapter"), os.path.dirname(base))
            if alt:
                p = alt
        if p and not os.path.isabs(p):
            # ⚠️ 相对路径要**两级**找：config 目录 → 其父目录（= 仓库根）。
            #    原来只按 config 目录解析，于是配置里写 `data/xxx` 会被解成
            #    `config/data/xxx`，静默加载成**空表**——而空表的后果不是"缺数据"，
            #    是"所有号都被判成查无"（比报错危险得多）。
            for b in (base, os.path.dirname(base)):
                cand = os.path.join(b, p)
                if os.path.exists(cand):
                    p = cand
                    break
            else:
                p = os.path.join(base, p)
        try:
            out[name] = cls(p)
        except Exception:
            pass
    # 通用机构目录：扫 data/catalogue_<代码>.tsv，有就注册为 numcat 源。
    # 放一份表 = 该机构的号从「无目录」升级为「可核验」，还能多判一档「号不存在」。
    try:
        for d in (os.path.join(os.path.dirname(base), "data"), os.path.join(base, "data")):
            nc = NumberCatalogue(d)
            if nc.n_files:
                if out.get("numcat") is not None:
                    out["numcat"].bycode.update(nc.bycode)
                    out["numcat"].n_files += nc.n_files
                else:
                    out["numcat"] = nc
    except Exception:
        pass
    return out

# ---------------- 4. NCBI Entrez：GenBank / Taxonomy ----------------
# ---------------------------------------------------------------- 网络兜底超时
_NET_TIMEOUT_SET = [False]


def _net_timeout():
    """给所有网络操作设一个兜底超时（秒）。可用 TYPEAUDIT_NET_TIMEOUT 覆盖，设 0 关闭。

    ⚠️ 为什么必须有（本机实测，2026-10-16）：
       **Bio.Entrez 自身一个 timeout 都不设** —— 在安装的 Bio/Entrez/__init__.py 里
       `timeout=` 出现 **0 次**。NCBI 偶尔会「接受 TCP 连接但不回数据」，
       此时 urllib 会**永远挂着**，不会抛异常、不会超时、也不会写缓存。
       实测症状：流水线停在 `[1/20] Colletotrichum abscissum`，
       16.8 分钟零缓存写入；netstat 显示一条到
       eutils.ncbi.nlm.nih.gov(34.107.134.59):443 的 ESTABLISHED 连接。
       而且此时「停止」按钮**无效** —— 它只在【种与种之间】生效
       （collect() 内部的网络调用无法硬中断）。

    为什么用 setdefaulttimeout 而不是逐个改调用：
      · 它是**每次 socket 操作**的超时，不是整个下载的时间上限 ——
        只要数据还在持续到达就不会触发，所以大文件（实测 86 MB 的 GenBank）安全；
      · 只影响**之后新建**的 socket。本地网页服务的监听 socket 在启动时就已建好，
        accept() 继承的是监听 socket 自己的 timeout，不受这个默认值影响。

    超时后的行为：该次查询失败 → 被 collect() 的 except 记进「notes」→
    该种照常出表、只是这一项标注失败原因。**比无声挂死好得多。**
    """
    if _NET_TIMEOUT_SET[0]:
        return
    _NET_TIMEOUT_SET[0] = True
    try:
        import socket as _s
        raw = (os.environ.get("TYPEAUDIT_NET_TIMEOUT") or "60").strip()
        sec = float(raw)
    except Exception:
        return
    if sec > 0:
        try:
            _s.setdefaulttimeout(sec)
        except Exception:
            pass


def _entrez(email=None, tool="typeaudit"):
    """惰性导入 biopython；要求提供邮箱（NCBI 规定）。"""
    _net_timeout()      # ← 必须在发起任何 Entrez 请求前设兜底超时，见上方说明
    import os as _os
    try:
        from Bio import Entrez
    except ImportError:
        raise SystemExit("需要 biopython：pip install biopython")
    Entrez.email = email or os.environ.get("NCBI_EMAIL")
    if not Entrez.email:
        raise SystemExit("NCBI 要求邮箱：设 NCBI_EMAIL 或用 --email。")
    Entrez.tool = tool
    return Entrez


def _rd(h):
    d = h.read()
    return d.decode("utf-8", "replace") if isinstance(d, bytes) else d


# ---------------------------------------------------------------- 查询日志
# 记录每一次 esearch 的 (库, 检索式, 窗口, 服务端命中数)。
# 为什么需要：`号不存在` 在本工具里其实是若干合取项同时成立，其中一项是
# 「我的窗口够不够大」。本机实测出现过「检索命中 1,917 条、我们只取前 40 条」
# 导致类型材料召回 0% —— 而输出里只有「号 + 标记」，看不出窗口问题。
# 把窗口和命中数记下来，这类漏检才可能被一眼看穿。
QUERY_LOG = []


def _log_query(db, term, retmax, count):
    QUERY_LOG.append({"db": db, "term": term, "retmax": retmax,
                      "count": count})
    if len(QUERY_LOG) > 20000:          # 防止长跑时无限增长
        del QUERY_LOG[:5000]


class _EntrezBase(Source):
    kind = "api"

    def __init__(self, cache_dir=".cache", email=None, pause=1.3):
        self.cache_dir = cache_dir
        self.email = email
        self.pause = pause
        if cache_dir:
            os.makedirs(cache_dir, exist_ok=True)

    def _cached(self, name, make, alt=()):
        """按 (接口, 参数) 落盘缓存——重跑零请求，这是可复现性的前提。

        `alt` = 旧缓存文件名兜底表。当缓存键换了（例如本次的「按请求内容寻址」），
        老文件仍然有效时先用着 —— 否则一次键变更会把上百 MB 已下载的数据全部作废、
        逼着重新下（本机实测：光 abscissum 一种就是 90.9 MB）。
        命中 alt 后会把内容复制到新键下，下次直接命中新键。
        """
        if not self.cache_dir:
            return make()
        p = os.path.join(self.cache_dir, name)
        if os.path.exists(p):
            return io.open(p, encoding="utf-8", errors="replace").read()
        for a in alt:
            q = os.path.join(self.cache_dir, a)
            if os.path.exists(q):
                s = io.open(q, encoding="utf-8", errors="replace").read()
                try:
                    io.open(p, "w", encoding="utf-8", errors="replace").write(s)
                except Exception:
                    pass
                return s
        s = make()
        # GenBank 文本落盘前剥掉 ORIGIN 序列块（只影响磁盘，不影响解析——
        # 全链路零处读序列，见 _strip_origin 的说明）。
        if name.endswith(".gb"):
            s = _strip_origin(s)
        io.open(p, "w", encoding="utf-8", errors="replace").write(s)
        return s

    def _read_stream(self, h, label=""):
        """流式读响应体，每 4 MB 打一行进度。

        ⚠️ 为什么需要（本机实测）：`_cached` 只在【整个响应读完之后】才写盘，
        所以一个 90.9 MB 的 efetch 在外部看来与「卡死」完全一样
        （零缓存写入 + 一条 ESTABLISHED 连接）。实测那次卡了 19 分钟，
        分不清是在下还是死了。有了进度行就能一眼分辨。
        """
        buf = []
        got = 0
        nxt = 4 * 1024 * 1024
        while True:
            b = h.read(65536)
            if not b:
                break
            if isinstance(b, bytes):
                got += len(b)
                buf.append(b.decode("utf-8", "replace"))
            else:
                got += len(b)
                buf.append(b)
            if got >= nxt:
                print("      [ncbi] 已收 {:.1f} MB …".format(got / 1048576.0))
                # ⚠️ 这里必须是 固定步长 递增。原来写成 `nxt += (nxt + 4MB)`
                #    → nxt 翻倍增长，进度行会越来越稀、最后几乎不报。
                nxt += 4 * 1024 * 1024
        return "".join(buf)

    # 单批的体量上限（碱基数）。超过就拆批，让大响应变成多个小响应。
    # 可用环境变量 TYPEAUDIT_MAX_BASES 覆盖；设 0 关闭拆批（保持旧行为）。
    def _max_bases(self):
        if getattr(self, "_mb", None) is None:
            try:
                self._mb = int(float(os.environ.get("TYPEAUDIT_MAX_BASES", "10000000")))
            except Exception:
                self._mb = 10000000
        return self._mb

    def _sizes(self, ids):
        """用 esummary 拿每条记录的碱基数（只用来决定怎么拆批，不参与判定）。"""
        if getattr(self, "_sz", None) is None:
            self._sz = {}
        missing = [i for i in ids if i not in self._sz]
        if missing:
            for k in range(0, len(missing), 300):
                part = missing[k:k + 300]
                key = "esz_{}_{}.xml".format(len(part), hashlib.sha1(
                    ",".join(part).encode()).hexdigest()[:20])
                def mk(c=part):
                    time.sleep(self.pause)
                    return _rd(_entrez(self.email).esummary(db="nucleotide",
                                                            id=",".join(c)))
                sx = self._cached(key, mk)
                for m in re.finditer(r"<DocSum>.*?</DocSum>", sx, re.S):
                    doc = m.group(0)
                    mid = re.search(r"<Id>(\d+)</Id>", doc)
                    # ⚠️ 字段名是 Length，**不是** Slen。
                    #    Slen 是蛋白 esummary 用的；nucleotide 用 Length。
                    #    写错不会报错，只会让每条记录的体量都是 0 ⇒ 拆批永不触发，
                    #    而拆批失效是静默的（表现为"大响应还是老样子"）。本机实测踩到。
                    sl = re.search(r'Name="Length"[^>]*>(\d+)<', doc)
                    if mid:
                        self._sz[mid.group(1)] = int(sl.group(1)) if sl else 0
                for i in part:                 # esummary 没回的按 0 处理（不倒扣）
                    self._sz.setdefault(i, 0)
        return {i: self._sz.get(i, 0) for i in ids}

    def _chunks_by_size(self, ids, batch):
        """先按 batch 切，再把「体量过大」的那批按碱基数继续拆。

        为什么不无脑用固定小批：拆批会改变 efetch 的 ID 列表 ⇒ 缓存键变化 ⇒
        已下载的几百 MB 全部作废。所以**只有确实过大的批才拆**，
        小记录的批保持原样（旧缓存继续命中）。
        """
        cap = self._max_bases()
        rough = [ids[k:k + batch] for k in range(0, len(ids), batch)]
        if not cap:
            return [(c, k) for k, c in enumerate(rough)]
        try:
            sz = self._sizes(ids)
        except Exception:
            return [(c, k) for k, c in enumerate(rough)]
        out = []
        for k, c in enumerate(rough):
            if sum(sz.get(i, 0) for i in c) <= cap or len(c) == 1:
                out.append((c, k))
                continue
            cur, acc = [], 0
            for i in c:
                b = sz.get(i, 0)
                if cur and acc + b > cap:
                    out.append((cur, k)); cur, acc = [], 0
                cur.append(i); acc += b
            if cur:
                out.append((cur, k))
        return out

    # ---------------------------------------------------------------- 别名索引
    # 为什么需要它：`gbx_<id列表哈希>.gb` 是「按请求内容寻址」的名字，而旧缓存用的是
    # 调用方起的名字（cc_/sv_/tm_/nts_）。两者内容相同时，**不能再复制一份** ——
    # 实测那样做会凭空多出 179 MB（abscissum 一份 90.9 MB、一份 88.1 MB 各存了两遍）。
    # 改用一张小索引表记「新键 → 已有文件名」，几 KB 就能替代上百 MB。
    def _idx_load(self):
        if getattr(self, "_gidx", None) is None:
            self._gidx = {}
            p = os.path.join(self.cache_dir or "", "gbx_index.json")
            if os.path.isfile(p):
                try:
                    self._gidx = json.load(io.open(p, encoding="utf-8")) or {}
                except Exception:
                    self._gidx = {}
        return self._gidx

    def _idx_put(self, key, fname):
        if not self.cache_dir:
            return
        idx = self._idx_load()
        if idx.get(key) == fname:
            return
        idx[key] = fname
        try:
            io.open(os.path.join(self.cache_dir, "gbx_index.json"), "w",
                    encoding="utf-8").write(
                        json.dumps(idx, ensure_ascii=False, indent=0))
        except Exception:
            pass

    def _seq_stop(self):
        """序列截断长度。0 = 关闭。可用环境变量 TYPEAUDIT_SEQ_STOP 覆盖。"""
        if getattr(self, "_ss", None) is None:
            try:
                self._ss = int(float(os.environ.get("TYPEAUDIT_SEQ_STOP", "10")))
            except Exception:
                self._ss = 10
        return self._ss

    def _fetch_gb_raw(self, ids):
        """取一批 GenBank 记录文本（唯一的 efetch 出口）。

        【为什么带 seq_start/seq_stop】
        它是在告诉 NCBI：「这条记录的序列，我只要第 1 到第 N 个碱基」。
        GenBank 平文件的格式是 LOCUS + FEATURES + ORIGIN + //，ORIGIN（碱基序列）
        是格式的一部分、NCBI 不给字段级选择 ⇒ 想要 FEATURES 就必须连 ORIGIN 一起收下。
        而这个项目只读 FEATURES 里的限定符，ORIGIN 一个字母都不用（收下后立刻剥掉）。
        实测（BK072762.1，声明 45,176 bp）：58,744 → 1,546 字节，省 97.4%。

        【代价，已逐条核对】
        NCBI 会把特征【位置列】和 LOCUS 行的长度改成 `1..10`：
              source   1..529   →  source   1..10
            LOCUS … 529 bp      →  LOCUS … 10 bp
        全仓库作用在 GenBank 文本上的 27 条正则里，**没有一条**读这两处；
        唯一用 `^LOCUS` 的地方（specimen.py:64）是 `len(re.findall(...))`，
        只数记录条数，与长度字段无关。四个关键字段的值实测完全相同。
        设 TYPEAUDIT_SEQ_STOP=0 可关闭（回到整袋取）。
        """
        n = self._seq_stop()
        kw = {"seq_start": 1, "seq_stop": n} if n > 0 else {}

        def one(extra=None):
            time.sleep(self.pause)
            k = dict(kw)
            if extra:
                k.update(extra)
            h = _entrez(self.email).efetch(db="nucleotide", id=",".join(ids),
                                           rettype="gb", retmode="text", **k)
            return self._read_stream(h)

        txt = one()
        # ⚠️ 回退：截断版万一没带回记录（连 LOCUS 都没有），就整袋重取一次。
        #    绝不能把空响应当结果返回 —— 那会被 _cached 当成"取到了"写进缓存，
        #    此后这个键永远命中一个空值，静默丢掉一整批记录。
        if n > 0 and "LOCUS" not in txt:
            print("      [ncbi] ⚠️ 截断取回无记录，整袋重取一次（{} 条）".format(len(ids)))
            time.sleep(self.pause)
            h = _entrez(self.email).efetch(db="nucleotide", id=",".join(ids),
                                           rettype="gb", retmode="text")
            txt = self._read_stream(h)
        return txt

    def _efetch_gb(self, ids, legacy=(), org="", batch=100):
        """按【请求内容】取一批 GenBank 记录。

        ⚠️ 为什么（本机实测）：`cc_` 与 `sv_` 用不同调用方名字缓存同一份数据 ——
        `cc_X_0.gb` 与 `sv_X_0.gb` 的 sha256 **完全相同**（falcatum/siamense 实测），
        因为两者 esearch 出的前 100 个 ID 逐项相同。于是同一批字节被下了两遍。
        对带全基因组记录的种，这一遍是 **90.9 MB**；abscissum 那次「停在 [1/20]
        19 分钟」就是它。

        缓存键 = 这批 ID 的哈希 ⇒ `cc_` / `sv_` / `tm_` 三个入口自然共享存储。
        **取哪些记录仍由 ID 列表决定**，所以口径不变。
        """
        key = "gbx_{}.gb".format(hashlib.sha1(",".join(ids).encode()).hexdigest()[:24])
        # ⚠️ 旧缓存必须在【拆批之前】查。顺序反了会造成最贵的浪费：
        #    abscissum 已有 cc_abscissum_0.gb（90.9 MB），若先拆批就会绕过它、
        #    把同一批数据按新粒度再下一遍。本机实测踩到过这个顺序问题。
        #    key 与 legacy 覆盖的是**同一个 ID 列表**，所以谁先命中都同样正确。
        if self.cache_dir:
            look = [key]
            _alias = self._idx_load().get(key)
            if _alias:
                look.append(_alias)
            look.extend(legacy)
            for nm in look:
                if not nm:
                    continue
                q0 = os.path.join(self.cache_dir, nm)
                if os.path.exists(q0):
                    # ⚠️ 直接读文件：写成 self._cached(key, lambda: "") 的话，
                    #    万一文件在两步之间消失，_cached 会把空字符串写进缓存，
                    #    此后这个键永远"命中"一个空值 —— 静默丢失一整批记录。
                    s0 = io.open(q0, encoding="utf-8", errors="replace").read()
                    # 旧缓存是加 _strip_origin 之前下的，带着大段碱基序列
                    # （实测 abscissum 一份 90.9 MB、一份 88.1 MB，序列占 74%）。
                    # 命中的时候顺手就地剥掉：磁盘立刻掉下来，解析结果不变。
                    if "\nORIGIN" in s0:
                        s1 = _strip_origin(s0)
                        if s1 != s0:
                            try:
                                io.open(q0, "w", encoding="utf-8",
                                        errors="replace").write(s1)
                                s0 = s1
                            except Exception:
                                pass
                    # 命中别名时只【记一笔索引】，不再复制整份文件
                    if nm != key:
                        self._idx_put(key, nm)
                    return s0
        chunks = self._chunks_by_size(ids, batch)
        if len(chunks) == 1 and chunks[0][0] == list(ids):
            def mk():
                return self._fetch_gb_raw(ids)
            return self._cached(key, mk, alt=legacy)
        parts = []
        for n, (c, kk) in enumerate(chunks, 1):
            ck = "gbx_{}.gb".format(hashlib.sha1(",".join(c).encode()).hexdigest()[:24])
            def mk(c=c):
                return self._fetch_gb_raw(c)
            print("      [ncbi] 拆批 {}/{}：{} 条".format(n, len(chunks), len(c)))
            parts.append(self._cached(ck, mk))
        joined = "".join(parts)
        if self.cache_dir:
            try:
                io.open(os.path.join(self.cache_dir, key), "w",
                        encoding="utf-8", errors="replace").write(joined)
            except Exception:
                pass
        return joined

    def _esearch(self, db, term, retmax=20, tag=""):
        # ⚠️ 缓存键**必须带 retmax**：`retmax=0`（只取计数）与 `retmax=20`（取 ID）
        #    返回的是不同的响应体。键若不区分，先跑的计数查询会把 ID 列表挤掉——
        #    表现为「count 正确、ids 为空」的静默失败（本项目踩过一次）。
        key = "es_{}_{}_{}.xml".format(db, retmax, re.sub(r"\W+", "_", term)[:80])
        def mk():
            import time as _t
            _t.sleep(self.pause)
            E = _entrez(self.email)
            h = E.esearch(db=db, term=term, retmax=retmax)
            return _rd(h)
        x = self._cached(key, mk)
        _m = re.search(r"<Count>(\d+)</Count>", x)
        _log_query(db, term, retmax, int(_m.group(1)) if _m else None)
        return x


class GenBankAPI(_EntrezBase):
    """GenBank / RefSeq（db=nucleotide）。

    提供：
      records_of_number(code, number)  该号在 GenBank 的记录与属分布 → J-03/J-18
      type_material_of(species)        该种记录的 /type_material 声明 → J-16
      retrievability(code, number)     逐字段检索命中数 → J-21
      note_type_status(accession)      /note 里的 Type status → J-20
    """
    id = "genbank"
    def __init__(self, path=None, email=None, pause=None):
        super(GenBankAPI, self).__init__(
            cache_dir=(path or ".cache"),
            email=(email or os.environ.get("NCBI_EMAIL") or None),
            pause=(float(pause) if pause else 1.3))

    @staticmethod
    def num_term(code, number):
        """把「机构代码 + 号」变成一条尽量准的检索式。

        实测四种写法的命中数（同一个号 ICMP 1783 / CBS 112996）：

            写法                     结果        性质
            "ICMP 1783"（带空格加引号）  20,561     ❌ 被拆成两个词分别匹配 → 假阳性爆炸
            "CBS:112996"（冒号加引号）   20,561     ❌ 同上
            CBS112996（连写）            0（CBS）    ⚠️ 对带点号 CBS 号召不回
            ICMP:1783（冒号不加引号）    53         ✅ 准
            CBS:221.81（冒号不加引号）   9          ✅ 准
            ATCC:56816 / STE-U:5292     16 / 4      ✅ 准

        结论：**冒号 + 不加引号**是唯一在各机构都稳的写法；连写形作补充分支
        （有些记录只写 ICMP1783 不带冒号）。
        """
        c = (code or "").strip()
        n = (number or "").strip()
        flat = re.sub(r"[^A-Za-z0-9]", "", n)
        # 只用冒号形。实测把连写形 OR 进来反而放大：
        # CBS22181 命中 180 条（真值 9）、ICMP1783 分支让总数回到 20,561。
        return "{}:{}".format(c, n)

    def _count(self, term, db="nucleotide"):
        """只取命中数（esearch retmax=0）。"""
        xml = self._esearch(db, term, retmax=0)
        m = re.search(r"<Count>(\d+)</Count>", xml)
        return int(m.group(1)) if m else 0

    def _records_for(self, term, retmax):
        """esearch 拿到 ID → 取 GenBank 记录文本。两条路径都已缓存。"""
        xml = self._esearch("nucleotide", term, retmax=retmax)
        ids = re.findall(r"<Id>(\d+)</Id>", xml)
        if not ids:
            return "", []
        return self._efetch_gb(ids), ids

    @staticmethod
    def _parse_records(body):
        """把 GenBank 平文件切成一条条记录（// 是记录分隔符）。"""
        return [r for r in (body or "").split("\n//\n") if r.strip()]

    def records_of_number(self, code, number, retmax=20):
        """该号在 GenBank 有多少条记录、涉及哪些属、标题是什么。

        返回 {"count","genera","titles","ids"}：
          count  = 命中总数（esearch retmax=0 的真实总数，不受 retmax 截断）
          genera = {属名: 条数}  → 用来判「这个号被别的属占用」（跨属错号）
        """
        # ⚠️ 检索式用【连写形】（CBS112996）。带空格的引号形会被 Entrez 拆成
        #    "CBS" 与 "112996" 两个词分别匹配 —— 实测命中 20,561 条无关记录。
        term = self.num_term(code, number)
        body, ids = self._records_for(term, retmax)
        # ⚠️ 必须逐条核对「这条记录真的写了这个号」。
        #    检索式 `ICMP:1783` 仍会命中同现两个词的其他记录（实测取到
        #    Piscirickettsia、Thunnus 这类非目标属）。不核对，属分布就是错的，
        #    而「跨属错号」正是靠属分布判的 —— 会直接把结论带偏。
        need = re.sub(r"[^A-Z0-9]", "", "{}{}".format(code, number).upper())
        genera, titles = {}, []
        for rec in self._parse_records(body):
            if need and need not in re.sub(r"[^A-Z0-9]", "", rec.upper()):
                continue
            org = re.search(r'/organism="([^"]+)"', rec)
            if org:
                g = org.group(1).split()[0]
                genera[g] = genera.get(g, 0) + 1
            d = re.search(r"(?ms)^DEFINITION\s+(.+?)(?=^ACCESSION|^VERSION|^KEYWORDS)", rec)
            if d:
                titles.append(re.sub(r"\s+", " ", d.group(1)).strip()[:160])
        try:
            n = self._count(term)
        except Exception:
            n = len(ids)
        return {"code": code, "number": number, "count": n,
                "n_fetched": len(ids), "genera": genera, "titles": titles, "ids": ids}

    def _qual_values(self, species, qual, retmax):
        """取某个种所有记录里指定限定符（/strain、/type_material）的值，去重保序。"""
        term = '"{}"[Organism]'.format(species)
        body, _ = self._records_for(term, retmax)
        seen, out = set(), []
        pat = re.compile(r'/%s="([^"]*)"' % re.escape(qual))
        for rec in self._parse_records(body):
            for v in pat.findall(rec):
                v = v.strip()
                if v and v not in seen:
                    seen.add(v)
                    out.append(v)
        return out

    def strain_values_of(self, species, retmax=200):
        """该种记录里的 /strain 值（论文内部编号，只作旁证，**不是**种名）。"""
        return self._qual_values(species, "strain", retmax)

    def culture_collections_of(self, species, retmax=300, extra=None):
        """该种记录里的 /culture_collection 值（保藏号）—— 第 1 列的主键。

        extra="type" 时只取「带 /type_material 的那些记录」里的保藏号：
        这是「类型材料记录声明的号」，与「序列库里随便一条记录提到的号」是两回事，
        不能混在一列里比。

        ⚠️ GenBank 里写作 CBS:112996（带冒号），本项目统一显示成 CBS 112996（空格），
        与其它列、与机构目录的写法一致。
        """
        term = '"{}"[Organism]'.format(species)
        # ⚠️ 扫描窗口要够大。实测两个反例：
        #    · acutatum 前 300 条记录里没有 CBS 112980（在第 300–1200 条之间）
        #    · aenigma  前 300 条记录里【一个保藏号都没有】，1200 条才有 2 个
        #    窗口开 300 会静默漏号（表看起来正常，只是少了几管菌）。
        retmax = max(int(retmax or 0), 1500)
        # 记下「总共有多少条 / 这次扫了多少」，供上层标「覆盖不全」。
        # 不这么做的话，大种（aenigma 18,346 条）会静默少号，读者看不出来。
        try:
            self.last_coverage = {"total": self._count(term), "scanned": int(retmax),
                                  "species": species}
        except Exception:
            self.last_coverage = None
        body, _ = self._records_for(term, retmax)
        seen, out = set(), []
        for rec in self._parse_records(body):
            if extra == "type" and '/type_material=' not in rec:
                continue
            for v in re.findall(r'/culture_collection="([^"]*)"', rec):
                v = v.strip().replace(":", " ")
                v = re.sub(r"\s+", " ", v)
                if v and v not in seen:
                    seen.add(v)
                    out.append(v)
        return out

    def type_material_cc_of(self, species, retmax=1500):
        """只在【类型材料记录】里取 /culture_collection —— 第 3 列的号码来源。

        为什么要单开一路：上面 culture_collections_of(extra="type") 只扫 retmax 条，
        而类型材料声明常藏在很靠后的记录里（实测某些种要扫到第几百条），
        所以这条给更大的扫描窗口。
        """
        return self.culture_collections_of(species, retmax=retmax, extra="type")

    def type_material_of(self, species, retmax=1500):
        """该种记录里的 /type_material 声明值（形如 culture from holotype of X）。"""
        return self._qual_values(species, "type_material", retmax)

    def retrievability(self, code, number):
        """逐字段探一探：这个号在不同检索字段里能不能查到。

        用来回答「是不是我检索字段选窄了」——例如号只出现在 /culture_collection 里、
        而我只查了 [Strain]。三个字段各自计数，互不替代。
        """
        out = {}
        base = self.num_term(code, number)
        for tag, term in (("all", base),
                          ("strain", "{} AND {}[Strain]".format(base, base)),
                          ("culture_collection",
                           "{} AND {}[Culture Collection]".format(base, base))):
            try:
                out[tag] = self._count(term)
            except Exception:
                out[tag] = None
        return out

    def note_type_status(self, accession):
        """某条记录的 /note 里有没有 "Type status"（J-20 用）。"""
        body, _ = self._records_for('{}[Accession]'.format(accession), 3)
        out = []
        for rec in self._parse_records(body):
            for v in re.findall(r'/note="([^"]*)"', rec):
                if "type status" in v.lower():
                    out.append(v.strip())
        return out


    def type_material_text(self, code, number, retmax=30):
        """该号在 GenBank 记录里的 /type_material 原文（带缓存）。

        ⚠️ 检索式必须用【连写形】（ICMP1783），不能用带空格的引号形
           （"ICMP 1783"）。实测后者被 Entrez 拆成两个词，命中 20,561 条无关记录，
           取前几十条自然找不到真记录 → 误判成"记录没有该字段"。
           用连写形：ICMP1783 只命中 6 条，全部带 culture from epitype of ...
        另外再按 /culture_collection 里是否含该号做一次核对，防串号。
        """
        import time

        def _norm(x):
            return re.sub(r"[^A-Z0-9]", "", str(x).upper())

        want = _norm("{}{}".format(code, number))
        tag = re.sub(r"[^A-Za-z0-9]", "_", "{}_{}".format(code, number))
        key = "tmtext2_{}.txt".format(tag)

        def mk():
            # ⚠️ 不要"取命中最多的一批"。实测同一批检索式：
            #     "ICMP 1783"（带引号空格）→ 20,561 条（拆词，全是无关记录）
            #     ICMP:1783（不带引号）    → 53 条（含真记录）
            #     ICMP1783                 → 6 条（全是真记录）
            #   取"最多"就会选中那个两万条的假集合，真记录一条都进不来。
            #   正确做法：用不带引号的冒号式（num_term）取，再逐条核对正文里有没有这个号。
            best_ids = []
            try:
                xml = self._esearch("nucleotide", self.num_term(code, number), retmax=retmax)
                best_ids = re.findall(r"<Id>(\d+)</Id>", xml)
            except Exception:
                best_ids = []
            if not best_ids:
                try:
                    time.sleep(self.pause)
                    xml = self._esearch("nucleotide",
                                        "{}{}".format(code, re.sub(r"[^A-Za-z0-9]", "", str(number))),
                                        retmax=retmax)
                    best_ids = re.findall(r"<Id>(\d+)</Id>", xml)
                except Exception:
                    best_ids = []
            if not best_ids:
                return ""
            time.sleep(self.pause)
            return _rd(_entrez(self.email).efetch(
                db="nucleotide", id=",".join(best_ids), rettype="gb", retmode="text"))

        body = self._cached(key, mk)
        if not body:
            return []
        out = []
        for rec in body.split("\n//\n"):
            org = re.search(r'/organism="([^"]+)"', rec)
            # 核对：这条记录真的提到这个号吗（culture_collection / strain / 正文）
            blob = _norm(rec)
            if want and want not in blob:
                continue
            for s in re.findall(r'/type_material="([^"]+)"', rec):
                out.append({"organism": org.group(1) if org else "", "text": s.strip()})
        return out


class CGMCCAPI(Source):
    """CGMCC 官网菌株目录。

    站点性质（摘自 probe_cgmcc.py 的实测记录）：
      www.cgmcc.net 是 Vue SPA，页面 HTML 对所有路由返回同一份空壳，数据全走 XHR。
      API 藏在 /static/js/app.*.js 里：
        base = https://www.cgmcc.net/cgmccapi
        端点 = /commom/resource/cataloguelist    ← 注意站点把 common 拼成了 commom
      参数名：cgmcc_strain_number / latinname / species_name / strain_type / paginate / pageSize
      返回关键字段：is_typestrain = "模式菌株" / "非模式菌株" / …

    用途：J-16 的目录一侧（记录声明 vs 官网标注）。
    实测例：`CGMCC 3.17361` 记录里当模式用，而官网记「非模式菌株」（本项目 E2）。

    ⚠️ 限流纪律：第三方站点 ≥1.2 秒/请求；结果落盘缓存，重跑零请求；诚实 User-Agent。
    """
    id = "cgmcc"
    kind = "api"
    BASE = "https://www.cgmcc.net/cgmccapi"
    ENDPOINT = "/commom/resource/cataloguelist"
    UA = {"User-Agent": "typeaudit/0.1 (academic metadata verification)",
          "Referer": "https://www.cgmcc.net/resources",
          "Accept": "application/json, text/plain, */*"}
    license_note = "CGMCC 官网公开查询；使用前请确认其条款"

    def __init__(self, cache_dir=".cache", pause=1.3):
        # ⚠️ 配置里这一项没有 path（只有 _端点），于是 load() 传进来的是空串或 URL。
        #    直接当目录用的话：空串 ⇒ 缓存静默失效（每次查询都打官网）；
        #    URL ⇒ makedirs 失败、整个源被 try/except 悄悄丢掉。
        #    所以这里把「不像路径」的值统一落到本地缓存目录。
        if (not cache_dir) or ("://" in str(cache_dir)):
            cache_dir = os.path.join(".cache", "cgmcc")
        self.cache_dir, self.pause = cache_dir, pause
        if cache_dir:
            os.makedirs(cache_dir, exist_ok=True)

    def _query(self, number, page_size=10):
        import json as _json, time, urllib.parse, urllib.request, urllib.error
        key = "cgmcc_{}.json".format(re.sub(r"\W+", "_", str(number)))
        if self.cache_dir:
            p = os.path.join(self.cache_dir, key)
            if os.path.exists(p):
                return _json.load(io.open(p, encoding="utf-8"))
        params = {"cgmcc_strain_number": number, "paginate": 1, "pageSize": page_size}
        url = self.BASE + self.ENDPOINT + "?" + urllib.parse.urlencode(params)
        try:
            time.sleep(self.pause)
            r = urllib.request.urlopen(urllib.request.Request(url, headers=self.UA), timeout=25)
            j = _json.loads(r.read().decode("utf-8", "replace"))
        except urllib.error.HTTPError as e:
            j = {"error": "HTTP{}".format(e.code)}
        except Exception as e:
            j = {"error": type(e).__name__}
        if self.cache_dir:
            io.open(os.path.join(self.cache_dir, key), "w", encoding="utf-8").write(
                _json.dumps(j, ensure_ascii=False))
        return j

    def strain_of(self, number):
        """按 CGMCC 编号查目录。number 形如 '3.17361' 或 'CGMCC 3.17361'。"""
        m = re.search(r"(\d+\.\d+)", str(number))
        num = m.group(1) if m else str(number)
        j = self._query(num)
        if not isinstance(j, dict) or j.get("error"):
            return {"number": num, "found": False, "error": (j or {}).get("error")}
        if j.get("status") != 0:
            return {"number": num, "found": False, "error": j.get("msg", "status!=0")}
        d = j.get("data") or {}
        lst = d.get("list") or []
        if not lst:
            return {"number": num, "found": False, "total": d.get("totalRow", 0)}
        it = lst[0]
        return {"number": num, "found": True, "total": d.get("totalRow", 0),
                "is_typestrain": (it.get("is_typestrain") or "").strip(),
                "latinname": (it.get("latinname") or "").strip(), "raw": it}

    def is_type_strain(self, number):
        """True / False / None（查不到）。"""
        r = self.strain_of(number)
        if not r.get("found"):
            return None
        t = r.get("is_typestrain") or ""
        if "非模式" in t:
            return False
        if "模式" in t:
            return True
        return None


# ==================== 注册表 ====================
# ⚠️ 必须放在**所有类定义之后**——本字典在导入时求值，若某个 adapter 类尚未定义会 NameError。
#    （本项目在 rules.py 的 RULES 与这里的 KINDS 各踩过一次，同一类坑。）
# ---------------- 1c. ICMP（新西兰 Landcare Research 的 SCD 公开库） ----------------
    def strains_of_name(self, latinname, page_size=20, only_type=False):
        """按【学名】查 CGMCC 官网目录，返回该名下的菌株列表。

        ⚠️ 必须用 latinname 参数（= 完整学名），不要用 species_name：
           实测 species_name=camelliae 会跨属命中 8 条（Aeromicrobium camelliae 等），
           因为那只是种加词；latinname 才是完整学名。

        only_type=True 时加 strain_type=模式菌株 —— 用来回答
        「这家保藏机构说这个种有模式菌株吗」。

        返回 [{"number","is_typestrain","accession","latinname"}]；查不到返回 []。
        """
        import json as _json, time, urllib.parse, urllib.request, urllib.error
        key = "cgmccnm_{}{}.json".format(re.sub(r"\W+", "_", str(latinname))[:60],
                                         "_type" if only_type else "")
        def _keep_type(rows):
            """只保留 is_typestrain 为「模式菌株」的行（注意排除「非模式菌株」）。"""
            return [x for x in rows
                    if "模式菌株" in (x.get("is_typestrain") or "")
                    and "非" not in (x.get("is_typestrain") or "")]

        if self.cache_dir:
            p = os.path.join(self.cache_dir, key)
            if os.path.exists(p):
                try:
                    _cached_rows = _json.load(io.open(p, encoding="utf-8"))
                    # ⚠️ 缓存里存的是【未筛】全量（同一个 key 供两种调用），
                    #    所以命中缓存后必须同样过筛 —— 否则第一次查有筛、第二次没筛。
                    return _keep_type(_cached_rows) if only_type else _cached_rows
                except Exception:
                    pass
        # ⚠️ 不要用 strain_type=模式菌株 这个参数：实测服务端不认，静默返回 0 条
        #    （siamense 明明有 CGMCC 3.14174 标着「模式菌株」，带该参数却 0 条）。
        #    改为取全量后在客户端按 is_typestrain 字段筛。
        params = {"latinname": latinname, "paginate": 1, "pageSize": page_size}
        url = self.BASE + self.ENDPOINT + "?" + urllib.parse.urlencode(params)
        out = []
        try:
            time.sleep(self.pause)
            r = urllib.request.urlopen(urllib.request.Request(url, headers=self.UA), timeout=30)
            d = (_json.loads(r.read().decode("utf-8", "replace")).get("data") or {})
            for x in (d.get("list") or []):
                out.append({"number": x.get("cgmcc_strain_number") or "",
                            "is_typestrain": x.get("is_typestrain") or "",
                            "accession": x.get("genbank_accession_no") or "",
                            "latinname": x.get("latinname") or ""})
        except Exception:
            out = []
        if self.cache_dir:
            try:
                os.makedirs(self.cache_dir, exist_ok=True)
                _json.dump(out, io.open(os.path.join(self.cache_dir, key), "w",
                                        encoding="utf-8", newline="\n"), ensure_ascii=False)
            except Exception:
                pass
        return _keep_type(out) if only_type else out


class ICMPAPI(Source):
    """
    ICMP（International Collection of Microorganisms from Plants）的公开库 SCD。

    【为什么值得单独接 —— 它给的是 CBS 目录给不了的东西】
        · **独立的 Type status 字段**（Holotype / …）。实测 CBS 目录 62,740 行里
          只有 15,098 行（24.1%）有模式地位标注，ICMP 这条轴能补上那个空洞。
        · Host / Collector / Collection date / 标准产地 —— 可查「宿主与种是否相符」这类错。
        · 它是**第三条独立权威**：同一株菌可以三方比名
          （GenBank 声明 ｜ CBS 目录 ｜ ICMP 记录）。实测 ICMP 18578 三方一致。
        许可：NZGOAL（新西兰政府开放获取与许可框架），七个国家级收藏统一提供。
        ⚠️ 该站条款自己声明「数据可能不准确、不完整、不精确」—— 那正是本工具的题材。

    【接口（本机实测，2026-10）】
        GET /Search/ResultsPartial?collectionId=ICMP&query=<号>
            &pageNumber=0&sortField=relevance&currentDisplayTab=list
        返回的是 **HTML 而不是 JSON**，列表表格的 <th> 表头稳定：
            Accession number / Taxonomic name / Type status / Country / … / Host
        详情页：/Specimen/ICMP%20<号>
        检索语法：空格=AND、双引号=短语、`18578*`=通配、字段检索 `taxonName:Colletotrichum`

    【限流】≥1.5 s/请求；结果落盘缓存，重跑零请求。

    ⚠️ 解析失败时**明确报错**（{"error": ...}），不静默返回空 —— 页面改版会断，
       而「静默返回空」会伪装成「号不存在」，那比报错危险得多。
    """
    id = "icmp"
    kind = "api"
    BASE = "https://scd.landcareresearch.co.nz"
    ENDPOINT = "/Search/ResultsPartial"
    UA = {"User-Agent": "colleto-type-audit/1.0 (academic metadata review; contact: audit@example.org)",
          "Accept": "*/*", "X-Requested-With": "XMLHttpRequest",
          "Referer": "https://scd.landcareresearch.co.nz/Search"}
    license_note = "Landcare Research SCD；按 NZGOAL 开放获取与许可框架提供"

    _RE_TABLE = re.compile(r"<table[^>]*>(.*?)</table>", re.S | re.I)
    _RE_TR = re.compile(r"<tr[^>]*>(.*?)</tr>", re.S | re.I)
    _RE_TH = re.compile(r"<th[^>]*>(.*?)</th>", re.S | re.I)
    _RE_TD = re.compile(r"<td[^>]*>(.*?)</td>", re.S | re.I)

    # ---- 跨机构号：Duplicates 字段（ICMP 这一侧的桥）----
    RE_DUP = re.compile(r"^(.*?):\s*(.+)$")

    def duplicates_of(self, code, number):
        """从 ICMP 详情页解析 Duplicates（重复材料）—— 同一份材料在其他机构的号。

        为什么必须加它：CBS 的 other_coll 并不是全的。实测 `CBS 125378` 在 CBS 一侧
        别名数为 **0**，但 ICMP 18642 的详情页把它和 `GZAAS 080001` 一起写在 Duplicates 里
        ⇒ **ICMP 能补 CBS 的缺口**，两座桥互补。

        页面结构（实测，别改成别的选择器）：
            <div class="span12 major-heading">Duplicates</div>
            之后每条形如 <div class="span12 field-heading">机构全名: 号码</div>
            （号码串自带机构代码，如 "CBS 130417"；机构名只作参考）

        返回 [{"institution","number","code"}]；没有该段返回 []。
        """
        import json as _json, time, urllib.request
        key = "icmpdup_{}_{}.json".format(re.sub(r"\W+", "_", str(code))[:20],
                                          re.sub(r"\W+", "_", str(number))[:30])
        cd = getattr(self, "cache_dir", None)
        if cd:
            p = os.path.join(cd, key)
            if os.path.exists(p):
                try:
                    return _json.load(io.open(p, encoding="utf-8"))
                except Exception:
                    pass
        url = self.BASE + "/Specimen/ICMP%20{}".format(str(number).strip())
        out = []
        try:
            time.sleep(1.5)                     # 限流：这是政府站点，别打疼人家
            with urllib.request.urlopen(urllib.request.Request(url, headers=self.UA),
                                        timeout=40) as r:
                html = r.read().decode("utf-8", "replace")
            parts = re.split(r'<div class="span12 major-heading[^"]*"\s*>', html)
            for seg in parts:
                head = re.sub(r"(?s)<[^>]+>", " ", seg.split("</div>")[0]).strip().lower()
                if head != "duplicates":
                    continue
                for m in re.finditer(r'<div class="span12 field-heading"\s*>(.*?)</div>',
                                     seg, re.S):
                    t = re.sub(r"(?s)<[^>]+>", " ", m.group(1))
                    t = re.sub(r"&nbsp;?", " ", t)
                    t = re.sub(r"\s+", " ", t).strip()
                    mm = self.RE_DUP.match(t)
                    if not mm:
                        continue
                    num = mm.group(2).strip()
                    cm = re.match(r"^([A-Za-z][A-Za-z\.\-]{0,12})\s*[:#]?\s*([\w\.\-]+)$", num)
                    # ⚠️ 必须校验：实测有条目只写了机构名没有号（如 "…: ICMP"），
                    #    不过滤就会产出裸 "ICMP" 这种脏号，进而污染连通分组。
                    if not cm or not re.search(r"\d", cm.group(2)) or len(cm.group(2)) < 1:
                        continue
                    # 号码部分不能只是机构代码本身
                    if cm.group(2).upper() == cm.group(1).upper():
                        continue
                    out.append({"institution": mm.group(1).strip(),
                                "number": num,
                                "code": cm.group(1).upper()})
                break
        except Exception:
            out = []
        if cd:
            try:
                os.makedirs(cd, exist_ok=True)
                _json.dump(out, io.open(os.path.join(cd, key), "w", encoding="utf-8",
                                        newline="\n"), ensure_ascii=False)
            except Exception:
                pass
        return out



    def __init__(self, cache_dir=".cache", pause=1.5):
        if (not cache_dir) or ("://" in str(cache_dir)):
            cache_dir = os.path.join(".cache", "icmp")
        self.cache_dir, self.pause = cache_dir, pause
        if cache_dir:
            os.makedirs(cache_dir, exist_ok=True)

    @staticmethod
    def _txt(frag):
        import html as _html
        return re.sub(r"\s+", " ", _html.unescape(re.sub(r"<[^>]+>", " ", frag))).strip()

    @staticmethod
    def _bare(x):
        """归属一号：去掉集合前缀与所有非字母数字，便于比对。"""
        return re.sub(r"[^0-9A-Za-z]", "", re.sub(r"^ICMP", "", str(x or ""), flags=re.I)).upper()

    def _parse_rows(self, doc):
        """从结果页里取记录：只认表头含 'Accession number' 的那张表。"""
        out = []
        for tbl in self._RE_TABLE.findall(doc):
            rows = self._RE_TR.findall(tbl)
            if not rows:
                continue
            heads = []
            for r in rows:
                hs = [self._txt(x) for x in self._RE_TH.findall(r)]
                if hs and "Accession number" in hs:
                    heads = hs
                    break
            if not heads:
                continue
            for r in rows:
                tds = [self._txt(x) for x in self._RE_TD.findall(r)]
                if not tds:
                    continue
                if len(tds) == len(heads):
                    out.append(dict(zip(heads, tds)))
                elif len(tds) == len(heads) - 1:
                    # 首列常是勾选框（表头有 th、数据行没有 td）
                    out.append(dict(zip(heads[1:], tds)))
        return out

    def _query(self, number, timeout=45):
        import gzip, json as _json, time, urllib.parse, urllib.request, urllib.error
        key = "icmp_{}.json".format(re.sub(r"\W+", "_", str(number)))
        cp = os.path.join(self.cache_dir, key) if self.cache_dir else None
        if cp and os.path.exists(cp):
            try:
                return _json.loads(io.open(cp, encoding="utf-8").read())
            except Exception:
                pass
        url = self.BASE + self.ENDPOINT + "?" + urllib.parse.urlencode({
            "collectionId": "ICMP", "query": str(number),
            "pageNumber": 0, "sortField": "relevance", "currentDisplayTab": "list"})
        time.sleep(self.pause)
        try:
            r = urllib.request.urlopen(urllib.request.Request(url, headers=self.UA), timeout=timeout)
            raw = r.read()
            if r.headers.get("Content-Encoding") == "gzip":
                raw = gzip.decompress(raw)
            doc = raw.decode("utf-8", "replace")
            rows = self._parse_rows(doc)
            j = {"rows": rows, "n": len(rows), "html_len": len(doc), "url": url}
            if not rows:
                # 明确记下「页面在但解析不出记录」—— 这可能是真没记录，也可能是页面改版
                j["warn"] = "页面已取到但解析不出记录行（可能真的查无，也可能页面结构变了）"
        except urllib.error.HTTPError as e:
            j = {"error": "HTTP{}".format(e.code), "url": url}
        except Exception as e:
            j = {"error": type(e).__name__, "url": url}
        if cp:
            try:
                io.open(cp, "w", encoding="utf-8").write(_json.dumps(j, ensure_ascii=False))
            except Exception:
                pass
        return j

    def claims(self, code):
        return (code or "").strip().upper().startswith("ICMP")

    def status_of(self, code, number):
        j = self._query(number)
        if j.get("error"):
            return {"rank": None, "name": None, "raw": "", "unannotated": True,
                    "api_error": j["error"]}
        rows = j.get("rows") or []
        if not rows:
            return None
        want = self._bare(number)
        hit = None
        for r in rows:
            if self._bare(r.get("Accession number")) == want:
                hit = r
                break
        if hit is None and len(rows) == 1:
            hit = rows[0]
        if hit is None:
            return None
        name = (hit.get("Taxonomic name") or "").strip()
        ts = (hit.get("Type status") or "").strip()
        return {"name": name, "type_status": ts, "raw": ts or name,
                "unannotated": (not name),
                "country": hit.get("Country") or "", "host": hit.get("Host") or "",
                "collector": hit.get("Collector") or "",
                "record": hit}


# ---------------- 8. 跨库号映射（CBSMap） ----------------
class CBSMap(Source):
    """跨库号映射：别家机构的号 → CBS 号。

    kind = offline：读一份 CBS 导出 TSV（含 other_coll 列，如 "BRIP 15842a"）。

    【为什么需要】一管菌常同时存放在几家机构（CBS 151555 = BRIP 15842a），
    GenBank 提交者往往只填手边那个号。若只查 CBS 目录，「GenBank 无声明 CBS 号」
    会产生大量假阳性（本机实测 6 条里 3 条）。

    【三态纪律】沿用基类 code_map：
        None          = 本源没有跨库映射能力（本类永远返回 dict，故不出现）
        {"matches": []} = 查了，没有这条映射
        {"matches":[...]} = 查到了

    【闸门】号必须【整体相等】：BRIP 63641 与 BRIP 63641a 只差一个尾缀，
    用"包含"判断会把两管菌当成一管（实测踩过）。
    """
    id = "cbsmap"
    kind = "offline"
    license_note = "Westerdijk Institute 公开查询页导出；使用前请确认其条款"

    def __init__(self, path):
        self.rows = []
        self.by_other = {}
        if path and os.path.exists(path):
            for ln in io.open(path, encoding="utf-8-sig", errors="replace").read().splitlines()[1:]:
                if not ln.strip():
                    continue
                c = ln.split("\t")
                if len(c) < 8:
                    continue
                row = {"cbs": c[0].strip(), "taxon": (c[4] or "").strip(),
                       "organism": (c[6] or "").strip(), "other_coll": (c[7] or "").strip()}
                self.rows.append(row)
                for part in re.split(r"[;|,]", row["other_coll"]):
                    part = part.strip()
                    if not part:
                        continue
                    m = re.match(r"^\s*([A-Za-z][A-Za-z\.\-]{0,12})\s*[:#]?\s*([\w\.\-]+)\s*$", part)
                    if not m:
                        continue
                    k = (m.group(1).upper(), m.group(2).upper().replace(" ", ""))
                    _num = row["cbs"].replace("CBS", "").strip()
                    self.by_other.setdefault(k, []).append(
                        {"code": "CBS", "number": _num,
                         "display": "CBS {}".format(_num),
                         "taxon": row["taxon"], "via": "other_coll"})

    def code_map(self, code, number):
        """别家号 → CBS 号。两条路：本地导出表（other_coll）+ 联网查询缓存。

        ⚠️ 只读本地表会漏。实测：本地 5,001 行的导出表里 **没有任何 MFLU 记录**，
        但 `.cache/cbsmap/q_MFLUCC_10_0277.json` 里明确存着
        "MFLUCC 10-0277 → CBS 130420"（那是原实现联网查 CBS 后缓存下来的）。
        两条路都要走，否则 F31/F32 这两条已知错误查不出来。
        """
        c = (code or "").strip().upper()
        n = (number or "").strip()
        out = []

        # ① 本地导出表（other_coll 列）
        k = (c, n.upper().replace(" ", ""))
        out += list(self.by_other.get(k, []))

        # ② 联网查询缓存 .cache/cbsmap/q_<CODE>_<NUM>.json
        try:
            tag = "q_{}_{}.json".format(re.sub(r"[^A-Za-z0-9]", "_", c),
                                        re.sub(r"[^A-Za-z0-9]", "_", n))
            cp = os.path.join(_CBSMAP_CACHE, tag)
            if os.path.exists(cp):
                for it in (json.loads(io.open(cp, encoding="utf-8", errors="replace").read()) or []):
                    cbs = str(it.get("cbs") or "").strip()
                    taxon = str(it.get("taxon") or "")
                    # taxon 有时是对象字符串（里面嵌着 'Value': 'Colletotrichum siamense'）
                    m = re.search(r"'Value':\s*'([^']+)'", taxon)
                    if m:
                        taxon = m.group(1)
                    else:
                        taxon = taxon.split("|")[0].strip()
                    num = re.sub(r"^CBS\s*", "", cbs).strip()
                    if num:
                        out.append({"code": "CBS", "number": num,
                                    "display": "CBS {}".format(num),
                                    "taxon": taxon, "via": "cbs_web_cache"})
        except Exception:
            pass

        # 去重（同一 CBS 号可能两条路都给）
        seen, uniq = set(), []
        for m in out:
            key = (m.get("code"), str(m.get("number")))
            if key in seen:
                continue
            seen.add(key)
            uniq.append(m)
        return {"matches": uniq}

    @property
    def n_rows(self):
        return len(self.rows)


# ---------------- 9. NCBI Taxonomy ----------------
class TaxonomyAPI(_EntrezBase):
    """NCBI Taxonomy（db=taxonomy）：查 taxid、种清单、类型材料声明。

    提供：
      taxid_of(name)                  学名 → taxid（精确匹配才返回，不做模糊）
      species_in_genus(genus)         该属下的种清单（只取双名，剔除 sp./cf. 之类噪声）
      type_material_of(species)       → (taxid, ["CBS 112996", ...]) 类型材料声明里的号
      taxon_info(name)                → {taxid, scientific_name, rank, authority}
    """
    id = "taxonomy"

    def __init__(self, path=None, email=None, pause=None):
        super(TaxonomyAPI, self).__init__(
            cache_dir=(path or ".cache"),
            email=(email or os.environ.get("NCBI_EMAIL") or None),
            pause=(float(pause) if pause else 1.3))

    def _fetch_taxon(self, taxid):
        """取一条分类学记录的 XML 原文（带缓存）。"""
        def mk():
            import time as _t
            _t.sleep(self.pause)
            return _rd(_entrez(self.email).efetch(
                db="taxonomy", id=str(taxid), retmode="xml"))
        return self._cached("tax_{}.xml".format(taxid), mk)

    def taxid_of(self, name):
        """学名 → taxid。只有【学名精确相同】才认，避免把近名种认错。"""
        name = (name or "").strip()
        if not name:
            return None
        try:
            xml = self._esearch("taxonomy", '"{}"'.format(name), retmax=5)
            ids = re.findall(r"<Id>(\d+)</Id>", xml)
        except Exception:
            return None
        want = re.sub(r"\s+", " ", name).strip().lower()
        for tid in ids:
            try:
                txt = self._fetch_taxon(tid)
            except Exception:
                continue
            sci = re.search(r"<ScientificName>([^<]+)</ScientificName>", txt)
            if sci and re.sub(r"\s+", " ", sci.group(1)).strip().lower() == want:
                return tid
        return ids[0] if len(ids) == 1 else None

    def taxon_info(self, name):
        """学名 → {taxid, scientific_name, rank, authority}；查不到返回 None。"""
        tid = self.taxid_of(name)
        if not tid:
            return None
        try:
            txt = self._fetch_taxon(tid)
        except Exception:
            return None
        sci = re.search(r"<ScientificName>([^<]+)</ScientificName>", txt)
        rank = re.search(r"<Rank>([^<]+)</Rank>", txt)
        auth = re.search(r"<Authority>([^<]+)</Authority>", txt)
        return {"taxid": tid,
                "scientific_name": (sci.group(1).strip() if sci else ""),
                "rank": (rank.group(1).strip() if rank else ""),
                "authority": (auth.group(1).strip() if auth else "")}

    def species_in_genus(self, genus, retmax=3000, limit=None):
        """该属下的种清单（只保留规范双名）。

        为什么要过滤：NCBI 里同一个属会有大量 `Colletotrichum sp.` / `cf.` / 菌株代号，
        直接拿来逐种跑会白跑上千个不存在的"种"（本机实测 2,277 个 taxon 里只有 408 个是双名）。
        """
        term = '{}[Organism] AND species[Rank]'.format(genus)
        try:
            xml = self._esearch("taxonomy", term, retmax=retmax)
            ids = re.findall(r"<Id>(\d+)</Id>", xml)
        except Exception:
            return []
        out, seen = [], set()
        for tid in ids:
            try:
                txt = self._fetch_taxon(tid)
            except Exception:
                continue
            sci = re.search(r"<ScientificName>([^<]+)</ScientificName>", txt)
            if not sci:
                continue
            nm = re.sub(r"\s+", " ", sci.group(1)).strip()
            parts = nm.split()
            if len(parts) != 2:
                continue
            if parts[1].lower() in ("sp", "sp.", "spp", "spp.", "cf", "cf.", "aff", "aff."):
                continue
            if parts[1][0].islower() is False:      # 种加词必须小写
                continue
            if nm.lower() in seen:
                continue
            seen.add(nm.lower())
            out.append(nm)
            if limit and len(out) >= int(limit):
                break
        return out

    def type_material_of(self, species, retmax=3):
        """该种在 NCBI Taxonomy 里登记的类型材料名（如 'CBS 127945'）。

        ⚠️ 真实 XML 结构是（**不是** <TypeMaterial> 标签）：
            <Name><ClassCDE>type material</ClassCDE><DispName>CBS 127945</DispName></Name>
          同一号常以空格形与冒号形各出现一次（'CBS 127945' / 'CBS:127945'），
          故返回值需去重。按 <TypeMaterial> 解析会**静默得到空列表**（本项目踩过）。
        """
        xml = self._esearch("taxonomy", '"{}"'.format(species), retmax=retmax)
        ids = re.findall(r"<Id>(\d+)</Id>", xml)
        if not ids:
            return None, []
        def mk():
            import time
            time.sleep(self.pause)
            return _rd(_entrez(self.email).efetch(db="taxonomy", id=",".join(ids), rettype="xml"))
        x = self._cached("tmt_{}.xml".format(re.sub(r"\W+", "_", species)), mk)
        names = list(dict.fromkeys(re.findall(
            r"<ClassCDE>type material</ClassCDE>\s*<DispName>([^<]+)</DispName>", x, re.I)))
        # ⚠️ XML 实体必须解码：实测有 `BCC&lt;TH…` 这种写法，
        #    不解码会与对账表里的 `BCC<TH…` 对不上（6 行"差异"全由此产生）。
        def _un(s):
            for a, b in (("&lt;", "<"), ("&gt;", ">"), ("&amp;", "&"), ("&quot;", '"'), ("&apos;", "'")):
                s = s.replace(a, b)
            return re.sub(r"&#x([0-9A-Fa-f]+);", lambda m: chr(int(m.group(1), 16)), s)
        names = [_un(n) for n in names]
        return (ids[0] if ids else None), names



# ---------------- 5. CGMCC 官网（保藏中心官网，非 NCBI） ----------------

    def organism_of_number(self, code, number, retmax=5):
        """该号在 NCBI 里被登记在哪些种名下（旁证，不参与判档）。"""
        try:
            xml = self._esearch("taxonomy", '"{} {}"'.format(code, number), retmax=retmax)
            ids = re.findall(r"<Id>(\d+)</Id>", xml)
        except Exception:
            return []
        names = []
        for tid in ids:
            try:
                d = self.taxon_info_by_id(tid)
            except Exception:
                continue
            if d:
                names.append(d)
        return names

    def taxon_info_by_id(self, taxid):
        try:
            txt = self._fetch_taxon(taxid)
        except Exception:
            return None
        sci = re.search(r"<ScientificName>([^<]+)</ScientificName>", txt)
        rank = re.search(r"<Rank>([^<]+)</Rank>", txt)
        return {"taxid": str(taxid),
                "scientific_name": (sci.group(1).strip() if sci else ""),
                "rank": (rank.group(1).strip() if rank else "")}


KINDS = {
    "ICMPAPI": ICMPAPI,
    "CBSCatalogue": CBSCatalogue,
    "CBSMap": CBSMap,            # 跨库号映射（需 out/cbs_live_*.tsv 这种带 other_coll 的导出）
    "Biocollections": Biocollections,
    "MycoBank": MycoBank,
    "GenBankAPI": GenBankAPI,
    "TaxonomyAPI": TaxonomyAPI,
    "CGMCCAPI": CGMCCAPI,
}


def authority_verdict(name_a, name_b, tx, mb):
    """两家权威对账：NCBI Taxonomy（taxid 是否同一个）vs MycoBank（接受名是否同一个）。

    为什么要这一步：异名判定过去只问 MycoBank 一家。当两家说法相反时
    （实测：NCBI 保留 Colletotrichum analogum / C. celtidis 为独立种，各自有 taxid；
      MycoBank 把它们分别并入 camelliae / karsti），单问一家会把分歧误判成「错值」。

    返回：
        {"ncbi": "same|diff|unknown", "mycobank": "same|diff|unknown",
         "verdict": "一致|分歧|无法对账",
         "ncbi_label": "...", "mycobank_label": "..."}

    ⚠️ 「unknown」（查无）必须与「diff」（不同）分开 —— 这与本模块 code_map 的三态是同一条纪律。
    """
    out = {"ncbi": "unknown", "mycobank": "unknown", "verdict": "无法对账",
           "ncbi_label": "查无", "mycobank_label": "查无"}
    try:
        ta = tx.taxid_of(name_a) if tx is not None else None
        tb = tx.taxid_of(name_b) if tx is not None else None
    except Exception:
        ta = tb = None
    if ta and tb:
        out["ncbi"] = "same" if ta == tb else "diff"
        out["ncbi_label"] = "同 taxid({})".format(ta) if ta == tb else "独立({}≠{})".format(ta, tb)
    try:
        cur = mb.current_name(name_a) if mb is not None else None
    except Exception:
        cur = None
    if cur:
        w = str(cur).split()
        acc = "{} {}".format(w[0], w[1]) if len(w) >= 2 else str(cur)
        same = acc.split()[-1].lower() == str(name_b).split()[-1].lower()
        out["mycobank"] = "same" if same else "diff"
        out["mycobank_label"] = "异名→{}".format(acc.split()[-1]) if same else "独立({})".format(acc.split()[-1])
    if out["ncbi"] != "unknown" and out["mycobank"] != "unknown":
        out["verdict"] = "一致" if out["ncbi"] == out["mycobank"] else "分歧"
    return out
# ---------------------------------------------------------------- 种级四方核验
def species_status(name, sources, use_cgmcc=True, weights=None):
    """【一个种，四处问状态，再对账】MycoBank / NCBI / CBS 目录 / CGMCC 官网。

    ⚠️ 权重与「空白不计错」是两条硬规则（用户定）：
        · CBS / NCBI / MycoBank 权重 3；**CGMCC 权重 1**（收录不全，少一条不能当证据）
        · **某一方"没有"永远是空白，不是反证** —— 只有「双方都正面声明且互斥」才算矛盾。
          因此不再因为「CGMCC 查不到」就报冲突；那只会进 notes。

    返回四类结论：
        number_shared   同一个号在几方出现（回答「号是否相同」的一半）
        type_conflict   模式菌株 vs 普通菌株 的正面冲突
        name_conflict   异名 vs 独立种 的正面冲突（NCBI taxid ↔ MycoBank 接受名）
        conflicts       上面两类里【权重≥2 的两方】对撞出来的，才算矛盾
        notes           空白/未收录/低权重差异 —— 一律不计错
    """
    W = {"myco": 3, "ncbi": 3, "cbs": 3, "cgmcc": 1}
    if weights:
        W.update(weights)
    def _orig_num(_out, _k):
        """把归一化号还原成「代码 号码」形式，供 code_map 使用。"""
        for x in (_out["cgmcc"].get("rows") or []):
            if _nk(x["number"]) == _k:
                return str(x["number"]).replace("CGMCC", "").strip()
        return _k

    out = {"name": name, "weights": W, "pending": [],
           "myco": {"ok": False}, "ncbi": {"ok": False},
           "cbs": {"ok": False}, "cgmcc": {"ok": False},
           "number_shared": {}, "type_conflict": [], "name_conflict": None,
           "conflicts": [], "notes": []}

    def _nk(x):
        return re.sub(r"[\s:<>\.]", "", str(x or "")).upper()

    # ---- MycoBank ----
    mb = sources.get("mycobank") if sources else None
    try:
        if mb is not None and getattr(mb, "n_rows", lambda: 0)() > 0:
            cur = mb.current_name(name) or ""
            # Name status 必须一起取：Orthographic variant（拼写变体）与真异名
            # 在「接受名不同」这件事上长得一样，只有 status 能把两者分开。
            _st = ""
            try:
                _st = (mb.status_of(name) or "").strip()
            except Exception:
                _st = ""
            out["myco"] = {"ok": bool(cur), "current": cur, "status": _st}
        else:
            out["notes"].append("MycoBank 未加载（空表）—— 不参与判定")
    except Exception:
        out["notes"].append("MycoBank 查询失败 —— 不参与判定")

    # ---- NCBI Taxonomy ----
    tx = sources.get("taxonomy") if sources else None
    try:
        if tx is not None:
            tid = tx.taxid_of(name)
            if tid:
                tm = []
                try:
                    r = tx.type_material_of(name, retmax=3)
                    tm = list(r[1]) if isinstance(r, tuple) else list(r or [])
                except Exception:
                    tm = []
                out["ncbi"] = {"ok": True, "taxid": tid, "taxid_label": "taxid({})".format(tid),
                               "type_material": tm}
    except Exception:
        pass

    # ---- CBS 目录（严格：名字 ∧ 类型词）----
    cbs = sources.get("cbs") if sources else None
    try:
        if cbs is not None and getattr(cbs, "n_rows", lambda: 0)() > 0:
            TW = (r"(?i)\b(type|holotype|isotype|neotype|epitype|paratype|syntype|"
                  r"lectotype|ex-?type|ex-?holotype|authentic)\b")
            typed = []
            for r in getattr(cbs, "rows", []):
                st = (r.get("status") or "").strip()
                if not st:
                    continue
                if not re.search(r"\b" + re.escape(name) + r"\b", st, re.I):
                    continue
                if not re.search(TW, st):
                    continue
                typed.append({"number": r.get("cbs_number") or "", "status": st})
            out["cbs"] = {"ok": True, "type_numbers": typed}
    except Exception:
        pass

    # ---- CGMCC 官网 ----
    if use_cgmcc:
        cg = sources.get("cgmcc") if sources else None
        try:
            if cg is not None and hasattr(cg, "strains_of_name"):
                rows = cg.strains_of_name(name)
                out["cgmcc"] = {"ok": True, "rows": rows,
                                "type_rows": [x for x in rows
                                              if "模式菌株" in (x.get("is_typestrain") or "")
                                              and "非" not in (x.get("is_typestrain") or "")]}
        except Exception:
            out["notes"].append("CGMCC 查询失败 —— 不参与判定")

    # ---- ① 号是否在多方出现 ----
    where = {}

    def _add(num, src):
        k = _nk(num)
        if k:
            where.setdefault(k, {"display": str(num).strip(), "sources": set()})
            where[k]["sources"].add(src)

    for x in (out["ncbi"].get("type_material") or []):
        _add(x, "ncbi")
    for x in (out["cbs"].get("type_numbers") or []):
        _add(x["number"], "cbs")
    for x in (out["cgmcc"].get("rows") or []):
        _add(x["number"], "cgmcc")
    out["number_shared"] = {k: {"display": v["display"], "sources": sorted(v["sources"])}
                            for k, v in where.items() if len(v["sources"]) >= 2}

    # ---- ①b 号比对明细：NCBI 声明的类型号 vs 各机构自己的账 ----
    # 回答「四家的菌种号是否相同」：不是比数字大小，而是比【谁把哪个号算作本种的类型材料】。
    def _nums(lst):
        return [str(x).strip() for x in (lst or []) if str(x).strip()]

    ncbi_tm = _nums(out["ncbi"].get("type_material"))
    cbs_tm = [x["number"] for x in (out["cbs"].get("type_numbers") or [])]
    cg_all = [x["number"] for x in (out["cgmcc"].get("rows") or [])]
    cg_typ = [x["number"] for x in (out["cgmcc"].get("type_rows") or [])]

    def _same_code(nums, code):
        return [n for n in nums if _nk(n).startswith(_nk(code)[:len(code)])]

    ncbi_cg = _same_code(ncbi_tm, "CGMCC")
    ncbi_cbs = _same_code(ncbi_tm, "CBS")
    out["number_check"] = {
        "ncbi_type_numbers": ncbi_tm,
        "cbs_type_numbers": cbs_tm,
        "cgmcc_all": cg_all,
        "cgmcc_type": cg_typ,
        # NCBI 说是 CGMCC 的模式号，但 CGMCC 该学名下没有这个号
        "ncbi_cgmcc_absent": [n for n in ncbi_cg if _nk(n) not in {_nk(x) for x in cg_all}],
        # NCBI 说是 CBS 号，但 CBS 目录按该名的类型声明里没有它（CBS 标注率不全，属弱证据）
        "ncbi_cbs_absent": [n for n in ncbi_cbs if _nk(n) not in {_nk(x) for x in cbs_tm}],
    }
    # 只进 notes：机构"没有"是空白，不是反证
    for n in out["number_check"]["ncbi_cgmcc_absent"]:
        out["notes"].append("NCBI 声明 {} 为本种模式号，CGMCC 该学名下未见此号"
                            "（空白·权重 1，不计错）".format(n))
    for n in out["number_check"]["ncbi_cbs_absent"]:
        out["notes"].append("NCBI 声明 {} 为本种类型材料，CBS 目录该名的类型声明里没有它"
                            "（CBS 仅 24.1% 条目有标注，弱证据，不计错）".format(n))

    # ---- ② 模式 vs 普通：只认【双方都正面声明】；空白一律不算 ----
    # 规则（用户定）：
    #   · CGMCC 的「空白」不计错（少一条不能当反证）
    #   · 但 CGMCC 的「正面声明」算数 —— 只要双方都给了明确说法，冲突就是错误，与权重无关
    #   · 冲突有两种形态：
    #       (a) 同一个号：一方说「模式」、另一方说「非模式」
    #       (b) 双方各自指定了【不同】的模式号
    claim = {}   # 归一化号 → [(来源, "模式"/"非模式", 原文)]

    def _claim(num, src, kind, raw=""):
        k = _nk(num)
        if k:
            claim.setdefault(k, []).append((src, kind, str(raw)[:60]))

    for x in (out["ncbi"].get("type_material") or []):
        _claim(x, "ncbi", "模式", "NCBI 类型材料")
    for x in (out["cbs"].get("type_numbers") or []):
        _claim(x["number"], "cbs", "模式", x["status"])
    for x in (out["cgmcc"].get("rows") or []):
        ist = x.get("is_typestrain") or ""
        if not ist:
            continue                      # 空标注 = 没说，不算声明
        _claim(x["number"], "cgmcc", "模式" if ("模式菌株" in ist and "非" not in ist) else "非模式", ist)

    # (a) 同一号上的正反冲突
    for k, lst in claim.items():
        kinds = {x[1] for x in lst}
        if len(kinds) < 2:
            continue
        msg = "模式/普通冲突 ｜ 号 {} ｜ ".format(lst[0][2] and k or k) + "；".join(
            "{} 说{}{}".format(s_, kd, "（{}）".format(r) if r else "") for s_, kd, r in lst)
        out["conflicts"].append(msg)
        out["type_conflict"].append({"kind": "同号正反", "number": k,
                                     "claims": lst, "level": "矛盾"})

    # (b) 双方各自指定了【不同】的模式号
    #     ⚠️ 这一档只能判「待核」，不能判错：
    #        一个种可以合法地同时有 holotype / ex-type / paratype；而且同一管菌在不同机构
    #        各有编号。实测 siamense：NCBI 列 BPDI2 / CBS 125378，CGMCC 标 3.14174，
    #        而 CBS 130417（真 holotype）只登记了 MFLU 090230 / ICMP 18578 两个别名
    #        —— **没有任何一方登记过 CGMCC 3.14174 与它们的对应关系**，所以无从判断
    #        它到底是不是同一管菌。既不能算错，也不该当作没问题。
    #     映射豁免：若能通过 CBS 的其他机构号把两边对上，则不算冲突。
    ncbi_typ_nums = [k for k, v in claim.items()
                     if any(src == "ncbi" and kd == "模式" for src, kd, _ in v)]
    cg_typ_nums = [k for k, v in claim.items()
                   if any(src == "cgmcc" and kd == "模式" for src, kd, _ in v)]
    if ncbi_typ_nums and cg_typ_nums and not (set(ncbi_typ_nums) & set(cg_typ_nums)):
        mp = sources.get("cbsmap") if sources else None
        linked = False
        if mp is not None and hasattr(mp, "code_map"):
            for cn in cg_typ_nums:
                try:
                    r = mp.code_map("CGMCC", _orig_num(out, cn))
                    hits = {_nk(m.get("display")) for m in ((r or {}).get("matches") or [])}
                    if hits & set(ncbi_typ_nums):
                        linked = True
                        break
                except Exception:
                    pass
        if linked:
            out["notes"].append("双方模式号不同，但可经 CBS 其他机构号对上 → 同一管菌，不计错")
        else:
            out["pending"].append(
                "模式号归属待核 ｜ NCBI 声明为类型材料的是 {} ｜ CGMCC 标为模式菌株的是 {} "
                "｜ 两者无交集，且无任何一方登记过对应关系（需查原始发表/机构 strain history）".format(
                    [x for x in (out["ncbi"].get("type_material") or []) if _nk(x) in ncbi_typ_nums][:3],
                    [x["number"] for x in (out["cgmcc"].get("type_rows") or [])][:3]))
            out["type_conflict"].append({"kind": "不同号各自称模式（待核）",
                                         "ncbi": ncbi_typ_nums, "cgmcc": cg_typ_nums,
                                         "level": "待核"})

    # ---- ③ 异名 vs 独立种（NCBI ↔ MycoBank）----
    try:
        cur = out["myco"].get("current") or ""
        if cur and out["ncbi"].get("ok"):
            w = str(cur).split()
            accept = "{} {}".format(w[0], w[1]) if len(w) >= 2 else str(cur)
            same = accept.split()[-1].lower() == name.split()[-1].lower()
            # 拼写变体先于「分歧」判定：parvae/parva、tuyouyou/eucommiae 这类，
            # 种加词确实不同，但 MycoBank 已明确标 Orthographic variant —— 同一个名字的
            # 两种写法，不是两个物种。判成「分歧」会让 127 条分歧里混进 7 条假分歧。
            _st = (out["myco"].get("status") or "").lower()
            _is_ortho = _st.startswith("orthographic")
            # 词尾变体：同一个种加词、因属名性别不同而词尾不同
            #   Trichoderma（中性）latizonatum ↔ Hypocrea（阴性）latizonata
            # 这是「一菌一名」后两世代属名并存的产物，不是两个物种。
            # 判据：去掉末尾 0–2 个字母后词干完全相同（词干至少 5 个字母）。
            def _epi_variant(x, y):
                for ca in (0, 1, 2):
                    for cb in (0, 1, 2):
                        sa = x[:len(x) - ca] if ca else x
                        sb = y[:len(y) - cb] if cb else y
                        if len(sa) >= 5 and sa == sb:
                            return True
                return False
            _is_ending = (not same) and _epi_variant(name.split()[-1].lower(),
                                                     accept.split()[-1].lower())
            nb = out["name_conflict"] = {
                "ncbi": "独立 taxid({})".format(out["ncbi"]["taxid"]),
                "mycobank": ("异名 → {}".format(accept.split()[-1]) if same
                             else "独立（{}）".format(accept.split()[-1])),
                "status": out["myco"].get("status", ""),
                "verdict": ("一致" if same
                            else ("拼写变体" if _is_ortho
                                  else ("词尾变体" if _is_ending else "分歧")))}
            if not same and not _is_ortho and not _is_ending:
                lvl = "矛盾" if min(W["ncbi"], W["myco"]) >= 2 else "提示"
                out["name_conflict"]["level"] = lvl
                note = ("NCBI 把 {} 当独立种（taxid {}），MycoBank 把它的接受名判为 {}"
                        .format(name, out["ncbi"]["taxid"], accept))
                (out["conflicts"] if lvl == "矛盾" else out["notes"]).append(note)
    except Exception:
        pass

    # ---- 空白：只进 notes，永不计错 ----
    if out["cbs"].get("ok") and not out["cbs"].get("type_numbers"):
        out["notes"].append("CBS 目录没有该名的类型声明（空白，不计错）")
    if out["cgmcc"].get("ok"):
        if not (out["cgmcc"].get("rows") or []):
            out["notes"].append("CGMCC 未收录该学名（空白，权重 1，不计错）")
        elif not (out["cgmcc"].get("type_rows") or []):
            out["notes"].append("CGMCC 有该种但未标模式菌株（空白，不计错）")
    return out


# ---------------------------------------------------------------- 跨机构号连通分组
def link_numbers(numbers, sources, icmp_budget=200, progress=None):
    """把一批菌株号按【同一管实物】聚成组（并查集）。**两座桥，缺一不可**。

    桥一 · CBS.other_coll    CBS 号 → 别家号（覆盖广，但会漏登：实测 CBS 125378 无别名）
    桥二 · ICMP.Duplicates   ICMP 号 → 别家号（覆盖窄，但能补上面那个洞）

    输入其它机构的号时，先用 CBSMap.code_map 反查它对应的 CBS 号，再并入。

    ⚠️ 三条边界（必须随结果一起交代）：
        1. 只有两家登记过对应关系的号才会被连起来；**没连上不等于不是同一管菌**，
           应视为「独立组」，不能当成「不同的菌」。
        2. 非 CBS / 非 ICMP 的号（如 CGMCC、MycoBank）没有自己的桥，
           只能靠上面两座桥被间接连上。
        3. ICMP 抓取是逐个详情页（1.5 秒/次），用 icmp_budget 封顶，
           **被跳过的会记在 stats.icmp_skipped**，不能当作「查过没有」。

    返回：
        {"groups": [[号...]] ,      # 只含 ≥2 个成员的连通组
         "singletons": [号...],     # 没连上任何别的号
         "edges": [{"a","b","via"}],# 每条连接及来源（用于审计）
         "stats": {...}}
    """
    W = [{"v": None}]

    def find(x):
        while W[x["v"]] if False else x["v"] != x["v"]:
            pass
        return x
    # —— 用简单字典实现并查集（避免上面那种易错写法）——
    parent = {}

    def find_(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union_(a, b):
        ra, rb = find_(a), find_(b)
        if ra != rb:
            parent[rb] = ra

    def nk(v):
        r"""归一化一个号。**号码部分必须含数字**，否则返回 None。

        ⚠️ 为什么必须卡这一条：正则的号码段是 [\w.\-]+，**字母也能匹配**，
           于是裸代码 "ICMP" 会被拆成 代码=ICM / 号码=P，归一化后仍是 "ICMP"，
           变成一个「看起来合法的脏号」，污染连通分组（实测出现过 CBS102189 = ICMP = ICMP13903）。
        """
        m = re.match(r"^\s*([A-Za-z][A-Za-z\.\-]{0,12}(?:<[^>]*>)?)\s*[:#]?\s*([\w\.\-]+)\s*$",
                     str(v or ""))
        if not m or not re.search(r"\d", m.group(2)):
            return None
        return norm_key(m.group(1), m.group(2))

    edges = []
    stats = {"cbs_queries": 0, "icmp_queries": 0, "icmp_skipped": 0,
             "code_map_queries": 0, "unknown_format": 0}
    cbsmap = sources.get("cbsmap") if sources else None
    icmp = sources.get("icmp") if sources else None

    seen = []
    for raw in numbers:
        k = nk(raw)
        if not k or not re.search(r"\d", k):
            stats["unknown_format"] += 1     # 裸机构代码等脏号：记为格式异常，不进图
            continue
        find_(k)
        if k not in seen:
            seen.append(k)

    for k in list(seen):
        m = None
        for raw in numbers:
            if nk(raw) == k:
                m = re.match(r"^\s*([A-Za-z][A-Za-z\.\-]{0,12}(?:<[^>]*>)?)\s*[:#]?\s*([\w\.\-]+)\s*$",
                             str(raw or ""))
                break
        if not m:
            continue
        code, num = m.group(1).upper(), m.group(2)
        # 桥一：CBS → 别家
        if code == "CBS" and cbsmap is not None and hasattr(cbsmap, "numbers_of"):
            stats["cbs_queries"] += 1
            try:
                for a in ((cbsmap.numbers_of(num) or {}).get("aliases") or []):
                    t = nk(a.get("display"))
                    if t and t != k:
                        union_(k, t); edges.append({"a": k, "b": t, "via": "CBS.other_coll"})
            except Exception:
                pass
        # 桥二：ICMP → 别家
        elif code == "ICMP" and icmp is not None and hasattr(icmp, "duplicates_of"):
            if stats["icmp_queries"] >= icmp_budget:
                stats["icmp_skipped"] += 1
                continue
            stats["icmp_queries"] += 1
            try:
                for d in (icmp.duplicates_of("ICMP", num) or []):
                    t = nk(d.get("number"))
                    if t and t != k:
                        union_(k, t); edges.append({"a": k, "b": t, "via": "ICMP.Duplicates"})
            except Exception:
                pass
            if progress:
                progress(stats["icmp_queries"])
        # 其它机构：反查 CBS 号
        elif cbsmap is not None and hasattr(cbsmap, "code_map"):
            stats["code_map_queries"] += 1
            try:
                r = cbsmap.code_map(code, num)
                for mm in ((r or {}).get("matches") or []):
                    t = nk(mm.get("display"))
                    if t and t != k:
                        union_(k, t); edges.append({"a": k, "b": t, "via": "CBSMap.code_map"})
            except Exception:
                pass

    comp = {}
    for k in list(parent.keys()):
        comp.setdefault(find_(k), []).append(k)
    groups = sorted([sorted(v) for v in comp.values() if len(v) >= 2], key=len, reverse=True)
    single = sorted([v[0] for v in comp.values() if len(v) == 1])
    return {"groups": groups, "singletons": single, "edges": edges, "stats": stats}


def inherit_from_groups(groups, judgements):
    """从连通组里推出「可借用的权威」——把邻居已核实的结论传给未核的号。

    参数
        groups     : link_numbers() 的 groups（每项是归一化号的列表）
        judgements : {归一化号: {"flags": [...], "display": ..., "source": ...}}

    返回 {归一化号: {"from": 锚点展示号, "flags": [...], "via": 连接来源}}

    ⚠️ 三条保守规则（写死，别放宽）：
        1. **只借「已核过」的锚点**：锚点必须来自 CBS 目录（唯一有公开机读目录的那家），
           或者本身就是已判定的 CBS 号。
        2. **绝不继承「号不存在」**：那等于替一个没核过的号断言它不存在 —— 是最危险的一类外推。
        3. 结果一律标成「同组映射（借自 …，未确认）」，**不能当成已核实**。
    """
    out = {}
    for g in groups or []:
        anchors = []
        for k in g:
            j = (judgements or {}).get(k) or {}
            fl = j.get("flags") or []
            disp = j.get("display") or k
            # 锚点必须【干净】—— 用白名单，不用黑名单。
            # ⚠️ 为什么必须是白名单：实测黑名单写法会漏掉「号不存在」，
            #    于是从"这个号不存在"的邻居那里借信用，等于替未核的号造结论。
            if not str(disp).upper().startswith("CBS"):
                continue
            SAFE = {"异名（非错）", "基名属", ""}
            if any(str(x).strip() not in SAFE for x in fl):
                continue
            anchors.append((k, disp, fl))
        if not anchors:
            continue
        for k in g:
            j = (judgements or {}).get(k) or {}
            fl = j.get("flags") or []
            # 只借这两档：
            #   无目录 —— 该机构没有公开目录，靠邻居拿到它是哪管菌、什么模式地位
            #   未核   —— 号存在但目录没标模式地位，邻居的「Holotype of X」正好补上
            # 别的一律不动（尤其「号不存在」：那是断言，绝不能从邻居外推）
            if not ({"无目录", "未核"} & set(fl)):
                continue
            a = anchors[0]
            out[k] = {"from": a[1], "via": "同组",
                      "flags": ["同组映射（借自 {}，未确认）".format(a[1])]}
    return out
