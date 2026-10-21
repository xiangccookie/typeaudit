# typeaudit · 分类学类型材料元数据审计工具箱

从公开数据库与分类学文献里，审计「类型材料保藏号」这一类元数据是否自洽。
每个模块都能独立使用，**纯 Python，不依赖任何本项目数据**。

## 两种用法：选一条

### 路线一（不想碰命令行）：下载免安装 exe

到本仓库 **Releases** 页面下载 **`类型校对台.exe`** → 双击 → 浏览器打开 <http://127.0.0.1:8765>

- 不需要安装 Python、不需要装任何依赖
- 它是**本机**网页服务：窗口关掉服务即停，数据不出本机
- Windows 若提示「未知发布者」：更多信息 → 仍要运行

### 路线二（想改代码 / 想审计实现）：命令行

需要 **Python 3.8 或更高**；核心功能**零第三方依赖**（不装任何包也能跑）。

```bash
git clone https://github.com/xiangccookie/typeaudit.git
cd typeaudit

# 网页界面（然后浏览器打开 http://127.0.0.1:8765）
python tools/pipeline_web.py

# 或直接命令行跑全属扫描（--limit 20 = 先跑前 20 个种看效果）
python -m typeaudit.pipeline --genus Colletotrichum --out out --limit 20
```

想知道装完是否正常，跑一次自检（不联网、不需要数据文件）：

```bash
python tests/smoke_test.py
```

## 跑之前：先准备外部数据（关键，别跳过）

仓库里**不带**第三方数据快照（CBS 目录、NCBI 分类表、MycoBank 导出）——体积大且涉及第三方许可。
不准备也能跑，但结果会缺少一部分源，表里会标「无公开数据库 / 未核」。
**那不是「没问题」，是「没查」**，请勿当成核查通过。

```bash
# 1) NCBI 分类表（联网）
python tools/setup_data.py --email 你的邮箱@example.com

# 2) CBS 保藏目录（先 --check 只问总数，确认通了再正式拉）
python tools/fetch_cbs_live.py --check
python tools/setup_data.py --fetch-cbs --genus Colletotrichum --email 你的邮箱@example.com

# 3) MycoBank：先在官网导出 MBlist 压缩包，再导入
python tools/import_mycobank.py MBlist.zip -o data/mycobank_full.tsv
```

产物默认落在 `data/`。这些文件属第三方数据、体积较大，**不入库**（`.gitignore` 已排除），
换机器时重新跑上面三步即可。

---

## 装上


```bash
git clone https://github.com/xiangccookie/typeaudit.git
cd typeaudit
pip install -e .            # 只用标准库
pip install -e ".[ncbi]"    # 需要联网抓 PMC 全文时（装 biopython）
```

跑一次自检（不联网、不需要任何数据文件）：

```bash
python tests/smoke_test.py
```

> 📖 **想审计整个属（如炭疽菌全属）**：见 [`docs/查炭疽菌全属_上手路径.md`](docs/查炭疽菌全属_上手路径.md)
> —— 六个任务的入口选择、参数含义、外部数据清单、以及「这套工具做不到什么」。


## 流水线：给一个属名，出三列表

```bash
# 联网：自动取该属的物种清单 + 各源声明，跑一遍
python -m typeaudit.pipeline --genus Colletotrichum --out out --limit 20

# 投喂框架论文（PDF 抽出的 txt）→ 额外产出框架对照表（J-01 四态）
python -m typeaudit.pipeline --genus Colletotrichum --framework SiM2022.txt --out out

# 离线：只给本地语料目录取种名
python -m typeaudit.pipeline --genus Colletotrichum --corpus local/ --out out
```

（未 `pip install -e .` 时加 `PYTHONPATH=src`。）

**主产出 `<genus>_三列表.tsv`**：

| 列 | 语义 | 取数 |
|---|---|---|
| 种名 | 该属的物种 | 框架论文 / 本地语料 / GenBank 检索 |
| **菌株号（含问题）** | 序列库记录里声明的号 | `GenBankAPI.culture_collections_of()`（保藏号，主）+ `strain_values_of()`（论文内部编号，旁证） |
| **模式菌（含问题）** | 类型材料字段里声明的号 | `TaxonomyAPI.type_material_of()` + `GenBankAPI.type_material_of()` |

**两列为什么必须分开**：菌株号是"序列库随手记的"，模式菌是"被声明为类型的"。
错误恰恰长在两者的落差上——本项目最值钱的一条错值（`CBS 130417` 漏登 / `CBS 125378` 混入）
就是"接受名的类型材料槽位被异名的模式号占据"，**只有并排才看得见**。

**问题标记**（grep 文字，勿改 emoji）：`号不存在` / `目录记为别种` / `属异名` / `一株多名` /
`位数孤例` / `仅数据库有` / `仅文献有` / `未核`。

另出 `<genus>_framework.tsv`（投喂框架论文时）与 `<genus>_metrics.json`（含 `conventions` 口径）。

## 五个模块

| 模块 | 作用 | 依赖 |
|---|---|---|
| `typeaudit.designations` | 培养物型：抽 `Typus:` 段与等价保藏号链（`A = B = C`） | 无 |
| `typeaudit.specimen` | 标本型：抽「标本馆代码 + 编号」并与 GenBank voucher 对账 | 无 |
| `typeaudit.gate` | 换属前的可行性门（三根轴） | 无（联网需 biopython） |
| `typeaudit.numbers` | 核数器：文档声明数字 ↔ 源表重算（规则外部化） | 无 |
| `typeaudit.gaps` | 缺环分类器：挂账项 → 可交接任务单（配置外部化） | 无 |

## 命令行

```bash
python -m typeaudit.gate --file paper.xml --genus Diaporthe
python -m typeaudit.gate --genus Diaporthe --batch --email you@example.org   # 联网
python -m typeaudit.numbers --rules rules/example_rules.json --report out/report.tsv
python -m typeaudit.gaps --config config/gaps_example.json --data-dir examples/data --out out/gaps.tsv
```

（未安装时用 `PYTHONPATH=src` 即可运行。）

## 两条来自实测的取数纪律

**其一，分类学清单走 Taxonomy，不要用 GenBank 标题抠。**
为拿「属里有哪些种」这个索引信息去 `efetch` 数百条 GenBank 记录（2–3 MB），
成本比 Taxonomy 路径高**两个数量级**（实测：旧路径分钟级、Taxonomy 秒级）。
正确做法三步，一步都不能省：

```python
tx.species_in_genus("Colletotrichum")     # 已封装
# ① 检索用 `{genus}[Subtree]`——**不要加 `AND species[Rank]`**，实测不生效
#    （Colletotrichum 返回 2277，把亚种/环境样本/未培养全算进来了）
# ② efetch 后本地按 `<Rank>species</Rank>` 过滤
# ③ **先剥 `<LineageEx>` 再取 `<ScientificName>`**——否则会把
#    cellular organisms / Eukaryota / Fungi 当成种名（实测第一版就这样）
```
（`cf.` / `aff.` / `sp.` 是疑似鉴定标记，非正式种名，已一并过滤。）

**其二，拿机构目录核号之前，先问「这个号属于哪个机构」。**
用 CBS 目录去核 `ICMP 18581`、`MFLU 090228` 这类号，会**一律判成「号不存在」**——
本项目实测第一版因此产出 **40 个纯假阳性**。**没有对应的机读目录就只能标「无目录」，
不能标「号不存在」。** 另有三个易错点：`CBS H-xxxx` 是**标本馆号**（不走培养物目录）；
目录里的基名属（`Vermicularia` → *Colletotrichum*）**不是错**；
`norm_key` **必须去点**（目录规则：`CBS 1546.1` → `CBS15461`，保留点会命中不了）。

## 三条来自实测的硬约束

这些不是风格偏好，是**实测翻车后固化下来的**：

**其一，切块规则不得预设作者名格式。**
要求「作者名以大写字母+点开头」（`F. Liu`）会整类漏掉完整姓氏写法
（`Diaporthe afzeliae Monkai and S. Lumyong` → 一块都切不出来）；
重组合名 `Xylaria phyllocharis (Mont.) Fr.` 的种加词后接的是 `(` 而不是大写字母，
也需一并覆盖。

**其二，标签变体表 + 命中数自检。**
类型指定写作拉丁 `Typus:`，不是英文 `Type:`。只搜 `Type:` 会得 0 处，
从而**误判「该文无类型指定」**——这是假阴性。故必须打印每个标签的命中数，
让 0 命中可见。

**其三，方法有适用边界，且边界要量化。**
等价链法**只适用于培养物型**类群：实测 166 条目 / 78 属，
培养物型抽出率 **63%**，标本型（正文只写 `Holotype`）抽出率 **0%**。
这是硬边界，不是调参问题——标本型需要馆藏号正规化层，或改用 GBIF 等第三方源。

## 已知的负结果（记录在案，避免重复踩）

- **标本型对账产出率为 0**：12 个候选种里仅 4 个在 GenBank 有序列，
  模式标本号被独立印证 **0 个**。三层根因：① 多数标本型新种尚无序列；
  ② 测序标本 ≠ 模式标本；③ `/specimen_voucher` 无受控词表
  （同一批里出现 `HAST 147515` / `528 (HAST, JF)` / `LIP:Bon 90132` 等写法）。
- **标题检索会命中同名异属乃至异界**：`"Bipolaris"[Title]` 命中一篇
  《Two New Chalcid Wasps…》——`Ophelimus bipolaris` 是**膜翅目小蜂**的种加词，
  正文里 `Bipolaris` 出现 0 次。故必须校验属名在正文的出现次数。

## 需要联网时

NCBI E-utilities 要求提供邮箱，且有限流（建议 ≥1 秒/请求）。
用 `--email` 或环境变量 `NCBI_EMAIL` 提供自己的邮箱；抓取结果按
`(接口, 参数)` 落盘缓存，重跑零请求。

## 验证记录：洁净环境能跑（2026-10-21）

目的：证明「别人拿到代码、在没有本机配置的机器上也能跑」，而不是只有开发机可用。

### 做法

1. 新建隔离环境（不带开发机 site-packages）：
   ```bash
   python -m venv <临时目录>
   ```
2. **清空 UTF-8 相关环境变量**（模拟中文 Windows 默认的 GBK 控制台）：
   ```bat
   set "PYTHONUTF8=" & set "PYTHONIOENCODING="
   ```
3. 用该环境的解释器跑自检：
   ```bash
   <venv>\Scripts\python.exe tests\smoke_test.py
   ```

### 结果

| 检查项 | 结果 |
|---|---|
| `tests/smoke_test.py` | **47 项全部 OK，失败项 = 0，退出码 0，无 Traceback** |
| `tools/check_known_errors.py --no-run` | 正常输出含 emoji 的报告，无 `UnicodeEncodeError`（退出码 2 是它自身的语料覆盖率语义，不是崩溃） |
| `python -m typeaudit.pipeline --help` | 正常打印帮助 |

### 这次验证抓出并修掉的一个真 bug

**症状**（修之前，在洁净环境下必现）：

```
UnicodeEncodeError: 'gbk' codec can't encode character '\u2194' in position 19: illegal multibyte sequence
```

**根因**：源码里大量用 `↔ ⇒ ⚠️ 🔴 🟠 🔵 ✅ ❌` 作为输出符号，而中文 Windows 的控制台默认是
GBK(cp936)。开发机设了 UTF-8 环境变量，所以**本机永远看不出这个问题**；换一台干净机器直接崩。
更要紧的是 `rules` / `reconcile` 的判定结果本身就用 🔴🟠🔵 当标记 —— 等于主功能也会崩。

**修法**：

- 新增 `src/typeaudit/_console.py`：把 `stdout`/`stderr` `reconfigure(encoding="utf-8", errors="replace")`；
- `typeaudit/__init__.py` 导入即调用 → 覆盖所有 `import typeaudit` 的代码（主流水线 + 13 个工具）；
- `tools/` 下**不**导入本包的 19 个独立脚本，各自顶部插入同样的引导代码。

> **教训**：**「我本机能跑」不等于「别人能跑」**。
> 只要输出里有非 ASCII 符号，就必须在**清空 UTF-8 环境变量的隔离环境**里验一遍；
> 只在本机跑通就宣布可用，等于把一个必崩的版本发出去。

---

## 许可

MIT。见 `LICENSE`。

## 全自动管道（一条命令）

```bash
# 离线：分析一个本地语料目录
python -m typeaudit.audit --genus Colletotrichum --corpus examples/corpus --out out

# 联网：按检索阶梯抓 PMC 全文再分析（需 biopython 与自己的邮箱）
python -m typeaudit.audit --genus Colletotrichum --fetch --email you@example.org --out out
```

六个阶段：`acquire → gate → extract → rules → metrics → report`。
产物：`<genus>_species.tsv`（物种级记录）、`<genus>_findings.tsv`（J 规则候选）、
`<genus>_metrics.json`（指标 + **口径** + 规则覆盖报告）。

## J 规则库（21 条，`typeaudit.rules`）

《类型材料元数据审计规范 v1》的 J-01…J-21 已编码为可执行检查，每条标注自动化程度：

| 程度 | 条数 | 规则 |
|---|---|---|
| `offline` 纯本地 | **8** | J-01 J-06 J-07 J-08 J-13 J-14 J-15 J-20 |
| `needs_table` 需外部表 | 4 | J-02 J-09 J-17 J-18 |
| `needs_net` 需联网/适配器 | 9 | J-03 J-04 J-05 J-10 J-11 J-12 J-16 J-19 J-21 |

### ⚠️ 候选质量：本版机械实现**假阳性高**，只可用于**缩小候选集**

实测（12 篇炭疽菌论文，`J-01 J-06 J-07 J-08 J-13 J-14 J-15` 全跑）：

| 规则 | 候选 | 抽检结论 |
|---|---|---|
| J-07 一号两主 | 12 → 修后 **4** | 修前几乎全是**块污染**（块跨表格 ⇒ 同一个号出现在多个物种块；`Table1` 不是号） |
| J-08 位数孤例 | 14 → 修后 **12** | 几乎全是假阳性：`LC 0951` vs `LC 13850` 只是流水号，**LC 本不定宽**；真异常（`CGMCC 3.2051`）成立是因为 CGMCC 固定 5 位 |
| J-13 近失配 | 0 | 本语料无跨来源同物种对 |
| J-01 跨来源对账 | 0（1 条为 `agree`，按约定不入表） | —— |

**两条必须一起用的护栏**：

1. **J-07 只可喂「该物种 `Typus` 段内声明的号」**，喂块内全部号会因块跨表格产生大量假阳性。
2. **J-08 需要「同批同构」前提**：代码里加了门槛（同批 ≥5 个且众数位数占比 ≥70%），
   但要真正提高精度，**使用方应按保藏机构补一张「编号是否定宽」表**
   （CGMCC 的 `3.xxxxx` 定宽；LC 的流水号不定宽）。

**另有两条口径约定（写进代码注释与规范）**：

- **J-01 四态不可合并**：`agree` / `disagree` / `only_a` / `only_b`。
  「对方没给值」与「对方给了别的值」在集合运算里同形，但含义相反；
  把后者混进前者会把"不可检验"当"可疑"（本项目曾误算 31 行）。
- **J-15 是元规则，报「检出 0」之前必跑**：实测 12 篇炭疽菌论文里，
  「类型材料字段」枚举到 287 个号，而全文含 2690 个类保藏号记号
  ⇒ **2508 个号从未被检查**。规范里 `CBS 147945` 就是这么漏掉的。

## 口径（数字自解释）

报告里的每个指标都随 `_metrics.json` 带出 `conventions`，三处口径**声明式固定**：

| 口径 | 取值 | 为什么 |
|---|---|---|
| 切块变体 | `B_paren` | 种加词后接大写字母**或左括号**；只要求大写字母会漏掉 `Xylaria phyllocharis (Mont.) Fr.` |
| 标签集 | `full(13)` | 全变体表；窄集会低估含类型指定块数 |
| 链计数 | `all` | 统计块内**全部**链；只取首链会低估多号链数 |

实测差异（同一批 12 篇炭疽菌论文）：**链计数 `all` = 105 条多号链，只取首链 = 54 条**。
口径不写死，「同一个数字」的含义就会随时间漂移——本项目曾因此漏改 3 处。

## 复现评估：本工具链能否自动检出已确认的 10 条错值

`examples/reproduce_colletotrichum.py` 是端到端检验（不是单元测试）：每条错值给出检测配方，
实跑后报告检出与否。

```bash
python examples/reproduce_colletotrichum.py \
  --cbs data/cbs_catalogue.tsv \
  --registry data/NCBI_biocollections_full.csv \
  --sim data/SiM2022.xml --protologue 6422929 \
  --live --cache .cache
```

**实测结果：8 / 10 可自动检出。**
> **口径**：合并错值条目时，「检出数」要同步减 —— E4/E5 原先各计 1 个检出点（都用 J-05），合并后只剩 1 个，故 9 → **8**；而未检出的 E7、E9 与合并无关，**仍为 2 条**。

| 错值 | 检出规则 | 用到的源 |
|---|---|---|
| E1 `CBS 147945` | J-04（目录无标注）+ J-03（属分布 `{Glomerella:5, Penicillium:1}`） | CBS 目录（离线）+ GenBank |
| E2 `CGMCC 3.17361` | J-16（官网记「非模式菌株」） | CGMCC 官网 API |
| E3 `CGMG 3.14187` | J-02（前缀不在登记表内） | Biocollections 登记表 |
| E4 `CBS 130417` 漏登 ／ `CBS 125378` 异名混入（**同一槽位两面**，原 E4+E5 已合并） | J-05（按学名反查目录） | CBS 目录 |
| E6 `CBS 132551` | J-04（目录记 Holotype of *Chrysosporium magnasporum*） | CBS 目录 |
| E8 `CGMCC 3.2051` | J-13（与正确值互为前缀） | 两侧数字 |
| E10 `BPI 892771` | J-03 + J-13 | GenBank |
| E11 `CBS 304.67` | J-20（`/note` 的 `Type status:`） | GenBank |

**未检出的 2 条及原因（都已定位，不是谜）**：

- **E9**（同源内部矛盾：图注 vs 正文 Typus）：已用 `<p>` 段落结构验证**可判**，
  但尚未并入正式流程。根因：压平文本上做正则切块会「取块取错」。
- **E7**（一号两主）：两处 `GMBC0209` 的模式材料**不在 `<p>` 段落里**（疑在表格/列表结构），
  段落法与文本法都取不到。⇒ 需要按 XML 的 table/list 结构再取一层。

## ⚠️ 一条重要的工程结论

**只要语料是结构化全文（PMC XML / JATS），就走「文档序单元」路径，不要在压平文本上做正则切块。**

```python
from typeaudit.textutil import segment
bl = segment(raw_xml, "Colletotrichum")     # 自动判别：XML 走单元，纯 txt 走 blocks
```

单元种类（按文档序）：`<title>` / `<p>` / `<caption>` / `<td>` / `<list-item>`。

**⚠️ `<title>` 绝不能漏。** JATS 的**种级标题**就在它里面：

```xml
<sec><title>Colletotrichum jishouense</title><p>…</p></sec>
```

漏掉 `<title>` 会让标题从未被看到，模式指定段因此被归到**上一个种**名下。
（另注：原文常把属名标成 `<italic>Colletotrichum</italic> jishouense`，
所以在原始 XML 里直接搜字面串 `"Colletotrichum jishouense"` **搜不到**——必须先按单元去标签。）

**实测收益**（同一批 13 篇炭疽菌论文）：

| 指标 | 压平文本 + 正则切块 | 文档序单元 |
|---|---|---|
| J-07 一号两主候选 | 12 条（几乎全假：块跨表格） | **2 条** |
| E7 的物种归属 | 错（报 `agaves`/`jishouense`） | **对**（`jishouense`/`tongrenense`） |

`blocks()` 保留，但**只作纯文本语料（PDF 抽出的 txt）的 fallback**。

## 一条工程教训：不要用字符串切片删代码

我在删除一份旧的 `coverage_report` 时，用「从 A 的起点切到 B 的起点」的方式删除，
结果**一起切掉了夹在中间的 4 个函数**（`units` / `is_heading` / `species_units` / `segment`），
`gate.py` 立刻 ImportError，冒烟测试报了 4 次才逐个找回来。

**正确做法**：用 AST 定位要删的函数边界，或至少先打印将被切掉的内容再动刀。
用正则/切片做结构性删除，边界一定按「函数定义到你确认的下一个函数定义」来算。

## 另一条工程教训：XML 实体必须解码

对账时有 6 行"值级不一致"，追查后**全部**是 XML 实体未解码：

| 抽到的 | 应为 |
|---|---|
| `BCC&lt;THA&gt;38876` | `BCC<THA>38876` |
| `MFLU&lt;THA&gt;LC0600` | `MFLU<THA>LC0600` |
| `CMM&lt;BRA&gt;4083` | `CMM<BRA>4083` |

带 `<馆代码>` 的复合号（如泰国 BCC 的 `BCC<THA>`）在 XML 里写成实体形式，
不解码就会与对账表逐字对不上。**凡从 XML 抽文本，一律先解实体**（`&lt; &gt; &amp; &quot;`
以及数值实体 `&#xNNNN;`）。

## 踩过的坑

- **`docs/bug判别特征手册.md`** —— 本项目实测的 bug 判别特征（症状 → 判别特征 → 根因 → 正解）。
  新写代码前扫一眼，尤其是「顺序敏感」和「静默失败」两节。
