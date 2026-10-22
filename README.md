
## 两种用法：

### 路线一（不想碰命令行、只想看结论）：双击 exe

**下载整个仓库**（GitHub 页面 `Code ▾` → `Download ZIP`）→ 解压 → 双击文件夹里的
**`类型校对台.exe`** → 浏览器自动打开 <http://127.0.0.1:8765>

打开后**不需要联网、不需要跑流水线**，两步看到结论：

1. 顶部「**三属审计结论速览**」——规模、发现构成、不一致率、分歧裁决、验证子样本，一屏看完
2. 左栏「**不想跑流水线？直接载入仓库里已有的结果**」→ 选 Colletotrichum / Trichoderma / Diaporthe
   → 切到审稿区，逐条查看每个种、每个号的判定与证据

- 不需要安装 Python、不需要装任何依赖
- 它是**本机**网页服务：进程关掉服务即停，数据不出本机
- Windows 若提示「未知发布者」：更多信息 → 仍要运行

> ⚠️ **exe 只是启动器，本身不含代码**（14 MB，内含 Python 运行时与 Biopython）。
> 必须与仓库里的 `src/`、`tools/` 处在同一目录树内。
> 单独拿走 exe 双击时，程序会弹窗说明，并在它旁边生成 `repo.txt` 模板——
> 填入仓库路径（或设环境变量 `TYPEPAUDIT_REPO`）即可。

### 路线二（想改代码 / 想审计实现）：命令行

需要 **Python 3.8 或更高**；核心功能**零第三方依赖**（不装任何包也能跑）。

```bash
git clone https://github.com/xiangccookie/typeaudit.git
cd typeaudit

# 网页界面（然后浏览器打开 http://127.0.0.1:8765）
python tools/pipeline_web.py

```

想知道装完是否正常，跑一次自检（不联网、不需要数据文件）：

```bash
python tests/smoke_test.py
```

## 只看结论 → 不用准备任何数据

仓库自带 `results/`（三属结果表，约 1.8 MB）。双击 exe → 顶部「结论速览」+
左栏「载入已含结果」，即可看到全部结论，**不联网、不装依赖、不用准备第三方数据**。
下一节只在你**要自己重跑**时才需要。

## 改代码前须知（会咬人的一条）

**任何改动 `tools/pipeline_web.py` 之后，先跑一次自检再启动：**

```bash
python tools/check_web.py
```

它做三件事：① 内联 JS 语法（node --check）② 重复声明 ③ 元素 id 对账。
**第一项不通过 → 页面上【所有按钮都会点不动】**（整个 `<script>` 被浏览器丢弃），
而现象只是「点了没反应」，很容易误判成后端问题。

> 本项目真实踩过：插入新函数时把原有的 `async function setTab` 改成了 `function setTab`，
> 导致重复 `async` → 语法错误 → 开始按钮失效。自检脚本当时就写着这条后果，但没跑。

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

不解码就会与对账表逐字对不上。**凡从 XML 抽文本，一律先解实体**（`&lt; &gt; &amp; &quot;`
以及数值实体 `&#xNNNN;`）。

## 踩过的坑

- **`docs/bug判别特征手册.md`** —— 本项目实测的 bug 判别特征（症状 → 判别特征 → 根因 → 正解）。
  新写代码前扫一眼，尤其是「顺序敏感」和「静默失败」两节。
