# MycoBank MBList 增强：同名异物候选清单

> 日期：2026-10-23 ｜ 用途：把「待消歧」从「只知道要人判断」升级为「列出候选，人只需指认」
> **统计口径未改动**：分析宇宙仍是 2026-10-10 快照的 426,834 行 / 406,521 个学名。

---

## 一、文件来源与内容

`C:\Users\Administrator\Desktop\MBList.zip`（144 MB）→ 内含 `MBList.xlsx`（146 MB，544,942 行 × 13 列）。

**13 列**：
```
ID ｜ Taxon name ｜ Authors ｜ Authors (abbreviated) ｜ Rank ｜ Year of effective publication
｜ Name status ｜ MyCoBank # ｜ Hyperlink to MB ｜ Classification ｜ Current MycoBank #
｜ Current name ｜ Synonymy
```

**与旧库的匹配**（实测）：

```
旧库 426,818 个 MyCoBank #  →  新导出中 100% 都能对上
新导出 544,942 行（比旧库多 28%，为更新的一版）
95.3% 的行有 Authors
```

因此本目录只做**信息增强**，不替换分析宇宙——论文里 426,834 行 / 4.3% 这些数字**不受影响**。

---

## 二、实测数字（本研究口径）

| 项 | 数量 |
|---|---|
| 旧库行数 / 学名数 | 426,834 / 406,521 |
| **多行学名** | **17,404（4.3%）** |
| 为这些名字取到的候选行 | **37,718 行** |
| 2 个候选的名字 | 15,268（87.7%） |
| 3 个候选 | 1,642 |
| ≥4 个候选 | 497 |
| **只有 1 行 Legitimate** | **12,825（73.7%）** |
| 无 Legitimate 行 | 916 |

**含义**：73.7% 的多行名只需在「其余非法行」中确认，人工成本从"全文检索"降到"二选一"。

---

## 三、本研究 4 例的候选（可直接填进论文 2.3）

| 学名 | 候选 | MyCoBank # | Authors | Year | Status | 接受名 |
|---|---|---|---|---|---|---|
| *D. biguttulata* | 1 | 827907 | Y.H. Gao & L. Cai | 2014 | Legitimate | *D. eres* |
| | 2 | 810579 | **F. Huang, K.D. Hyde & Hong** | **2015** | Illegitimate | *D. pseudobiguttulata* |
| *D. longispora* | 1 | 802938 | **(Wehm.) R.R. Gomes et al.** | **2013** | Legitimate | *D. longispora*（本身） |
| | 2 | 827908 | Y.H. Gao & L. Cai | 2014 | Illegitimate | *D. eres* |
| *D. salinicola* | 1 | 250222 | **Speg.** | **1912** | Legitimate | *D. salinicola*（本身） |
| | 2 | 556589 | Dayar. | 2020 | Illegitimate | *D. eucommiae* |

**这解释了当初为什么判错**：机械规则「取 Legitimate 那行」在 *biguttulata* 上取到 Gao & Cai 2014 的 *eres*，
而 NCBI 指的是 F. Huang 2015 那个；*longispora* 则被取到 2014 那条（Illegitimate），而正确行是 2013 的本身。

### 关于 *D. glabrae* 的更正

我最初误以为 *D. glabrae* 是第 4 个同名异物。**实际不是**——它是 *D. chinensis* 指向的接受名。

真正的第 4 个同名异物是 ***D. chinensis***，在新导出中有 2 行：

| MyCoBank # | Authors | Year | Status | 接受名 |
|---|---|---|---|---|
| 570777 | H. Dong, J.W. Xia & X.G. Zhang | 2021 | Legitimate | *D. chinensis*（本身） |
| 556135 | **N.I. de Silva, Lumyong & K.D. Hyde** | 2021 | Illegitimate | **D. glabrae** |

**结论：4 条同名异物（biguttulata / longispora / salinicola / chinensis）全部成立，数量不变。**

---

## 四、为什么仍然不能全自动

**缺的是 NCBI 那一侧的命名人。**

抽 20 个种查 NCBI Taxonomy 的 `<Authority>` 字段：

```
有 authority 的: 0/20 = 0%
```

NCBI 对这些种根本没登记命名人，因此无法「拿 NCBI 的作者去匹配 MycoBank 的候选行」。
**能自动的是"列出候选"，不能自动的是"指认哪一行"。**

---

## 五、对论文 2.3 的建议

**现在写的是**：
> 本地导出仅 6 列、**无作者列**，故同名异物无法自动消歧，本研究的 4 条为人工核对所得。

**建议改为**：
> MycoBank 的完整导出含 Authors 与 Year of effective publication，可据此为多行学名自动列出候选
> （本研究涉及 17,404 个多行学名、73.7% 只有 1 行合法）；但 NCBI Taxonomy 对这些种未登记命名人，
> **指认哪一行仍须人工判断**，本研究涉及的 4 条即如此。

并补一句**导出日期**（本次导出为 2026-10 版），否则"是否同名异物"会随导出时点变化（如上述 *D. glabrae*）。

---

## 六、产物

| 文件 | 内容 |
|---|---|
| `mycobank_authors.tsv` | 全部 544,942 行 × 8 列（MB# / 学名 / 作者 / 作者缩写 / 年份 / 状态 / 接受名 / 分类） |
| `同名异物_候选清单_17404.tsv` | 本研究口径的 17,404 个多行学名 × 37,718 行候选（每个名字按年份排序，标注是否在快照内） |
| `同名异物_本研究4例.tsv` | 上述 3 个名字（*glabrae* 已非多行）的候选明细 |
| 本文件 | 说明与建议 |
