# 三张图的 Origin 作图说明

数据文件与本说明在同一目录。把 `.tsv` 直接拖进 Origin 即可（Origin 认 TSV，中文表头正常）。

---

## 通用设置（三张图都先设这个）

**① 期刊风格**
- 字体：**Arial**，坐标轴标签 **8 pt**、刻度数字 **7 pt**、图例 **7 pt**
- 坐标轴线宽 **0.75 pt**；刻度**朝内**（Tick Labels → Bottom/Left → In）
- 去掉图例边框与背景（Legend → Properties → Frame: None, Background: None）
- 图宽：单栏 **8.5 cm**，1.5 栏 **12 cm**，双栏 **17.5 cm**（Graph → Plot Details → 或 File → Page Setup 设 Width/Height）

**② 色盲友好配色（Origin 里手动输入 HEX）**

| 用途 | 颜色 | HEX |
|---|---|---|
| CBS | 蓝 | `#0072B2` |
| ICMP | 橙 | `#E69F00` |
| CGMCC | 绿 | `#009E73` |
| 不可核 / 其他 | 灰 | `#999999` |
| 修正前 | 灰 | `#666666` |
| 修正后 | 蓝 | `#0072B2` |

> **不要用 Origin 默认的彩虹色或红绿对比**——灰度打印会糊，色盲读者分不开。

**③ 导出**
File → Export Graphs → 格式 **TIFF** 或 **PDF**，分辨率 **600 dpi**（位图）或矢量（PDF），
尺寸勾 **Use Page Size**（不要勾 "Auto"）。

---

## 图 1　审计流程与三个分母（漏斗图）

**数据**：`图1_漏斗.tsv`

```
阶段                    数量
全部号级单元格            11,103
可核号（CBS/ICMP/CGMCC）   3,192
号级不一致（修正前）          558
号级不一致（修正后）          402
```

**Origin 做法（Origin 没有现成漏斗图，用水平条形图模拟）**

1. 选中 A、B 两列 → Plot → Bar, Pie, Surface → **Bar**（水平条形）
2. 双击条形 → Plot Details → Spacing 标签：
   - **Gap Between Bars = 0**（条之间不留缝，才像漏斗）
   - **Bar Width = 1**
3. 双击 Y 轴 → **Scale → Reverse**（勾上），让"全部号"在最上面
4. 双击 X 轴 → Scale 设 **0 到 12000**；**隐藏 X 轴的刻度与标签**（Tick Labels → Hide）
5. 加数值标签：双击条形 → Label 标签 → 勾 **Enable** → Position 选 **Center**，字体白/黑视颜色
6. 用 **Rectangle Tool** 画三条水平分隔虚线（或 Insert → Line），把四段隔开
7. 每条右侧加文字说明（Text Tool）：`28.7%`、`17.5%`、`12.6%`

**替代做法（更省事）**：直接做成**四根水平条 + 数值标签**，不追求"漏斗"形状，
在正文里说明"逐级筛选"，读者一样能懂。

---

## 图 2　三属可核率与机构构成（堆叠柱）

**数据**：`图2_机构构成.tsv`

```
属                 CBS   ICMP  CGMCC  其他或不可核   合计
Colletotrichum    1085    286    182      2080       3633
Trichoderma        850     21     99      3230       4200
Diaporthe          542     36     91      2601       3270
```

**做法（推荐 100% 堆叠，突出"构成比例"）**

1. 选中 B～E 四列（CBS / ICMP / CGMCC / 其他）→ Plot → Bar, Pie, Surface → **Stacked Column**
2. 若要百分比：双击图形 → Plot Details → Stack 标签 → **勾 "Normalize to 100%"**
   （或直接用 **Stacked Column 100%** 模板）
3. 配色：按上表 HEX 逐个设；"其他或不可核"用灰 `#999999`
4. 给 CBS/ICMP/CGMCC 三段加数值标签（Label → Position: Center，只加这三段）
5. X 轴标签改成属名（双击 → Tick Labels → 手动输入）

**注意**：ICMP 在 Trichoderma 只有 21 个号（0.5%）——**在 100% 堆叠图里会细到看不见**。
两种处理：
- 若重点是"可核率差异"→ 用**绝对数堆叠**，读者能看到 3,230 的灰色块；
- 若重点是"构成"→ 用 100% 堆叠，但要在图注里写明 ICMP 在另两属几乎为零。

---

## 图 3　三属号级不一致率与 95% 置信区间（点 + 误差线）

**数据**：`图3_率与CI.tsv`（每属 × 三种口径 × 率 + CI 上下限）

**做法**

1. 先做**保守口径**这一组（论文主结果）：
   - 选中"率%"列 → Plot → Basic 2D → **Scatter**
   - X 用属名（或 1/2/3），Y 用率
2. 加误差线：双击数据点 → Plot Details → **Error Bar 标签**
   - Direction: **Y**
   - 勾 **Minus** 和 **Plus**，分别绑定 **CI下限** 与 **CI上限** 列
   - **⚠️ 关键**：Origin 默认是"误差值"（±值），这里给的是**绝对上下限**，
     所以要选 **"Ends of Error Bars"** 而不是 "Plus/Minus"
3. 三种口径并列（推荐）：
   - 把"修正前 / 保留 / 乐观"设成三个 Group（Worksheet → Set As → Group），
     或直接用 **Grouped Scatter**
   - 给三组设不同**形状**（实心圆 / 空心方 / 三角形），不要只靠颜色区分
4. 参考线：给合计的 12.6% 画一条水平虚线（Insert → Line，或双击 Y 轴 → Reference Lines）
5. Y 轴范围建议 **8–22%**，刻度 2%

**图注里必须写**：误差线 = Wilson 95% 置信区间；分母为各属可核号数（1,553 / 970 / 669）。

---

## 数据溯源（供图注引用）

| 图 | 数据来源 |
|---|---|
| 图 1 | 表 1、表 5（2.2、3.1） |
| 图 2 | 各属 `results/<属>/<属>_三列表.tsv` 的机构前缀统计 |
| 图 3 | 表 5 + Wilson 区间（本目录 `图3_率与CI.tsv` 已算好） |

---

## 一个提醒

**图 2 与图 3 的"合计"不要画进图里**——它们是加权结果，与三属并列会误导
（读者会以为是第四个属）。合计只放在表 5 与正文里。
