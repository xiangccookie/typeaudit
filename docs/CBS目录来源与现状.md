# CBS 目录是怎么来的（证据链 + 现状）

> 记于 2026-10-16。**目的**：避免"这文件哪来的"再被忘掉，也避免误以为"抓不了就没救"。

## 一、当时的获取方式（从 10-11 会话日志还原）

日志：`~/.dsh/dsh-bio-genie/log/2026-10-11.jsonl`，其中三条 bio_python 记录：

| 行 | 耗时 | 代码首行（preview 仅存 200 字符） |
|---|---|---|
| 12 | 1.2 s | `# 目的：探 westerdijk / wi.knaw.nl 的真菌目录页面，看后端接口在哪。` |
| 17 | 9.9 s | `# 目的：确认 BioAware 后端可用，并找出「菌株检索」的接口。`<br>`# 基地址来自 https://wi.knaw.nl/config.json：https://webservices.bio-aware.com/cbsdatabase`<br>`# 必需请求头：WebsiteId: 88` |
| 21 | 13.7 s | `# 目的：用正确的参数名 p_PageName 取 fungal_table 页面配置…` |

**⇒ 结论：不是抓 HTML，而是走 BioAware（BioloMICS）的 REST 后端。**

链路：

```
wi.knaw.nl（Angular SPA，HTML 抓不到）
   └─ config.json  →  {"baseUrl":"https://webservices.bio-aware.com/cbsdatabase",
                        "websiteId":88}
        └─ 带 WebsiteId:88 调 BioAware 后端
             └─ cbs_catalogue_raw.jsonl   62,772 行
                  └─ cbs_catalogue.tsv    （norm_key = 去点去空格大写）
```

**原始数据字段**（`cbs_catalogue_raw.jsonl` 首行）：
```json
{"name": "CBS 14", "id": 1, "status": "Type of Rhodotorula gracilis Rennerfelt", "taxon": null, "organism": "Yeast"}
```

## 二、现状（2026-10-16 13:26 实测）

| 环节 | 状态 |
|---|---|
| `wi.knaw.nl/config.json` | ✅ 活着，仍返回 baseUrl + websiteId:88 |
| BioAware 后端 `/` | ✅ 活着（Swagger UI 页面） |
| Swagger 文档 `/swagger/v1|v2/swagger.json` | ❌ **HTTP 500**（多种请求头组合都不行） |
| 常见端点探测（19 个 GET） | ❌ 全部 404 |
| 当时的完整抓取代码 | ❌ **未留存**（日志与记忆只存 200 字符预览 + code_hash） |

**⇒ 自动重抓目前做不到**：后端在，但端点未知、文档不可得、代码无存。

## 三、怎么办

| 需求 | 做法 |
|---|---|
| **继续用**（主流） | 直接用现有 `cbs_catalogue.tsv`（62,761 行） |
| **补个别新号** | `python tools/fetch_refdata.py cbs --add CBS <号> "<状态>"` |
| **体检覆盖** | `python tools/fetch_refdata.py cbs --check` |
| **想彻底恢复自动抓取** | 需要找回当时的完整代码（可能在其它机器/会话），或从 `swagger/v1/swagger.json` 恢复（当前 500）；实在要自己逆，就从 `index.js` 里那两个文档地址与 `config.json` 入手 |

> ⚠️ **不要**为了重抓去跑 6.3 万次请求——现有副本已覆盖 62,761 行，重抓换来的是同一份数据。
