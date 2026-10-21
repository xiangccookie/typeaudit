# -*- coding: utf-8 -*-
"""NCBI 结果表：把「联网拿到的结果」按属存成一张表，重跑零下载。

为什么要有它（本机实测的数据）：
  · .cache/*.gb 能防重下，但是 747 个碎片文件、最坏 8.48 GiB、且**无法分享**；
  · 表只存【判据真正消费的结果】：从 19 个 .gb（6,182,429 B）抽出的 575 条值，
    写成表只有 31,020 B —— 比原文小 99.5%。

设计要点（**可复现性的关键**，改动前请先读这三条）：
  一、**方法级记忆，不是字段级**。键 = (方法名, 规范化后的完整实参)。
      存的是方法的**返回值本身** ⇒ 重放按构造就与原结果逐字节相同，
      不存在「字段白名单漏一项 → 某个判据静默失效」的风险。
      （那个坑本项目踩过：ICMP 源没加载时全表被标「无目录」，不报错、只是不判。）
  二、返回值一律以 JSON 存表；tuple 用 {"__tuple__": [...]} 包裹，
      读回时还原成 tuple（TaxonomyAPI.type_material_of 返回的就是 tuple，
      直接存 JSON 会变成 list，调用方按元组解包就会炸）。
  三、未命中 → 照常联网，并把结果**补写进表** ⇒ 表随使用自动长大。
      meta 里记录代码指纹，指纹变了只是**提示可能过期**，不自动作废（由人决定是否重抓）。

表结构（TSV）：
    方法    参数    值(JSON)
值之所以是 JSON：不同方法的返回形状不同（list / tuple / dict / int / None），
JSON 是保证**精确重放**的最简载体。要看人读的形态请用 tools/ncbi_table_info.py。
"""
import functools
import hashlib
import inspect
import io
import json
import os
import time

FORMAT = 1
_TUPLE = "__tuple__"
_MISS = object()          # 「表里没有」的哨兵 —— 必须与「存了 None」区分开


# ---------------------------------------------------------------- JSON 往返
def _to_json(v):
    """转成可 JSON 化的形状。不可序列化时返回 None（调用方据此放弃记忆）。"""
    if v is None or isinstance(v, (str, int, float, bool)):
        return v
    if isinstance(v, tuple):
        return {_TUPLE: [_to_json(x) for x in v]}
    if isinstance(v, list):
        return [_to_json(x) for x in v]
    if isinstance(v, dict):
        return {str(k): _to_json(x) for k, x in v.items()}
    return None


def _from_json(v):
    if isinstance(v, dict):
        if _TUPLE in v and len(v) == 1:
            return tuple(_from_json(x) for x in v[_TUPLE])
        return {k: _from_json(x) for k, x in v.items()}
    if isinstance(v, list):
        return [_from_json(x) for x in v]
    return v


# ---------------------------------------------------------------- 键
def _norm(x):
    if x is None or isinstance(x, (str, int, float, bool)):
        return x
    return repr(x)


def make_key(method, arguments, owner=""):
    """键 = 【源类名.方法名】+ 完整实参（含默认值）。

    为什么要带类名：`type_material_of` 在 GenBankAPI 和 TaxonomyAPI 上**同名**，
    只按方法名 + 参数组键的话，两个源会共用一条记录 ——
    A 源会拿到 B 源的结果，而且**不报错**。
    （本次实测侥幸没撞上：一个 retmax=40、一个 retmax=3；靠巧合不是靠设计。）

    为什么要含默认值：让 `f(x)` 与 `f(x, retmax=300)` 落到同一个键，
    否则同一次调用换个写法就会重下一次。
    """
    name = "{}.{}".format(owner, method) if owner else method
    return "{}|{}".format(name, json.dumps(
        {k: _norm(v) for k, v in arguments.items()},
        sort_keys=True, ensure_ascii=False))


# ---------------------------------------------------------------- 表
class MemoTable(object):
    def __init__(self, path, genus="", record=True):
        self.path = path
        self.genus = genus
        self.record = record
        self.data = {}
        self.hits = 0
        self.misses = 0
        self.skipped = 0
        self.load_errors = []
        self._pending = 0
        self._load()

    # ---- 读写
    def _load(self):
        if not self.path or not os.path.exists(self.path):
            return
        try:
            with io.open(self.path, encoding="utf-8") as fh:
                head = fh.readline().rstrip("\n").split("\t")
                for i, ln in enumerate(fh, 2):
                    parts = ln.rstrip("\n").split("\t")
                    if len(parts) < 3:
                        continue
                    try:
                        # 存 JSON 形态即可（与 put 保持同一约定；get 负责转回）
                        self.data["{}\t{}".format(parts[0], parts[1])] = json.loads(parts[2])
                    except Exception as e:
                        self.load_errors.append("第{}行：{}".format(i, e))
        except Exception as e:
            self.load_errors.append("整体：{}: {}".format(type(e).__name__, e))
        # 载入时把「键」还原成 make_key 的形态（方法|参数JSON）
        fixed = {}
        for k, v in self.data.items():
            m, _, p = k.partition("\t")
            fixed["{}|{}".format(m, p)] = v
        self.data = fixed

    def get(self, key):
        # ⚠️ 内存里一律存【JSON 形态】，取的时候才转回 Python 形态。
        #    原来这里直接返回 data[key]，于是 tuple 会以 dict 的形态交回去
        #    （TaxonomyAPI.type_material_of 返回 (taxid, list)，调用方按元组解包就炸）。
        #    这是单测抓出来的：二次调用返回类型 dict ≠ 首次的 tuple。
        if key in self.data:
            return _from_json(self.data[key])
        return _MISS

    def put(self, key, value):
        if not self.record:
            return
        j = _to_json(value)
        if j is None and value is not None:
            # 形状不可 JSON 化 —— 宁可不记，也绝不写入 None
            # （写 None 会被当成「合法结果」在下次命中，属于最坏的一类静默错误）
            self.skipped += 1
            return
        if key not in self.data:
            self._pending += 1
        self.data[key] = j
        if self._pending >= 25:
            self.save()

    def save(self, force=False):
        if not self.path or (not self._pending and not force):
            return False
        d = os.path.dirname(os.path.abspath(self.path))
        if d:
            try:
                os.makedirs(d, exist_ok=True)
            except Exception:
                pass
        tmp = self.path + ".tmp"
        try:
            with io.open(tmp, "w", encoding="utf-8", newline="\n") as fh:
                fh.write("方法\t参数\t值(JSON)\n")
                for k in sorted(self.data):
                    m, _, p = k.partition("|")
                    fh.write("{}\t{}\t{}\n".format(
                        m, p, json.dumps(self.data[k], ensure_ascii=False,
                                         sort_keys=True)))
            os.replace(tmp, self.path)
            self._pending = 0
        except Exception:
            return False
        self._write_meta()
        return True

    def _write_meta(self):
        try:
            import collections
            by = collections.Counter(k.split("|", 1)[0] for k in self.data)
            meta = {
                "format": FORMAT,
                "genus": self.genus,
                "updated": time.strftime("%Y-%m-%d %H:%M:%S"),
                "entries": len(self.data),
                "by_method": dict(by),
                "hits": self.hits,
                "misses": self.misses,
                "skipped_unserializable": self.skipped,
                "code_fingerprint": code_fingerprint(),
                "_说明": ("本表按【方法级记忆】存放：键 = 方法名|完整实参，值 = 该方法的返回值(JSON)。"
                          "重放按构造与原结果逐字节相同。"
                          "code_fingerprint 与当前代码不一致 = 表可能已过期，"
                          "需要重抓时请显式删除本表（不会自动作废）。"),
            }
            mp = os.path.splitext(self.path)[0] + ".meta.json"
            with io.open(mp, "w", encoding="utf-8", newline="\n") as fh:
                fh.write(json.dumps(meta, ensure_ascii=False, indent=2))
        except Exception:
            pass

    def stats(self):
        return {"entries": len(self.data), "hits": self.hits,
                "misses": self.misses, "skipped": self.skipped}


# ---------------------------------------------------------------- 代码指纹
def code_fingerprint():
    h = hashlib.sha256()
    base = os.path.dirname(os.path.abspath(__file__))
    for fn in ("sources.py", "ncbitable.py"):
        p = os.path.join(base, fn)
        try:
            with io.open(p, "rb") as fh:
                h.update(fh.read())
        except Exception:
            pass
    return h.hexdigest()[:16]


# ---------------------------------------------------------------- 挂载
_SKIP = {"save", "stats", "n_rows"}


def attach(obj, table, skip=()):
    """给源对象的**所有公开方法**套上记忆层。

    之所以不逐个方法手写：将来新增判据方法会自动被覆盖，
    不会出现「加了新方法但忘了记忆 → 又重下」。
    """
    if table is None or not getattr(table, "record", True) and table.data is None:
        return 0
    n = 0
    for name in dir(type(obj)):
        if name.startswith("_") or name in _SKIP or name in skip:
            continue
        if not callable(getattr(type(obj), name, None)):
            continue
        # ⚠️ 必须取【绑定方法】`getattr(obj, name)`，**不能**取 `getattr(type(obj), name)`：
        #    后者是未绑定函数、没有 self。它一旦被 setattr 成实例属性，Python 就
        #    不再自动绑定 —— 调用时第一个实参会被塞进 self，
        #    表现为 `missing 1 required positional argument: 'organism'`。
        #    更坏的是这个 TypeError 会被 collect() 的 except 吞掉：
        #    结果退化却【不报错、也不联网】，看起来像"表生效了"。
        fn = getattr(obj, name)
        try:
            sig = inspect.signature(fn)     # 绑定方法的签名已排除 self
        except (TypeError, ValueError):
            continue
        setattr(obj, name, _wrap(type(obj).__name__, name, fn, sig, table))
        n += 1
    return n


def _wrap(owner, name, fn, sig, table):
    @functools.wraps(fn)
    def wrapper(*a, **kw):
        try:
            b = sig.bind(*a, **kw)
            b.apply_defaults()
            key = make_key(name, b.arguments, owner=owner)
        except Exception:
            # 绑不上参数（*args 之类）→ 老老实实走原路，不冒险
            return fn(*a, **kw)
        got = table.get(key)
        if got is not _MISS:
            table.hits += 1
            return got
        v = fn(*a, **kw)
        table.misses += 1
        table.put(key, v)
        return v
    return wrapper
