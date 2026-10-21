# -*- coding: utf-8 -*-
"""控制台编码兜底：把 stdout/stderr 强制成 UTF-8。

为什么需要它：
  Windows 中文版的控制台默认是 GBK(cp936)。本项目大量输出 `↔ ⇒ ⚠️ 🔴 🟠 🔵 ✅ ❌`
  这类符号，在 GBK 下 `print()` 会抛
      UnicodeEncodeError: 'gbk' codec can't encode character '\u2194'
  直接让程序中途崩掉 —— 在开发者本机（设了 UTF-8 环境变量）看不出来，
  换一台干净机器就必崩。这是"别人能不能跑"的头号坑。

做法：
  用 Python 3.7+ 的 `stream.reconfigure()`，把编码改成 UTF-8、把无法编码的字符
  替换掉而不是抛异常。改不动时静默跳过（比如输出被重定向到不支持的对象）。

用法（三种情形都已覆盖）：
  1. 任何 `import typeaudit` 的代码 —— 由 `typeaudit/__init__.py` 自动调用；
  2. `tools/` 里不导入 typeaudit 的独立脚本 —— 脚本顶部已插入同样的两行；
  3. 想手动兜底：`from typeaudit._console import setup; setup()`
"""

import sys

_DONE = False


def setup():
    """把标准输出/错误改成 UTF-8（幂等，可重复调用）。返回是否改动成功。"""
    global _DONE
    if _DONE:
        return True
    ok = False
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.reconfigure(encoding="utf-8", errors="replace")
            ok = True
        except Exception:
            # 输出被重定向成不支持 reconfigure 的对象时，保持原样即可
            pass
    _DONE = True
    return ok


# 导入即生效：这样 `from typeaudit._console import setup` 与
# `import typeaudit._console` 两种写法都能兜住。
setup()
