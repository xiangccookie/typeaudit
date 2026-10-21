# -*- coding: utf-8 -*-
# --- console utf8 bootstrap（Windows 控制台默认 GBK，输出 ↔ ⚠️ 🔴 会崩；详见 src/typeaudit/_console.py）
import sys as _sys
for _s in (_sys.stdout, _sys.stderr):
    try:
        _s.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass
r"""对账：结果表里的种，有没有从所有选项卡里消失的。

    python tools\check_conservation.py --out out_audit --genus Colletotrichum

⚠️ 断言是【覆盖】不是【划分】：一个种可以同时出现在多个选项卡里
   （一个号干净 + 一个号有问题），所以**不能**要求各选项卡种数之和等于种数。
   退出码：0 = 无种遗漏；1 = 有遗漏；2 = 找不到结果表。
"""
import argparse
import os
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(REPO, "src"))
sys.path.insert(0, HERE)

import pipeline_web as PW      # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="out_audit")
    ap.add_argument("--genus", default="Colletotrichum")
    a = ap.parse_args()
    od = a.out if os.path.isabs(a.out) else os.path.join(REPO, a.out)
    rep = PW.conservation_check(od, a.genus)
    if not rep:
        print("找不到结果表：{}".format(os.path.join(od, "{}_三列表.tsv".format(a.genus.lower()))))
        return 2
    print(PW.conservation_line(rep))
    print()
    print("说明：五种去向对【种】不是划分 —— 一个种可以同时出现在多个选项卡里，")
    print("      所以断言是「每个种至少出现在一个去向里」（覆盖），不是各表计数相加。")
    return 0 if rep["ok"] else 1


if __name__ == "__main__":
    sys.exit(main())
