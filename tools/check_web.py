#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
check_web.py —— 网页版内联 JS 的语法自检 + 元素 id 对账。

【为什么必须有它】
    pipeline_web.py 的页面是一整块 Python 字符串里塞的 HTML+JS。
    JS 只要有一个语法错误，**整段脚本都不执行 ⇒ 页面上所有按钮全部点不动**，
    而服务本身正常启动、页面也正常返回 200 —— 从服务端完全看不出来。
    本机实测踩过：`let TAB` 被声明了两次（原脚本一次 + 新加的全局变量一次），
    于是 SyntaxError，全部按钮失效。

    ⚠️ 判别特征：**页面能打开、能滚动、但点任何按钮都没反应** ⇒ 先怀疑 JS 语法错误，
       不要去看服务端日志（那里什么都不会有）。

用法
    python tools/check_web.py            # 语法检查 + id 对账
    python tools/check_web.py --quiet    # 只报结论
"""
import sys
try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

import argparse, io, os, re, shutil, subprocess, tempfile
from collections import Counter

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def extract(path):
    x = io.open(path, encoding="utf-8").read()
    i = x.index('PAGE = r"""') + len('PAGE = r"""')
    j = x.index('"""', i)
    page = x[i:j]
    m = re.search(r"<script>(.*)</script>", page, re.S)
    return page, (m.group(1) if m else "")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--file", default=os.path.join(REPO, "tools", "pipeline_web.py"))
    ap.add_argument("--quiet", action="store_true")
    a = ap.parse_args()

    page, js = extract(a.file)
    ok = True
    print("=" * 88)
    print("网页内联 JS 自检")
    print("=" * 88)
    print("  PAGE {:,} 字符；内联 JS {:,} 字符".format(len(page), len(js)))

    # ① 语法检查（有 node 就用 node，没有就退回括号配平）
    node = shutil.which("node")
    if node:
        fd, tmp = tempfile.mkstemp(suffix=".js")
        os.close(fd)
        io.open(tmp, "w", encoding="utf-8", newline="\n").write(js)
        r = subprocess.run([node, "--check", tmp], capture_output=True, text=True,
                           encoding="utf-8", errors="replace")
        os.remove(tmp)
        if r.returncode == 0:
            print("  ① JS 语法：✅ 通过（node --check）")
        else:
            ok = False
            print("  ① JS 语法：❌ 失败 —— 这会让【所有按钮点不动】")
            print("     " + r.stderr.strip().split("\n")[-4:][0][:160])
    else:
        bal = {c: js.count(c) for c in "{}()[]"}
        bad = [c for c in "{}()" if bal[c] % 2]
        print("  ① 没装 node，退回括号配平：{}".format(
            "✅ 通过" if not bad else "❌ 不配平 {}".format(bad)))
        if bad:
            ok = False

    # ② 重复声明（let/const 重复就是 SyntaxError，正是踩过的坑）
    # ⚠️ 只数【顶层】声明（行首无缩进）。函数内部的同名 const 是合法的，
    #    把它们也算进来会报出一堆假重复（本机实测：r/j/body 各十几处，全是函数内局部的）。
    decl = re.findall(r"^(?:let|const)\s+([A-Za-z_$][\w$]*)\s*=", js, re.M)
    dup = {k: v for k, v in Counter(decl).items() if v > 1}
    if dup:
        ok = False
        print("  ② 重复声明：❌ {}".format(dup))
        print("     let/const 同名重复声明是 SyntaxError，整段脚本都会失效。")
    else:
        print("  ② 重复声明：✅ 无")

    # ③ JS 引用的元素 id 是否都在 HTML 里
    # ⚠️ 必须先剥掉注释再扫 —— 注释里写个 $('xxx') 会报假阳性
    #    （本机实测：我在注释里解释「不要用 $('biocol')」，检查器就把 biocol 报成缺失元素）。
    #    剥行注释时要避开 `://`（URL 里的双斜杠）。
    code = re.sub(r"/\*[\s\S]*?\*/", " ", js)
    code = re.sub(r"(^|[^:])//[^\n]*", r"\1", code, flags=re.M)
    ids = set(re.findall(r'\$\("([^"]+)"\)', code)) | set(re.findall(r"\$\('([^']+)'\)", code))
    html_ids = set(re.findall(r'id="([^"]+)"', page))
    missing = sorted(ids - html_ids)
    if missing:
        ok = False
        print("  ③ 元素 id 对账：❌ 缺 {}".format(missing))
        print("     \\(id) 取不到元素会在运行时抛错，中断该事件处理。")
    else:
        print("  ③ 元素 id 对账：✅ 全部存在（{} 个）".format(len(ids)))

    print()
    print("  结论：{}".format("✅ 通过" if ok else "❌ 有问题 —— 修完再启动页面"))
    return 0 if ok else 2


if __name__ == "__main__":
    sys.exit(main())
