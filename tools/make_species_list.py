#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""make_species_list.py —— 取一个属的有效双名清单。

为什么要单独一个脚本：`TaxonomyAPI.species_in_genus()` 要为每个种级 taxon 取一次
分类学记录（Trichoderma 是 2,574 个），耗时可达数十分钟，不能塞在一次工具调用里。
所以做成可后台跑、可续的脚本 —— 取过的记录已缓存在 .cache，重跑很快。

过滤规则与炭疽菌清单一致：只保留规范双名（两个词、种加词全小写），
剔除 sp./cf./aff. 这类非种级写法。
"""
import argparse, io, os, sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "src"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genus", required=True)
    ap.add_argument("--out", default="")
    a = ap.parse_args()
    os.environ.setdefault("NCBI_EMAIL", "audit@example.org")

    from typeaudit import sources as S
    srcs = S.load(os.path.join(REPO, "config", "sources.json"))
    tx = srcs.get("taxonomy")
    if tx is None:
        print("❌ 没有 taxonomy 源"); return 2

    print("[{}] 开始取种清单（预计耗时较长，逐条取分类学记录）".format(a.genus), flush=True)
    names = tx.species_in_genus(a.genus, retmax=5000, limit=None)
    print("[{}] 取到有效双名 {} 个".format(a.genus, len(names)), flush=True)

    out = a.out or os.path.join(REPO, "data", "species_{}_all.txt".format(a.genus.lower()))
    io.open(out, "w", encoding="utf-8", newline="\n").write("\n".join(names) + "\n")
    print("已写 → {}".format(out), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
