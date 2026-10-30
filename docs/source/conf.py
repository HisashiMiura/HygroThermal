"""Sphinx の設定（計算方法のサイト）。

手元でサイトを作るとき（リポジトリの一番上のフォルダで）：
    python -m sphinx -b html docs/source docs/build/html
できたサイトは docs/build/html/index.html。Read the Docs では .readthedocs.yaml の設定で同じように作る。
"""

project = 'HydroThermal'
author = 'Hisashi Miura'
copyright = '2026, Hisashi Miura'
language = 'ja'

# Markdown（MyST）で書く。数式は $…$・$$…$$ で書き、$$…$$ (名前) で式に番号を付けて参照できるようにする
extensions = ['myst_parser']
myst_enable_extensions = ['dollarmath', 'amsmath']
source_suffix = {'.md': 'markdown'}

# 見た目（hc-energy と同じ Read the Docs のテーマ）
html_theme = 'sphinx_rtd_theme'
html_title = 'HydroThermal 計算方法'

# 式・表に番号を付ける（式の番号はページごとに 1 から）
numfig = True
