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
# colon_fence・deflist は、記号及び添え字の一覧（:::{container} symbols と、記号 / : 説明 の書き方）に使う
myst_enable_extensions = ['dollarmath', 'amsmath', 'colon_fence', 'deflist']
source_suffix = {'.md': 'markdown'}

# 見た目（hc-energy と同じ Read the Docs のテーマ）
html_theme = 'sphinx_rtd_theme'
html_title = 'HydroThermal 計算方法'

# 見た目の調整（式の番号の位置、長い式、記号の一覧）
html_static_path = ['_static']
html_css_files = ['custom.css']

# 表・図に番号を付ける
numfig = True

# 式の番号はページごとに 1 から振る（True だとサイト全体の通し番号になる）
math_numfig = False
