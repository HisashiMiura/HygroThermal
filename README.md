# HydroThermal

建築外皮の熱・水分同時移動解析のプログラム。齋藤宏昭氏の Fortran 版「Whamcavity」（2021 年 8 月）を、Python で作り直したもの。

多層の壁を 1 次元で扱い、各質点の温度と水分化学ポテンシャルを時間を追って求める。通気層の換気、外表面の雨水、層の境目にたまる水（水膜）、材料の中の自由水、木材腐朽を扱う。

## 必要なもの

- Python 3.11 以上
- numpy、pandas、scipy、openpyxl（`pyproject.toml` の依存）

リポジトリの一番上のフォルダで、次のように入れる。

```
pip install -e ".[dev]"
```

`[dev]` は、テスト（pytest）と図（matplotlib）のためのもの。計算方法の文書を作るときは `[docs]` も入れる。

## テスト

```
python -m pytest
```

## 例の計算

7 層の壁（外側から、サイディング、通気層、せっこうボード、構造用合板、グラスウール、構造用合板、せっこうボード）を、例の気象データ（`modules/weather.csv`）と室内の余弦の条件で 2 年間計算し、結果の図を描く。

```
python examples/one_year_main_wall.py            # 時間刻み 3600 秒
python examples/one_year_main_wall.py --dt 900   # 時間刻み 900 秒
```

結果の数値（npz）と図（png）は `examples/output/` に書き出す（git には入れない）。

## 計算方法の文書

計算方法（支配方程式、評価方法、根拠）は `docs/source/` に Markdown で書き、Sphinx でサイトにする。

```
python -m sphinx -b html docs/source docs/build/html
```

`docs/build/html/index.html` を開くと読める。Read the Docs でも同じものを作る（`.readthedocs.yaml`）。

## フォルダの構成

| フォルダ | 中身 |
|---|---|
| `modules/` | プログラムの本体と、材料の物性値・気象データのファイル |
| `test/` | テスト |
| `examples/` | 例の計算のスクリプト |
| `docs/source/` | 計算方法の文書（Sphinx） |
| `docs/basis_checks/` | 計算方法の文書の根拠のページに載せている試算を出すスクリプト |
| `docs/design/` | プログラムの組み立ての設計と、決めたことの記録 |
| `docs/reference/` | 材料の物性値の確かめなど、調べたことの記録 |
| `docs/note_*.md` | Fortran 版との違いと、作成者に確かめたいこと |
| `_fortran/` | もとの Fortran 版と、その入力のファイル（参考として残している） |
| `notes/` | 調べものに使ったノートブック |
