# HygroThermal

[![tests](https://github.com/HisashiMiura/HygroThermal/actions/workflows/tests.yml/badge.svg)](https://github.com/HisashiMiura/HygroThermal/actions/workflows/tests.yml)
[![docs](https://img.shields.io/readthedocs/hygrothermal/latest)](https://hygrothermal.readthedocs.io/ja/latest/)

建築外皮の熱・水分同時移動解析のプログラム。齋藤宏昭氏の Fortran 版「Whamcavity」（2021 年 8 月）を、Python で作り直したもの。

多層の壁を 1 次元で扱い、各質点の温度と水分化学ポテンシャルを時間を追って求める。通気層の換気、外表面の雨水、層の境目にたまる水（水膜）、材料の中の自由水、木材腐朽を扱う。

## 必要なもの

- Python 3.11 以上
- numpy、pandas、scipy、openpyxl、h5py（`pyproject.toml` の依存）

### 入れ方

パッケージの版は、ロックファイル `requirements-lock.txt` で固定している（いつ、どの PC で入れても同じ版がそろうように）。conda を使う場合は、このプログラム専用の環境を作り、リポジトリの一番上のフォルダで次のように入れる。

```
conda create -n hygrothermal python=3.13
conda activate hygrothermal
python -m pip install -r requirements-lock.txt
python -m pip install -e . --no-deps
```

3 行目でロックファイルの版のパッケージ（テスト用の pytest・matplotlib と、文書用の Sphinx を含む）を入れ、4 行目でこのプログラム自体を入れる。コマンド `hygrothermal` も使えるようになる。`-e` を付けるので、`modules/` のファイルを直すと、入れ直さなくても反映される。

版を固定しなくてよいときは、`python -m pip install -e ".[dev]"` でも入れられる（その時点の最新の版が入る）。

### ロックファイルの版を上げる

ロックファイルは、pip-tools で `pyproject.toml` から作る。版を上げるときは、次のようにして作り直し、入れ直す。

```
python -m pip install pip-tools
pip-compile --upgrade --extra dev --extra docs --strip-extras -o requirements-lock.txt pyproject.toml
python -m pip install -r requirements-lock.txt
```

上げたあとは、テストを流し、計算例の壁の結果が前の版と比べて変わっていないか（変わったなら、どれだけか）を確かめてから、ロックファイルをコミットする。

## 使い方

入力のファイル（JSON）を 1 つ用意して、コマンドに渡す。結果（CSV と HDF5）は、入力のファイルと同じフォルダの `入力のファイルの名前_output/` に書き出す。

```
hygrothermal examples/main_wall.json
```

計算例の `examples/main_wall.json` は、7 層の壁（外側から、サイディング、通気層、せっこうボード、構造用合板、グラスウール、構造用合板、せっこうボード）を、付属の気象データ（`modules/weather.csv`）と室内の余弦の条件で 2 年間計算する入力である。計算に 3 分ほどかかる。

入力のファイル・気象データのファイル・結果のファイルの形は、計算方法の文書の「プログラムの仕様」に示す。

https://hygrothermal.readthedocs.io/ja/latest/contents/02_00_spec.html

同じ壁を計算して結果の図も描くスクリプトもある。図（png）と数値は `examples/output/` に書き出す。

```
python examples/one_year_main_wall.py            # 時間刻み 3600 秒
python examples/one_year_main_wall.py --dt 900   # 時間刻み 900 秒
```

## テスト

```
python -m pytest
```

## 計算方法の文書

計算方法（支配方程式、評価方法、根拠）は、次のサイトで読める。

https://hygrothermal.readthedocs.io/

文書は `docs/source/` に Markdown で書き、Sphinx でサイトにする。手元で作るときは、次のようにする。

```
python -m sphinx -b html docs/source docs/build/html
```

`docs/build/html/index.html` を開くと読める。上のサイトは、GitHub に push すると Read the Docs が作り直す（設定は `.readthedocs.yaml`）。

## フォルダの構成

| フォルダ・ファイル | 中身 |
|---|---|
| `modules/` | プログラムの本体と、材料の物性値・気象データのファイル |
| `test/` | テスト |
| `examples/` | 計算例のスクリプト |
| `docs/source/` | 計算方法の文書（Sphinx） |
| `docs/basis_checks/` | 計算方法の文書の根拠のページに載せている試算を出すスクリプト |
| `docs/design/` | プログラムの組み立ての設計と、決めたことの記録 |
| `docs/reference/` | 材料の物性値の確かめなど、調べたことの記録 |
| `docs/note_*.md` | Fortran 版との違いと、作成者に確かめたいこと |
| `_fortran/` | もとの Fortran 版と、その入力のファイル（参考として残している） |
| `notes/` | 調べものに使ったノートブック |
| `requirements-lock.txt` | パッケージの版を固定したロックファイル（pip-tools で `pyproject.toml` から作る） |
| `.github/workflows/` | GitHub Actions のテスト（`tests.yml`：push のたびにロックファイルの版で、`latest.yml`：週に 1 回最新の版で） |
