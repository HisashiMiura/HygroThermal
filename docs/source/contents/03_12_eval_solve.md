# 1 ステップの解き方

## 1 はじめに

### 1.1 概要

本節では、「{doc}`熱・水分の収支 <03_11_eval_balance>`」の収支の式を満たすステップ $n + 1$ の温度と水分化学ポテンシャルを、ニュートン・ラフソン法（以下、ニュートン法）で求める方法を定める。熱と水分の収支の式は、温度と水分化学ポテンシャルの両方を含むので、すべての質点の熱と水分を一緒に解く。

### 1.2 入力

:::{container} symbols
$T_{i,n}$
: ：ステップ $n$ の質点 $i$ の温度（K）

$\mu_{i,n}$
: ：ステップ $n$ の質点 $i$ の水分化学ポテンシャル（J/kg）

$\varphi_{i,n}$
: ：ステップ $n$ の質点 $i$ の相対湿度。 $T_{i,n}$ と $\mu_{i,n}$ から、「{doc}`湿り空気と水の性質 <03_02_eval_air_water>`」の式 (3) で求める（%）
:::

また、「{doc}`熱・水分の収支 <03_11_eval_balance>`」の式 (8) と式 (9) を使う。

### 1.3 出力

:::{container} symbols
$T_{i,n+1}$
: ：ステップ $n + 1$ の質点 $i$ の温度（K）

$\mu_{i,n+1}$
: ：ステップ $n + 1$ の質点 $i$ の水分化学ポテンシャル（J/kg）
:::

## 2 記号及び添え字

### 2.1 記号

:::{container} symbols
$\boldsymbol{e}_q$
: ： $q$ 番目の成分だけが 1 で、ほかの成分が 0 のベクトル（-）

$F_p$
: ：未知数の並び $\boldsymbol{x}$ の、式 $p$ の残差（-）

$h_q$
: ：ヤコビ行列を求めるときに、未知数 $q$ をずらす量（K、J/kg）

$J_{pq}$
: ：ヤコビ行列の $p$ 行 $q$ 列（-/K、-/(J/kg)）

$L_{\mathrm{v}}$
: ：水の蒸発潜熱（J/kg）

$N$
: ：質点の数（-）

$s$
: ：修正の幅を半分にした回数（-）

$\boldsymbol{x}$
: ：未知数の並び（ベクトル）。成分は $x_p$（K、J/kg）

$x_p$
: ：未知数 $p$（K、J/kg）

$x_{\max,p}$
: ：未知数 $p$ の上限（K、J/kg）

$\Delta x_p$
: ：未知数 $p$ の修正の量（K、J/kg）

$\varepsilon$
: ：計算機イプシロン。 $\varepsilon = 2^{-52}$（-）

$\varepsilon_{\mathrm{h}}$
: ：熱の残差の許容値。 $\varepsilon_{\mathrm{h}} = 10^{-6}$（W/m²）

$\varepsilon_{\mathrm{m}}$
: ：水分の残差の許容値（kg/(m² s)）

$\varphi_{i,n}$
: ：ステップ $n$ の質点 $i$ の相対湿度（%）
:::

### 2.2 添え字

:::{container} symbols
$i$
: ：質点の番号

$n$
: ：時刻のステップ

$p$
: ：残差の式の番号。 $0, 1, \ldots, 2N - 1$

$q$
: ：未知数の番号。 $0, 1, \ldots, 2N - 1$

$(m)$
: ：ニュートン法の修正の回数。右上に括弧を付けて書く
:::

## 3 未知数と残差

未知数は、すべての質点の、ステップ $n + 1$ の温度と水分化学ポテンシャルとし、式 {eq}`sv-x` のように、質点ごとに温度と水分化学ポテンシャルを交互に並べる。

$$
x_{2i} = T_{i,n+1}, \quad x_{2i+1} = \mu_{i,n+1}
$$ (sv-x)

残差は、「{doc}`熱・水分の収支 <03_11_eval_balance>`」の熱の収支の式 (9) と水分の収支の式 (8) の、左辺から右辺を引いたものとし、それぞれの許容値で割って単位をそろえ、式 {eq}`sv-f` のように未知数と同じ順に並べる。熱の残差と水分の残差の大きさは、許容値の何倍かで比べる。

$$
F_{2i} = \frac{\left( \text{式 (9) の左辺} \right) - \left( \text{式 (9) の右辺} \right)}{\varepsilon_{\mathrm{h}}}, \quad
F_{2i+1} = \frac{\left( \text{式 (8) の左辺} \right) - \left( \text{式 (8) の右辺} \right)}{\varepsilon_{\mathrm{m}}}
$$ (sv-f)

水分の残差の許容値 $\varepsilon_{\mathrm{m}}$ は、潜熱に換算すると熱の残差の許容値と同じ大きさになるように、式 {eq}`sv-eps-m` により計算される。

$$
\varepsilon_{\mathrm{m}} = \frac{\varepsilon_{\mathrm{h}}}{L_{\mathrm{v}}}
$$ (sv-eps-m)

熱の残差の許容値 $\varepsilon_{\mathrm{h}} = 10^{-6}$ W/m² は、温度にすると、およそ $10^{-7}$ K の精度である。

## 4 ニュートン法

### 4.1 最初の仮の値

未知数の最初の仮の値 $x_p^{(0)}$ は、ステップ $n$ の値とする。

$$
x_{2i}^{(0)} = T_{i,n}, \quad x_{2i+1}^{(0)} = \mu_{i,n}
$$ (sv-x0)

### 4.2 未知数の上限

水分化学ポテンシャルの上限は 0（相対湿度 100%）とし、温度には上限を設けない。

$$
x_{\max,2i} = \infty, \quad x_{\max,2i+1} = 0
$$ (sv-xmax)

水分化学ポテンシャルが 0 より大きい範囲では、相対湿度が 100% で頭打ちになり、残差が水分化学ポテンシャルによらなくなる。上限は、ニュートン法の修正がその範囲へ行き過ぎないようにするためのものである。

### 4.3 ヤコビ行列

$m$ 回目の修正のヤコビ行列 $J_{pq}$ は、式 {eq}`sv-j` により、未知数 $q$ を $h_q$ だけずらしたときの残差の変化から求める。

$$
J_{pq} = \frac{\partial F_p}{\partial x_q} \approx \frac{F_p \left( \boldsymbol{x}^{(m)} + h_q \boldsymbol{e}_q \right) - F_p \left( \boldsymbol{x}^{(m)} \right)}{h_q}
$$ (sv-j)

ここで、 $\boldsymbol{e}_q$ は $q$ 番目の成分だけが 1 で、ほかの成分が 0 のベクトルである。ずらす量 $h_q$ は、式 {eq}`sv-h` により計算される。ずらすと上限を超える未知数は、反対の向きにずらす（上限の外の、残差が変わらない範囲で傾きを求めないため）。

$$
h_q =
\begin{cases}
\sqrt{\varepsilon} \max \left( \left| x_q^{(m)} \right|, 1 \right) & \left( x_q^{(m)} + \sqrt{\varepsilon} \max \left( \left| x_q^{(m)} \right|, 1 \right) \le x_{\max,q} \right) \\
- \sqrt{\varepsilon} \max \left( \left| x_q^{(m)} \right|, 1 \right) & (\text{それ以外})
\end{cases}
$$ (sv-h)

質点 $i$ の残差は、質点 $i - 1$・ $i$・ $i + 1$ の温度と水分化学ポテンシャルにしかよらない。そのため、ヤコビ行列は、対角の下に 3 本、上に 3 本の帯の形になる。

### 4.4 修正

修正の量 $\Delta x_p$ は、式 {eq}`sv-dx` の連立 1 次方程式を解いて求める。

$$
\sum_{q} J_{pq} \thinspace \Delta x_q = - F_p \left( \boldsymbol{x}^{(m)} \right)
$$ (sv-dx)

次の仮の値 $x_p^{(m+1)}$ は、式 {eq}`sv-update` により計算される。

$$
x_p^{(m+1)} = \min \left( x_p^{(m)} + \frac{\Delta x_p}{2^s}, \thinspace x_{\max,p} \right)
$$ (sv-update)

$s$ は、 $s = 0, 1, \ldots, 10$ の順に試し、式 {eq}`sv-accept` のどちらかを満たす最初の値とする（修正の幅を半分ずつ小さくする直線探索）。 $\boldsymbol{x}^{(m+1)}$ で残差が計算できない（0 で割るなど）ときは、満たさないものとする。 $s = 10$ でも満たさなければ、収束しなかったものとする。

$$
\sum_{p} \left( F_p \left( \boldsymbol{x}^{(m+1)} \right) \right)^2 < \sum_{p} \left( F_p \left( \boldsymbol{x}^{(m)} \right) \right)^2
\quad \text{または} \quad
\max_{p} \left| F_p \left( \boldsymbol{x}^{(m+1)} \right) \right| \le 1
$$ (sv-accept)

答えの近くでは $s = 0$ で式 {eq}`sv-accept` を満たすので、ふつうのニュートン法と同じになる。

### 4.5 収束の判定

修正を 1 回以上行い、すべての残差が許容値以下（ $\max_p \left| F_p \left( \boldsymbol{x}^{(m)} \right) \right| \le 1$）になったら、収束したものとし、そのときの $\boldsymbol{x}^{(m)}$ を答えとする。

$$
T_{i,n+1} = x_{2i}^{(m)}, \quad \mu_{i,n+1} = x_{2i+1}^{(m)}
$$ (sv-ans)

最初の仮の値の残差が許容値以下でも、1 回は修正する。平衡の近くで、1 ステップの変化が小さいときに、状態が進まなくなるのを防ぐためである。

修正の回数の上限は 20 回とする。

## 5 収束しないとき

20 回修正しても収束しないとき、ヤコビ行列の連立 1 次方程式が解けないとき、または 4.4 の直線探索で残差が小さくならないときは、計算を止める。

このとき、ステップ $n$ に相対湿度 $\varphi_{i,n}$ が 99.9% 以上の質点があれば、結露したもの（水分が入りすぎて、相対湿度 100% 以下では水分の収支が合わない）として、その質点を示して止める。結露した水を扱うモデルは、まだ扱わない。

## 6 根拠

未知数の並べ方、許容値の決め方、数値微分のずらす量、水分化学ポテンシャルの上限と直線探索、修正の回数の下限と上限、収束しないときに止める理由は、「{doc}`数値解法の安定性 <04_05_basis_solver>`」に示す。
