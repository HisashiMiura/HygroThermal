# 室内の条件

## 1 はじめに

### 1.1 概要

本節では、ステップ $n$ の室内の条件（温度・相対湿度）を求める方法を定める。室温は 1 年周期の余弦で変わり、相対湿度は一定とする。室内の条件は、室内側の境界の条件に使う（「熱・水分の収支」）。

### 1.2 入力

:::{container} symbols
$\theta_{\mathrm{in,ave}}$
: ：年平均の室温（deg.C）

$\theta_{\mathrm{in,amp}}$
: ：室温の振幅。0 以上とする（K）

$d_{\mathrm{peak}}$
: ：室温がいちばん高い日の通算日。1 月 1 日を 1 とし、1〜365 の整数とする（d）

$\varphi_{\mathrm{in}}$
: ：室内の相対湿度。0 より大きく 100 以下とする（%）

$t_n$
: ：ステップ $n$ の時刻。1 月 1 日 0 時からの経過時間（s）
:::

### 1.3 出力

:::{container} symbols
$\theta_{\mathrm{in},n}$
: ：ステップ $n$ の室温（deg.C）

$\varphi_{\mathrm{in},n}$
: ：ステップ $n$ の室内の相対湿度（%）
:::

## 2 記号及び添え字

### 2.1 記号

:::{container} symbols
$d_{\mathrm{peak}}$
: ：室温がいちばん高い日の通算日。1 月 1 日を 1 とする（d）

$t_n$
: ：ステップ $n$ の時刻。1 月 1 日 0 時からの経過時間（s）

$t_{\mathrm{peak}}$
: ：室温がいちばん高い時刻。1 月 1 日 0 時からの経過時間（s）

$\theta_{\mathrm{in,amp}}$
: ：室温の振幅（K）

$\theta_{\mathrm{in,ave}}$
: ：年平均の室温（deg.C）

$\theta_{\mathrm{in},n}$
: ：ステップ $n$ の室温（deg.C）

$\varphi_{\mathrm{in}}$
: ：室内の相対湿度（%）

$\varphi_{\mathrm{in},n}$
: ：ステップ $n$ の室内の相対湿度（%）
:::

### 2.2 添え字

:::{container} symbols
$n$
: ：時刻のステップ
:::

## 3 室温

室温がいちばん高い時刻 $t_{\mathrm{peak}}$ は、通算日 $d_{\mathrm{peak}}$ の日の 12 時とし、式 {eq}`ic-peak` により計算される。

$$
t_{\mathrm{peak}} = \left( d_{\mathrm{peak}} - 1 + 0.5 \right) \times 86400
$$ (ic-peak)

ステップ $n$ の室温 $\theta_{\mathrm{in},n}$ は、式 {eq}`ic-theta` により計算される。1 年は 365 日とする。余弦は 1 年で元に戻るので、2 年目以降も、 $t_n$ を 1 年で割った余りをとらずにそのまま使ってよい。

$$
\theta_{\mathrm{in},n} = \theta_{\mathrm{in,ave}} + \theta_{\mathrm{in,amp}} \cos \left( 2 \pi \frac{t_n - t_{\mathrm{peak}}}{365 \times 86400} \right)
$$ (ic-theta)

## 4 相対湿度

ステップ $n$ の室内の相対湿度 $\varphi_{\mathrm{in},n}$ は、式 {eq}`ic-rh` により計算される。

$$
\varphi_{\mathrm{in},n} = \varphi_{\mathrm{in}}
$$ (ic-rh)
