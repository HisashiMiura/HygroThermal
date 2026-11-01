# 外気の条件

## 1 はじめに

### 1.1 概要

本節では、1 時間ごと 1 年分の気象データから、ステップ $n$ の外気の条件（温度・相対湿度・日射・夜間放射）を求める方法を定める。温度と相対湿度は、外気側の境界の条件に使う（「熱・水分の収支」）。日射と夜間放射は、外表面に外から直接入る熱の計算に使う（「外表面の日射・夜間放射」）。太陽の位置も、同じ時刻 $t_n$ について、「太陽位置」により求めて出力に含める。

### 1.2 入力

:::{container} symbols
$\theta_{\mathrm{ex},h}$
: ：気象データの行 $h$ の外気の温度（deg.C）

$x_{\mathrm{ex},h}$
: ：気象データの行 $h$ の外気の絶対湿度（g/kg(DA)）

$P_{\mathrm{ex},h}$
: ：気象データの行 $h$ の気圧（Pa）

$\hat{I}_{\mathrm{dn},h}$
: ：気象データの行 $h$ に書かれた法線面直達日射量（W/m² または MJ/(m² h)）

$\hat{I}_{\mathrm{sky},h}$
: ：気象データの行 $h$ に書かれた水平面天空日射量（W/m² または MJ/(m² h)）

$\hat{R}_{\mathrm{eff},h}$
: ：気象データの行 $h$ に書かれた夜間放射量（W/m² または MJ/(m² h)）

$\phi_{\mathrm{lat}}$
: ：緯度。北緯を正とする（deg）

$\lambda_{\mathrm{lon}}$
: ：経度。東経を正とする（deg）

$t_n$
: ：ステップ $n$ の時刻。1 月 1 日 0 時からの経過時間（s）
:::

また、気象データの日射・夜間放射について、次の 3 つを与える。

- 単位：W/m² または MJ/(m² h)（3 節）
- 値の種類：瞬時値、または 1 時間の値（積算値・平均値）（5 節・6 節）
- 1 時間の区間：行 $h$ の値が表す 1 時間。前・中央・後のいずれか。値の種類が 1 時間の値のときだけ使う（6 節）

### 1.3 出力

:::{container} symbols
$\theta_{\mathrm{ex},n}$
: ：ステップ $n$ の外気の温度（deg.C）

$\varphi_{\mathrm{ex},n}$
: ：ステップ $n$ の外気の相対湿度（%）

$I_{\mathrm{dn},n}$
: ：ステップ $n$ の法線面直達日射量（W/m²）

$I_{\mathrm{sky},n}$
: ：ステップ $n$ の水平面天空日射量（W/m²）

$R_{\mathrm{eff},n}$
: ：ステップ $n$ の夜間放射量（W/m²）。負の値もある

$\sin h_{\mathrm{sun},n}$, $\cos h_{\mathrm{sun},n}$
: ：ステップ $n$ の太陽高度の正弦と余弦（-）

$\sin a_{\mathrm{sun},n}$, $\cos a_{\mathrm{sun},n}$
: ：ステップ $n$ の太陽方位角の正弦と余弦（-）
:::

## 2 記号及び添え字

### 2.1 記号

:::{container} symbols
$a_{\mathrm{sun},n}$
: ：ステップ $n$ の太陽方位角。南を 0、西回りを正とする（rad）

$f_n$
: ：ステップ $n$ の、前の行の時刻から測った 1 時間の中の割合（-）

$f_{\mathrm{u}}$
: ：日射・夜間放射を W/m² に換算する係数（-）

$h_{0,n}$, $h_{1,n}$
: ：ステップ $n$ の時刻の前の行と次の行の番号（-）

$h_n$
: ：ステップ $n$ の時刻を含む 1 時間の値が書かれた行の番号（-）

$h_{\mathrm{sun},n}$
: ：ステップ $n$ の太陽高度（rad）

$I_{\mathrm{dn},h}$
: ：行 $h$ の法線面直達日射量（W/m²）

$\hat{I}_{\mathrm{dn},h}$
: ：気象データの行 $h$ に書かれた法線面直達日射量（W/m² または MJ/(m² h)）

$I_{\mathrm{dn},n}$
: ：ステップ $n$ の法線面直達日射量（W/m²）

$I_{\mathrm{sky},h}$
: ：行 $h$ の水平面天空日射量（W/m²）

$\hat{I}_{\mathrm{sky},h}$
: ：気象データの行 $h$ に書かれた水平面天空日射量（W/m² または MJ/(m² h)）

$I_{\mathrm{sky},n}$
: ：ステップ $n$ の水平面天空日射量（W/m²）

$P_{\mathrm{ex},h}$
: ：行 $h$ の気圧（Pa）

$p_{\mathrm{v,ex},h}$
: ：行 $h$ の外気の水蒸気圧（Pa）

$p_{\mathrm{v,ex},n}$
: ：ステップ $n$ の外気の水蒸気圧（Pa）

$p_{\mathrm{vs}}(T)$
: ：温度 $T$ の飽和水蒸気圧（Pa）。「水分の性質」で定める

$R_{\mathrm{eff},h}$
: ：行 $h$ の夜間放射量（W/m²）

$\hat{R}_{\mathrm{eff},h}$
: ：気象データの行 $h$ に書かれた夜間放射量（W/m² または MJ/(m² h)）

$R_{\mathrm{eff},n}$
: ：ステップ $n$ の夜間放射量（W/m²）

$t_n$
: ：ステップ $n$ の時刻。1 月 1 日 0 時からの経過時間（s）

$t_{\mathrm{y},n}$
: ：ステップ $n$ の、1 年の中での時刻。1 月 1 日 0 時を 0 とする（h）

$x_{\mathrm{ex},h}$
: ：行 $h$ の外気の絶対湿度（g/kg(DA)）

$X_h$
: ：行 $h$ の値。 $\theta_{\mathrm{ex},h}$、 $p_{\mathrm{v,ex},h}$、 $I_{\mathrm{dn},h}$、 $I_{\mathrm{sky},h}$、 $R_{\mathrm{eff},h}$ のいずれか

$X_n$
: ：ステップ $n$ の値。 $X_h$ に対応する、ステップ $n$ の値

$\theta_{\mathrm{ex},h}$
: ：行 $h$ の外気の温度（deg.C）

$\theta_{\mathrm{ex},n}$
: ：ステップ $n$ の外気の温度（deg.C）

$\lambda_{\mathrm{lon}}$
: ：経度。東経を正とする（deg）

$\phi_{\mathrm{lat}}$
: ：緯度。北緯を正とする（deg）

$\varphi_{\mathrm{ex},n}$
: ：ステップ $n$ の外気の相対湿度（%）
:::

### 2.2 添え字

:::{container} symbols
$h$
: ：気象データの行。1 月 1 日 0 時の行を 0 とし、1 時間ごとに 1 つずつ増える（0, 1, …, 8759）

$n$
: ：時刻のステップ
:::

## 3 気象データの値

気象データは、1 月 1 日 0 時から 12 月 31 日 23 時まで、1 時間ごと 8760 行の値とする。行 $h$ は、1 月 1 日 0 時から $h$ 時間後の行である。

行 $h$ の外気の水蒸気圧 $p_{\mathrm{v,ex},h}$ は、式 {eq}`oc-pv` により計算される。相対湿度が 100 % を超える行は、水蒸気圧を飽和水蒸気圧とする。

$$
p_{\mathrm{v,ex},h} = \min \left( \frac{P_{\mathrm{ex},h} \, x_{\mathrm{ex},h}}{622 + x_{\mathrm{ex},h}}, \ p_{\mathrm{vs}} \left( \theta_{\mathrm{ex},h} + 273.15 \right) \right)
$$ (oc-pv)

行 $h$ の法線面直達日射量 $I_{\mathrm{dn},h}$、水平面天空日射量 $I_{\mathrm{sky},h}$、夜間放射量 $R_{\mathrm{eff},h}$ は、式 {eq}`oc-idn`〜式 {eq}`oc-reff` により計算される。

$$
I_{\mathrm{dn},h} = f_{\mathrm{u}} \hat{I}_{\mathrm{dn},h}
$$ (oc-idn)

$$
I_{\mathrm{sky},h} = f_{\mathrm{u}} \hat{I}_{\mathrm{sky},h}
$$ (oc-isky)

$$
R_{\mathrm{eff},h} = f_{\mathrm{u}} \hat{R}_{\mathrm{eff},h}
$$ (oc-reff)

換算の係数 $f_{\mathrm{u}}$ は、気象データの日射・夜間放射の単位に応じて次のとおりとする。

| 単位 | $f_{\mathrm{u}}$ |
|---|---|
| W/m²（平均値・瞬時値） | 1 |
| MJ/(m² h)（1 時間の積算値） | $10^6 / 3600$ |

## 4 1 年の中での時刻

1 年は 365 日とし、2 年目以降も同じ気象データを繰り返す。ステップ $n$ の、1 年の中での時刻 $t_{\mathrm{y},n}$ は、式 {eq}`oc-ty` により計算される。

$$
t_{\mathrm{y},n} = \frac{t_n \bmod (365 \times 86400)}{3600}
$$ (oc-ty)

## 5 瞬時値の補い方

瞬時値は、行 $h$ の値を $h$ 時ちょうどの値とし、前後の行の間を直線で補う。ステップ $n$ の時刻の前の行の番号 $h_{0,n}$ と次の行の番号 $h_{1,n}$、前の行の時刻から測った割合 $f_n$ は、式 {eq}`oc-h0`〜式 {eq}`oc-f` により計算される。12 月 31 日 23 時の行の次は、1 月 1 日 0 時の行とする。

$$
h_{0,n} = \left\lfloor t_{\mathrm{y},n} \right\rfloor
$$ (oc-h0)

$$
h_{1,n} = \left( h_{0,n} + 1 \right) \bmod 8760
$$ (oc-h1)

$$
f_n = t_{\mathrm{y},n} - h_{0,n}
$$ (oc-f)

瞬時値のステップ $n$ の値 $X_n$ は、式 {eq}`oc-instant` により計算される。

$$
X_n = X_{h_{0,n}} + f_n \left( X_{h_{1,n}} - X_{h_{0,n}} \right)
$$ (oc-instant)

## 6 1 時間の値の取り出し方

1 時間の値（1 時間の積算値・平均値）は、その 1 時間の中では一定とする。行 $h$ の値が表す 1 時間は、気象データに応じて、次の 3 つのいずれかとする。1 時間の区間は、始まりの時刻を含まず、終わりの時刻を含む。

| 区間 | 行 $h$ の値が表す 1 時間 |
|---|---|
| 前 | $h - 1$ 時 〜 $h$ 時 |
| 中央 | $h - 0.5$ 時 〜 $h + 0.5$ 時 |
| 後 | $h$ 時 〜 $h + 1$ 時 |

ステップ $n$ の時刻を含む 1 時間の値が書かれた行の番号 $h_n$ は、式 {eq}`oc-hi` により計算される。

$$
h_n =
\begin{cases}
\left\lceil t_{\mathrm{y},n} \right\rceil \bmod 8760 & (\text{前}) \\
\left\lceil t_{\mathrm{y},n} - 0.5 \right\rceil \bmod 8760 & (\text{中央}) \\
\left( \left\lceil t_{\mathrm{y},n} \right\rceil - 1 \right) \bmod 8760 & (\text{後})
\end{cases}
$$ (oc-hi)

1 時間の値のステップ $n$ の値 $X_n$ は、式 {eq}`oc-interval` により計算される。

$$
X_n = X_{h_n}
$$ (oc-interval)

## 7 ステップ n の外気の条件

ステップ $n$ の外気の温度 $\theta_{\mathrm{ex},n}$ と水蒸気圧 $p_{\mathrm{v,ex},n}$ は、瞬時値として、式 {eq}`oc-instant` により計算される（ $X_h$ をそれぞれ $\theta_{\mathrm{ex},h}$、 $p_{\mathrm{v,ex},h}$ とする）。

ステップ $n$ の外気の相対湿度 $\varphi_{\mathrm{ex},n}$ は、式 {eq}`oc-rh` により計算される。行の間を補った時刻でも 100 % を超えないようにする。

$$
\varphi_{\mathrm{ex},n} = \min \left( 100 \frac{p_{\mathrm{v,ex},n}}{p_{\mathrm{vs}} \left( \theta_{\mathrm{ex},n} + 273.15 \right)}, \ 100 \right)
$$ (oc-rh)

ステップ $n$ の法線面直達日射量 $I_{\mathrm{dn},n}$、水平面天空日射量 $I_{\mathrm{sky},n}$、夜間放射量 $R_{\mathrm{eff},n}$ は、値の種類が瞬時値のときは式 {eq}`oc-instant` により、1 時間の値のときは式 {eq}`oc-interval` により計算される（ $X_h$ をそれぞれ $I_{\mathrm{dn},h}$、 $I_{\mathrm{sky},h}$、 $R_{\mathrm{eff},h}$ とする）。

## 8 太陽の位置

ステップ $n$ の太陽高度の正弦 $\sin h_{\mathrm{sun},n}$ と余弦 $\cos h_{\mathrm{sun},n}$、太陽方位角の正弦 $\sin a_{\mathrm{sun},n}$ と余弦 $\cos a_{\mathrm{sun},n}$ は、緯度 $\phi_{\mathrm{lat}}$、経度 $\lambda_{\mathrm{lon}}$、時刻 $t_n$ から、「太陽位置」により計算される。

「太陽位置」には、1 年の中での時刻 $t_{\mathrm{y},n}$ ではなく、時刻 $t_n$ を与える。通算日とその日の標準時は、「太陽位置」の中で $t_n$ から求める。
