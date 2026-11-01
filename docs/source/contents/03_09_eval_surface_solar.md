# 外表面の日射・夜間放射

## 1 はじめに

### 1.1 概要

本節では、壁の外表面に外から直接入る熱（日射の吸収から夜間放射を引いたもの）の計算方法を定める。この熱は、外表面の質点の熱の収支に入れる（「熱・水分の収支」）。

### 1.2 入力

:::{container} symbols
$\alpha_{\mathrm{w}}$
: ：外表面の方位角。外表面の方位から 3 節で定める（deg）

$\beta_{\mathrm{w}}$
: ：外表面の傾斜角（deg）

$\rho_{\mathrm{g}}$
: ：地面の日射反射率（アルベド）（-）

$a_{\mathrm{s}}$
: ：外表面の日射吸収率（-）

$\varepsilon$
: ：外表面の長波長放射率（-）

$\sin h_{\mathrm{sun},n}$, $\cos h_{\mathrm{sun},n}$
: ：ステップ $n$ の太陽高度の正弦と余弦。「{doc}`外気の条件 <03_07_eval_outdoor>`」で求める（-）

$\sin a_{\mathrm{sun},n}$, $\cos a_{\mathrm{sun},n}$
: ：ステップ $n$ の太陽方位角の正弦と余弦。「{doc}`外気の条件 <03_07_eval_outdoor>`」で求める（-）

$I_{\mathrm{dn},n}$
: ：ステップ $n$ の法線面直達日射量。「{doc}`外気の条件 <03_07_eval_outdoor>`」で求める（W/m²）

$I_{\mathrm{sky},n}$
: ：ステップ $n$ の水平面天空日射量。「{doc}`外気の条件 <03_07_eval_outdoor>`」で求める（W/m²）

$R_{\mathrm{eff},n}$
: ：ステップ $n$ の夜間放射量。「{doc}`外気の条件 <03_07_eval_outdoor>`」で求める（W/m²）
:::

### 1.3 出力

:::{container} symbols
$q_{\mathrm{sol},n}$
: ：ステップ $n$ の、外表面に外から直接入る熱（W/m²）
:::

## 2 記号及び添え字

### 2.1 記号

:::{container} symbols
$a_{\mathrm{sun},n}$
: ：ステップ $n$ の太陽方位角。南を 0、西回りを正とする（rad）

$a_{\mathrm{s}}$
: ：外表面の日射吸収率（-）

$F_{\mathrm{sky}}$
: ：外表面の天空に対する形態係数（-）

$h_{\mathrm{sun},n}$
: ：ステップ $n$ の太陽高度（rad）

$I_{\mathrm{D},n}$
: ：ステップ $n$ の、外表面に当たる直達日射量（W/m²）

$I_{\mathrm{dn},n}$
: ：ステップ $n$ の法線面直達日射量（W/m²）

$I_{\mathrm{R},n}$
: ：ステップ $n$ の、外表面に当たる地面の反射日射量（W/m²）

$I_{\mathrm{S},n}$
: ：ステップ $n$ の、外表面に当たる天空日射量（W/m²）

$I_{\mathrm{sky},n}$
: ：ステップ $n$ の水平面天空日射量（W/m²）

$I_{\mathrm{T},n}$
: ：ステップ $n$ の水平面全天日射量（W/m²）

$q_{\mathrm{sol},n}$
: ：ステップ $n$ の、外表面に外から直接入る熱（W/m²）

$R_{\mathrm{eff},n}$
: ：ステップ $n$ の夜間放射量（水平面の実効放射量。地表面からの上向きの長波長放射と、下向きの大気放射の差）（W/m²）

$\alpha_{\mathrm{w}}$
: ：外表面の方位角。南を 0、西回りを正とする（deg）

$\beta_{\mathrm{w}}$
: ：外表面の傾斜角。上向きの水平を 0、垂直を 90、下向きの水平を 180 とする（deg）

$\varepsilon$
: ：外表面の長波長放射率（-）

$\rho_{\mathrm{g}}$
: ：地面の日射反射率（アルベド）（-）

$\psi_n$
: ：ステップ $n$ の、外表面への日射の入射角（rad）
:::

### 2.2 添え字

:::{container} symbols
$n$
: ：時刻のステップ
:::

## 3 外表面の方位角

外表面の方位角 $\alpha_{\mathrm{w}}$ は、外表面の方位に応じて次のとおりとする。

| 方位 | 南 | 南西 | 西 | 北西 | 北 | 北東 | 東 | 南東 |
|---|---|---|---|---|---|---|---|---|
| 方位角, deg | 0 | 45 | 90 | 135 | 180 | 225 | 270 | 315 |

上向きの水平面（傾斜角 0 deg）と下向きの水平面（傾斜角 180 deg）は方位を持たないので、方位角を 0 deg とする。このとき $\sin \beta_{\mathrm{w}} = 0$ なので、方位角は計算の結果に影響しない。

式の三角関数に入れる方位角と傾斜角は、rad に換算した値とする。

## 4 外表面に当たる日射量

ステップ $n$ の、外表面への日射の入射角の余弦 $\cos \psi_n$ は、式 {eq}`ss-cos` により計算される。

$$
\cos \psi_n = \cos \beta_{\mathrm{w}} \sin h_{\mathrm{sun},n} + \sin \beta_{\mathrm{w}} \cos h_{\mathrm{sun},n} \left( \cos a_{\mathrm{sun},n} \cos \alpha_{\mathrm{w}} + \sin a_{\mathrm{sun},n} \sin \alpha_{\mathrm{w}} \right)
$$ (ss-cos)

ステップ $n$ の、外表面に当たる直達日射量 $I_{\mathrm{D},n}$ は、式 {eq}`ss-direct` により計算される。太陽が地平線より下にある場合（ $\sin h_{\mathrm{sun},n} \le 0$ ）と、太陽が外表面の裏側にある場合（ $\cos \psi_n \le 0$ ）は 0 とする。

$$
I_{\mathrm{D},n} =
\begin{cases}
I_{\mathrm{dn},n} \cos \psi_n & (\sin h_{\mathrm{sun},n} > 0 \ \text{かつ} \ \cos \psi_n > 0) \\
0 & (\text{それ以外})
\end{cases}
$$ (ss-direct)

なお、垂直の面では、太陽が地平線のすぐ下にあっても $\cos \psi_n$ が正になる。そのため、入射角だけでなく太陽高度も確かめる。

外表面の天空に対する形態係数 $F_{\mathrm{sky}}$ は、式 {eq}`ss-fsky` により計算される。

$$
F_{\mathrm{sky}} = \frac{1 + \cos \beta_{\mathrm{w}}}{2}
$$ (ss-fsky)

ステップ $n$ の、外表面に当たる天空日射量 $I_{\mathrm{S},n}$ は、式 {eq}`ss-sky` により計算される。

$$
I_{\mathrm{S},n} = F_{\mathrm{sky}} I_{\mathrm{sky},n}
$$ (ss-sky)

ステップ $n$ の水平面全天日射量 $I_{\mathrm{T},n}$ は、式 {eq}`ss-total` により計算される。直達日射の分は、太陽が地平線より上にあるときだけ足す。

$$
I_{\mathrm{T},n} = I_{\mathrm{dn},n} \max \left( \sin h_{\mathrm{sun},n}, 0 \right) + I_{\mathrm{sky},n}
$$ (ss-total)

ステップ $n$ の、外表面に当たる地面の反射日射量 $I_{\mathrm{R},n}$ は、式 {eq}`ss-reflect` により計算される。

$$
I_{\mathrm{R},n} = \rho_{\mathrm{g}} \left( 1 - F_{\mathrm{sky}} \right) I_{\mathrm{T},n}
$$ (ss-reflect)

下向きの面（ $\beta_{\mathrm{w}} = 180$ deg）では、 $F_{\mathrm{sky}} = 0$ なので、直達日射量と天空日射量は 0 になり、地面の反射日射量だけが当たる。

## 5 外表面に外から直接入る熱

ステップ $n$ の、外表面に外から直接入る熱 $q_{\mathrm{sol},n}$ は、式 {eq}`ss-qsol` により計算される。

$$
q_{\mathrm{sol},n} = a_{\mathrm{s}} \left( I_{\mathrm{D},n} + I_{\mathrm{S},n} + I_{\mathrm{R},n} \right) - \varepsilon F_{\mathrm{sky}} R_{\mathrm{eff},n}
$$ (ss-qsol)

右辺の第 1 項は外表面が吸収する日射、第 2 項は外表面から天空へ出ていく夜間放射である。
