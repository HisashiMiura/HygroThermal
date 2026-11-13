# 熱・水分の収支

## 1 はじめに

### 1.1 概要

本節では、質点ごとの、ステップ $n$ から $n + 1$ までの熱・水分の収支の式を定める。収支は、各質点のコントロールボリュームでとる。時間は後退差分で離散化し、流れと外から入る熱・水分はステップ $n + 1$ の値を使う。ただし、換気量はステップ $n$ の値を使う（「{doc}`通気層の換気 <03_10_eval_ventilation>`」）。ここで定める収支の式を満たすステップ $n + 1$ の温度と水分化学ポテンシャルを、1 ステップの解き方で求める（「1 ステップの解き方」）。

### 1.2 入力

:::{container} symbols
$\Delta t$
: ：時間刻み（s）

$T_{i,n}$, $\mu_{i,n}$
: ：ステップ $n$ の質点 $i$ の温度（K）と水分化学ポテンシャル（J/kg）

$T_{i,n+1}$, $\mu_{i,n+1}$
: ：ステップ $n + 1$ の質点 $i$ の温度（K）と水分化学ポテンシャル（J/kg）

$w_i$
: ：質点 $i$ のコントロールボリュームの幅。「{doc}`壁の構成と質点 <03_04_eval_wall_nodes>`」で求める（m）

$q_{k,n+1}$, $j_{\mathrm{v},k,n+1}$, $j_{\mathrm{l},k,n+1}$
: ：ステップ $n + 1$ の、つながり $k$ の熱・水蒸気・液水の流れ。「{doc}`熱・水分の流れ <03_05_eval_flows>`」で、ステップ $n$ を $n + 1$ に読み替えて求める（W/m²、kg/(m² s)）

$q_{\mathrm{v},i,n}$
: ：ステップ $n$ の、質点 $i$ の壁 1 m² あたりの換気量。「{doc}`通気層の換気 <03_10_eval_ventilation>`」で求める（m³/(m² s)）

$q_{\mathrm{sol},n+1}$
: ：ステップ $n + 1$ の、外表面に外から直接入る熱。「{doc}`外表面の日射・夜間放射 <03_09_eval_surface_solar>`」で、ステップ $n$ を $n + 1$ に読み替えて求める（W/m²）

$T_{\mathrm{ex},n+1}$
: ：ステップ $n + 1$ の外気の温度（K）
:::

材料の乾燥密度 $\rho_{\mathrm{s}}$、平衡含水率 $u$、乾燥した材料の容積比熱 $c_{\mathrm{s}} \rho_{\mathrm{s}}$ は、「{doc}`材料の性質 <03_03_eval_materials>`」で求める。

### 1.3 出力

質点ごとの、熱の収支の式と水分の収支の式。

## 2 記号及び添え字

### 2.1 記号

:::{container} symbols
$C_{i,n}$
: ：ステップ $n$ の質点 $i$ の熱容量（J/(m² K)）

$c_{\mathrm{a}}$
: ：空気の比熱（J/(kg K)）

$c_{\mathrm{s}} \rho_{\mathrm{s}}$
: ：乾燥した材料の容積比熱（J/(m³ K)）

$c_{\mathrm{w}}$
: ：水の比熱（J/(kg K)）

$j_{\mathrm{l},k,n}$
: ：ステップ $n$ の、つながり $k$ の液水の流れ（kg/(m² s)）

$j_{\mathrm{v},k,n}$
: ：ステップ $n$ の、つながり $k$ の水蒸気の流れ（kg/(m² s)）

$j_{\mathrm{vent},i,n}$
: ：ステップ $n$ の、質点 $i$ に換気で入る水分（kg/(m² s)）

$L_{\mathrm{v}}$
: ：水の蒸発潜熱（J/kg）

$N$
: ：質点の数（-）

$q_{k,n}$
: ：ステップ $n$ の、つながり $k$ の熱の流れ（W/m²）

$q_{\mathrm{cond},i,n}$
: ：ステップ $n$ の、質点 $i$ の凝縮熱。負なら蒸発の熱（W/m²）

$q_{\mathrm{sol},n}$
: ：ステップ $n$ の、外表面に外から直接入る熱（W/m²）

$q_{\mathrm{sol},i,n}$
: ：ステップ $n$ の、質点 $i$ に外から直接入る熱（W/m²）

$q_{\mathrm{v},i,n}$
: ：ステップ $n$ の、質点 $i$ の壁 1 m² あたりの換気量（m³/(m² s)）

$q_{\mathrm{vent},i,n}$
: ：ステップ $n$ の、質点 $i$ に換気で入る熱（W/m²）

$T_{i,n}$
: ：ステップ $n$ の質点 $i$ の温度（K）

$T_{\mathrm{ex},n}$
: ：ステップ $n$ の外気の温度（K）

$u_{i,n}$
: ：ステップ $n$ の質点 $i$ の平衡含水率。 $\varphi_{i,n}$ から「{doc}`材料の性質 <03_03_eval_materials>`」で求める（kg/kg）

$W_{i,n}$
: ：ステップ $n$ の質点 $i$ の水分量（kg/m²）

$w_i$
: ：質点 $i$ のコントロールボリュームの幅（m）

$\Delta t$
: ：時間刻み（s）

$\delta_i$
: ：質点 $i$ が材料の質点なら 1、通気層の空気の質点なら 0（-）

$\mu_{i,n}$
: ：ステップ $n$ の質点 $i$ の水分化学ポテンシャル（J/kg）

$\rho_{\mathrm{a},i,n}$
: ：ステップ $n$ の質点 $i$ の空気の密度（kg/m³）

$\rho_{\mathrm{a,ex},n}$
: ：ステップ $n$ の外気の密度（kg/m³）

$\rho_{\mathrm{s}}$
: ：材料の乾燥密度（kg/m³）

$\rho_{\mathrm{v},i,n}$
: ：ステップ $n$ の質点 $i$ の水蒸気の密度（kg/m³）

$\rho_{\mathrm{v,ex},n}$
: ：ステップ $n$ の外気の水蒸気の密度（kg/m³）

$\varphi_{i,n}$
: ：ステップ $n$ の質点 $i$ の相対湿度。 $T_{i,n}$ と $\mu_{i,n}$ から、「{doc}`湿り空気と水の性質 <03_02_eval_air_water>`」の式 (3) で求める（%）
:::

### 2.2 添え字

:::{container} symbols
$i$
: ：質点の番号

$k$
: ：つながりの番号

$n$
: ：時刻のステップ
:::

## 3 水分量と熱容量

ステップ $n$ の材料の質点 $i$ の水分量 $W_{i,n}$ は、式 {eq}`bl-w-mat` により、通気層の空気の質点 $i$ の水分量は、式 {eq}`bl-w-air` により計算される。水蒸気の密度 $\rho_{\mathrm{v},i,n}$ は、 $T_{i,n}$ と $\mu_{i,n}$ から、「{doc}`湿り空気と水の性質 <03_02_eval_air_water>`」の式 (5) で求める。材料の質点の空隙の中の水蒸気は、吸着水・液水に比べてごくわずかなので、水分量に含めない。

$$
W_{i,n} = \rho_{\mathrm{s}} \thinspace u_{i,n} \thinspace w_i
$$ (bl-w-mat)

$$
W_{i,n} = \rho_{\mathrm{v},i,n} \thinspace w_i
$$ (bl-w-air)

ステップ $n$ の材料の質点 $i$ の熱容量 $C_{i,n}$ は、式 {eq}`bl-c-mat` により、通気層の空気の質点 $i$ の熱容量は、式 {eq}`bl-c-air` により計算される。空気の密度 $\rho_{\mathrm{a},i,n}$ は、 $T_{i,n}$ から、「{doc}`湿り空気と水の性質 <03_02_eval_air_water>`」の式 (6) で求める。通気層の空気の熱容量には、水蒸気の分を含めない。

$$
C_{i,n} = c_{\mathrm{s}} \rho_{\mathrm{s}} \thinspace w_i + c_{\mathrm{w}} W_{i,n}
$$ (bl-c-mat)

$$
C_{i,n} = c_{\mathrm{a}} \thinspace \rho_{\mathrm{a},i,n} \thinspace w_i
$$ (bl-c-air)

材料の質点の熱容量には、含まれる水の分を含める。移動する水が運ぶ顕熱は、量が非常に小さいので含めない。

## 4 外から入る熱・水分

### 4.1 換気で入る熱・水分

通気層の空気の質点 $i$ に、ステップ $n + 1$ に換気で入る水分 $j_{\mathrm{vent},i,n+1}$ と熱 $q_{\mathrm{vent},i,n+1}$ は、式 {eq}`bl-j-vent` と式 {eq}`bl-q-vent` により計算される。入ってくる空気は外気とする。外気の水蒸気の密度 $\rho_{\mathrm{v,ex},n+1}$ と外気の密度 $\rho_{\mathrm{a,ex},n+1}$ は、「{doc}`湿り空気と水の性質 <03_02_eval_air_water>`」の式 (5) と式 (6) で求める。材料の質点では、換気量が 0 なので、どちらも 0 である。

$$
j_{\mathrm{vent},i,n+1} = q_{\mathrm{v},i,n} \left( \rho_{\mathrm{v,ex},n+1} - \rho_{\mathrm{v},i,n+1} \right)
$$ (bl-j-vent)

$$
q_{\mathrm{vent},i,n+1} = c_{\mathrm{a}} \thinspace \rho_{\mathrm{a,ex},n+1} \thinspace q_{\mathrm{v},i,n} \left( T_{\mathrm{ex},n+1} - T_{i,n+1} \right)
$$ (bl-q-vent)

### 4.2 外から直接入る熱

外表面の質点（質点 0）には、外から直接入る熱（日射の吸収 − 夜間放射）が入る。ステップ $n + 1$ の質点 $i$ に外から直接入る熱 $q_{\mathrm{sol},i,n+1}$ は、質点 0 では $q_{\mathrm{sol},n+1}$、それ以外の質点では 0 とする。

### 4.3 凝縮熱

ステップ $n + 1$ の質点 $i$ の凝縮熱 $q_{\mathrm{cond},i,n+1}$ は、式 {eq}`bl-q-cond` により計算される。凝縮量は、質点に入ってくる水蒸気から出ていく水蒸気を引いたものとする（材料の質点では、水蒸気のままたまる量はごくわずかなので無視する）。通気層の空気の質点では、たまる水分はすべて水蒸気のままで相変化は起きないので、凝縮熱を 0 とする（ $\delta_i = 0$）。

$$
q_{\mathrm{cond},i,n+1} = \delta_i \thinspace L_{\mathrm{v}} \left( j_{\mathrm{v},i,n+1} - j_{\mathrm{v},i+1,n+1} \right)
$$ (bl-q-cond)

## 5 収支の式

質点 $i$ の左のつながりはつながり $i$、右のつながりはつながり $i + 1$ である（「{doc}`壁の構成と質点 <03_04_eval_wall_nodes>`」）。

質点 $i$ の水分の収支は、式 {eq}`bl-m` で表される。左辺はコントロールボリュームにたまる水分、右辺は左のつながりから入る水分、右のつながりから出ていく水分、換気で入る水分の和である。たまる量を水分量の変化で書くので、水分量が厳密に保存される。

$$
\frac{W_{i,n+1} - W_{i,n}}{\Delta t}
= \left( j_{\mathrm{v},i,n+1} + j_{\mathrm{l},i,n+1} \right) - \left( j_{\mathrm{v},i+1,n+1} + j_{\mathrm{l},i+1,n+1} \right) + j_{\mathrm{vent},i,n+1}
$$ (bl-m)

質点 $i$ の熱の収支は、式 {eq}`bl-h` で表される。熱容量はステップ $n + 1$ の値を使う。

$$
C_{i,n+1} \frac{T_{i,n+1} - T_{i,n}}{\Delta t}
= q_{i,n+1} - q_{i+1,n+1} + q_{\mathrm{cond},i,n+1} + q_{\mathrm{sol},i,n+1} + q_{\mathrm{vent},i,n+1}
$$ (bl-h)

式 {eq}`bl-m` と式 {eq}`bl-h` を、 $i = 0, 1, \ldots, N - 1$ のすべての質点について満たす $T_{i,n+1}$ と $\mu_{i,n+1}$ を求める。
