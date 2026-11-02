# 湿り空気と水の性質

## 1 はじめに

### 1.1 概要

本節では、温度・相対湿度・水蒸気圧・水分化学ポテンシャルの関係と、計算に使う空気と水の値を定める。ここで定める式と値は、ほかの節で共通に使う。

### 1.2 入力

:::{container} symbols
$T$
: ：温度（K）

$\varphi$
: ：相対湿度（%）。0 より大きく 100 以下とする

$\mu$
: ：水分化学ポテンシャル（J/kg）
:::

ある状態の温度と、相対湿度または水分化学ポテンシャルのどちらか一方を与える。

### 1.3 出力

:::{container} symbols
$p_{\mathrm{vs}}$
: ：飽和水蒸気圧（Pa）

$\mu$
: ：水分化学ポテンシャル（J/kg）。相対湿度を与えたとき

$\varphi$
: ：相対湿度（%）。水分化学ポテンシャルを与えたとき

$p_{\mathrm{v}}$
: ：水蒸気圧（Pa）

$\rho_{\mathrm{v}}$
: ：水蒸気の密度（kg/m³）

$\rho_{\mathrm{a}}$
: ：空気の密度（kg/m³）
:::

## 2 記号及び添え字

### 2.1 記号

:::{container} symbols
$c_{\mathrm{a}}$
: ：空気の比熱（J/(kg K)）

$c_{\mathrm{w}}$
: ：水の比熱（J/(kg K)）

$L_{\mathrm{v}}$
: ：水の蒸発潜熱（J/kg）

$p_{\mathrm{v}}$
: ：水蒸気圧（Pa）

$p_{\mathrm{vs}}$
: ：飽和水蒸気圧（Pa）

$R_{\mathrm{v}}$
: ：水蒸気の気体定数（J/(kg K)）

$T$
: ：温度（K）

$\mu$
: ：水分化学ポテンシャル（J/kg）

$\rho_{\mathrm{a}}$
: ：空気の密度（kg/m³）

$\rho_{\mathrm{v}}$
: ：水蒸気の密度（kg/m³）

$\rho_{\mathrm{w}}$
: ：水の密度（kg/m³）

$\varphi$
: ：相対湿度（%）
:::

### 2.2 添え字

本節では添え字を使わない。

## 3 空気と水の値

計算に使う空気と水の値は、次のとおりとする。

| 記号 | 意味 | 値 |
|---|---|---|
| $R_{\mathrm{v}}$ | 水蒸気の気体定数 | 461.5 J/(kg K) |
| $L_{\mathrm{v}}$ | 水の蒸発潜熱 | 2.512 × 10⁶ J/kg |
| $c_{\mathrm{w}}$ | 水の比熱 | 4200 J/(kg K) |
| $\rho_{\mathrm{w}}$ | 水の密度 | 998 kg/m³ |
| $c_{\mathrm{a}}$ | 空気の比熱 | 1005 J/(kg K) |

## 4 飽和水蒸気圧

飽和水蒸気圧 $p_{\mathrm{vs}}$ は、温度 $T$ の関数として、式 {eq}`aw-pvs` により計算される。

$$
p_{\mathrm{vs}} = \exp \left( - \frac{6096.938}{T} + 21.2409642 - 2.711193 \times 10^{-2} \, T + 1.673952 \times 10^{-5} \, T^2 + 2.433502 \ln T \right)
$$ (aw-pvs)

ほかの節では、温度 $T$ の飽和水蒸気圧を $p_{\mathrm{vs}}(T)$ と書く。

## 5 水分化学ポテンシャルと相対湿度

相対湿度 $\varphi$ を与えたときの水分化学ポテンシャル $\mu$ は、式 {eq}`aw-mu` により計算される。

$$
\mu = R_{\mathrm{v}} T \ln \frac{\varphi}{100}
$$ (aw-mu)

水分化学ポテンシャル $\mu$ を与えたときの相対湿度 $\varphi$ は、式 {eq}`aw-rh` により計算される。 $\mu > 0$ のときも、相対湿度は 100 % を超えないものとする。

$$
\varphi = \min \left( 100 \exp \frac{\mu}{R_{\mathrm{v}} T}, \ 100 \right)
$$ (aw-rh)

## 6 水蒸気圧と水蒸気の密度

水蒸気圧 $p_{\mathrm{v}}$ は、式 {eq}`aw-pv` により計算される。

$$
p_{\mathrm{v}} = \frac{\varphi}{100} \, p_{\mathrm{vs}}(T)
$$ (aw-pv)

水蒸気の密度 $\rho_{\mathrm{v}}$ は、式 {eq}`aw-rhov` により計算される。

$$
\rho_{\mathrm{v}} = \frac{p_{\mathrm{v}}}{R_{\mathrm{v}} T}
$$ (aw-rhov)

## 7 空気の密度

空気の密度 $\rho_{\mathrm{a}}$ は、温度 $T$ の関数として、式 {eq}`aw-rhoa` により計算される。空気の容積比熱は、空気の比熱 $c_{\mathrm{a}}$ とこの密度の積 $c_{\mathrm{a}} \rho_{\mathrm{a}}$ とする。

$$
\rho_{\mathrm{a}} = \frac{353}{T}
$$ (aw-rhoa)
