# 熱・水分の流れ

## 1 はじめに

### 1.1 概要

本節では、ステップ $n$ の、つながりごとの熱・水蒸気・液水の流れの求め方を定める。外表面・内表面の伝達、層と層の接触、材料と通気層の間の伝達、塗膜・防湿層などは、境目の抵抗として扱う。ここで求める流れは、熱・水分の収支に使う（「熱・水分の収支」）。

### 1.2 入力

:::{container} symbols
$G_{\mathrm{h},k}$
: ：境目の抵抗が入るつながり $k$ の、境目の熱コンダクタンス（W/(m² K)）

$G_{\mathrm{m},k}$
: ：境目の抵抗が入るつながり $k$ の、境目の湿気コンダクタンス（kg/(m² s Pa)）

$T_{i,n}$
: ：ステップ $n$ の質点 $i$ の温度（K）

$\mu_{i,n}$
: ：ステップ $n$ の質点 $i$ の水分化学ポテンシャル（J/kg）

$T_{\mathrm{ex},n}$, $T_{\mathrm{in},n}$
: ：ステップ $n$ の外気・室内の温度。「{doc}`外気の条件 <03_07_eval_outdoor>`」「{doc}`室内の条件 <03_08_eval_indoor>`」の温度に $T_0$ を足したもの（K）

$\mu_{\mathrm{ex},n}$, $\mu_{\mathrm{in},n}$
: ：ステップ $n$ の外気・室内の水分化学ポテンシャル。外気・室内の相対湿度と温度から、「{doc}`湿り空気と水の性質 <03_02_eval_air_water>`」の式 (2) で求める（J/kg）

$d_{\mathrm{L},i}$, $d_{\mathrm{R},i}$
: ：質点 $i$ からコントロールボリュームの左・右の端までの距離。「{doc}`壁の構成と質点 <03_04_eval_wall_nodes>`」で求める（m）
:::

材料の質点の熱伝導率 $\lambda_{\mathrm{h}}$、湿気伝導率 $\lambda'_{\mathrm{m}}$、液相水分伝導率 $\lambda'_{\mu \mathrm{l}}$ は、「{doc}`材料の性質 <03_03_eval_materials>`」で求める。

### 1.3 出力

:::{container} symbols
$q_{k,n}$
: ：ステップ $n$ の、つながり $k$ の熱の流れ。外気側から室内側へ流れる向きを正とする（W/m²）

$j_{\mathrm{v},k,n}$
: ：ステップ $n$ の、つながり $k$ の水蒸気の流れ。向きは熱と同じ（kg/(m² s)）

$j_{\mathrm{l},k,n}$
: ：ステップ $n$ の、つながり $k$ の液水の流れ。向きは熱と同じ（kg/(m² s)）
:::

## 2 記号及び添え字

### 2.1 記号

:::{container} symbols
$G_{\mathrm{h},k}$
: ：つながり $k$ の境目の熱コンダクタンス（W/(m² K)）

$G_{\mathrm{m},k}$
: ：つながり $k$ の境目の湿気コンダクタンス（kg/(m² s Pa)）

$j_{\mathrm{l},k,n}$
: ：ステップ $n$ の、つながり $k$ の液水の流れ（kg/(m² s)）

$j_{\mathrm{v},k,n}$
: ：ステップ $n$ の、つながり $k$ の水蒸気の流れ（kg/(m² s)）

$N$
: ：質点の数（-）

$p_{\mathrm{v},i,n}$
: ：ステップ $n$ の質点 $i$ の水蒸気圧（Pa）

$q_{k,n}$
: ：ステップ $n$ の、つながり $k$ の熱の流れ（W/m²）

$R_{\mathrm{h},k}$
: ：つながり $k$ の熱抵抗（m² K/W）

$R_{\mathrm{h,int},k}$
: ：つながり $k$ の境目の熱抵抗（m² K/W）

$R_{\mathrm{l},k,n}$
: ：ステップ $n$ の、つながり $k$ の液水の抵抗（m² s (J/kg)/kg）

$R_{\mathrm{m},k,n}$
: ：ステップ $n$ の、つながり $k$ の透湿抵抗（m² s Pa/kg）

$R_{\mathrm{m,int},k}$
: ：つながり $k$ の境目の透湿抵抗（m² s Pa/kg）

$r_{\mathrm{h,L},i}$, $r_{\mathrm{h,R},i}$
: ：質点 $i$ からコントロールボリュームの左・右の端までの熱抵抗（m² K/W）

$r_{\mathrm{l,L},i,n}$, $r_{\mathrm{l,R},i,n}$
: ：ステップ $n$ の、質点 $i$ からコントロールボリュームの左・右の端までの液水の抵抗（m² s (J/kg)/kg）

$r_{\mathrm{m,L},i,n}$, $r_{\mathrm{m,R},i,n}$
: ：ステップ $n$ の、質点 $i$ からコントロールボリュームの左・右の端までの透湿抵抗（m² s Pa/kg）

$T_{i,n}$
: ：ステップ $n$ の質点 $i$ の温度（K）

$\lambda_{\mathrm{h}}$
: ：材料の熱伝導率（W/(m K)）

$\lambda'_{\mathrm{m}}$
: ：材料の湿気伝導率（kg/(m s Pa)）

$\lambda'_{\mu \mathrm{l}}$
: ：材料の液相水分伝導率（kg/(m s (J/kg))）

$\mu_{i,n}$
: ：ステップ $n$ の質点 $i$ の水分化学ポテンシャル（J/kg）

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

## 3 境目の抵抗

境目の抵抗が入るつながり（つながり 0、層の変わるつながり、つながり $N$）の、境目の熱抵抗 $R_{\mathrm{h,int},k}$ と透湿抵抗 $R_{\mathrm{m,int},k}$ は、式 {eq}`fl-rh-int` と式 {eq}`fl-rm-int` により計算される。同じ層の中のつながりでは、境目の抵抗は 0 とする。

$$
R_{\mathrm{h,int},k} = \frac{1}{G_{\mathrm{h},k}}
$$ (fl-rh-int)

$$
R_{\mathrm{m,int},k} = \frac{1}{G_{\mathrm{m},k}}
$$ (fl-rm-int)

境目では、液水は移動しないものとする。

## 4 質点から端までの抵抗

材料の質点 $i$ から、コントロールボリュームの左の端までの熱抵抗 $r_{\mathrm{h,L},i}$、透湿抵抗 $r_{\mathrm{m,L},i,n}$、液水の抵抗 $r_{\mathrm{l,L},i,n}$ は、式 {eq}`fl-rh-half`〜式 {eq}`fl-rl-half` により計算される。右の端までの抵抗 $r_{\mathrm{h,R},i}$・ $r_{\mathrm{m,R},i,n}$・ $r_{\mathrm{l,R},i,n}$ も、 $d_{\mathrm{L},i}$ を $d_{\mathrm{R},i}$ に替えた同じ式で計算される。湿気伝導率と液相水分伝導率は、その質点のステップ $n$ の状態で求める。

$$
r_{\mathrm{h,L},i} = \frac{d_{\mathrm{L},i}}{\lambda_{\mathrm{h}}}
$$ (fl-rh-half)

$$
r_{\mathrm{m,L},i,n} = \frac{d_{\mathrm{L},i}}{\lambda'_{\mathrm{m}} \left( \varphi_{i,n} \right)}
$$ (fl-rm-half)

$$
r_{\mathrm{l,L},i,n} = \frac{d_{\mathrm{L},i}}{\lambda'_{\mu \mathrm{l}} \left( \varphi_{i,n}, T_{i,n} \right)}
$$ (fl-rl-half)

液相水分伝導率が 0 の材料（液水が移動しない材料）では、液水の抵抗は無限大とする。

通気層の質点は、よく混ざった空気とみなし、端までの熱抵抗と透湿抵抗を 0 とする。

## 5 つながりの抵抗

つながり $k$ の熱抵抗 $R_{\mathrm{h},k}$、透湿抵抗 $R_{\mathrm{m},k,n}$、液水の抵抗 $R_{\mathrm{l},k,n}$ は、左の質点 $k - 1$ から右の端までの抵抗、境目の抵抗、右の質点 $k$ から左の端までの抵抗の和として、式 {eq}`fl-rh`〜式 {eq}`fl-rl` により計算される。外気（ $k = 0$ の左）と室内（ $k = N$ の右）の側は、端までの抵抗を 0 とする。

$$
R_{\mathrm{h},k} = r_{\mathrm{h,R},k-1} + R_{\mathrm{h,int},k} + r_{\mathrm{h,L},k}
$$ (fl-rh)

$$
R_{\mathrm{m},k,n} = r_{\mathrm{m,R},k-1,n} + R_{\mathrm{m,int},k} + r_{\mathrm{m,L},k,n}
$$ (fl-rm)

$$
R_{\mathrm{l},k,n} =
\begin{cases}
r_{\mathrm{l,R},k-1,n} + r_{\mathrm{l,L},k,n} & (\text{同じ層の中のつながり}) \\
\infty & (\text{境目の抵抗が入るつながり})
\end{cases}
$$ (fl-rl)

## 6 流れ

外気と室内の値は、質点の番号の外側の $-1$ と $N$ に置き、 $T_{-1,n} = T_{\mathrm{ex},n}$、 $T_{N,n} = T_{\mathrm{in},n}$ とする（ $\mu$ と $p_{\mathrm{v}}$ も同じ）。

ステップ $n$ の、つながり $k$ の熱の流れ $q_{k,n}$ は、式 {eq}`fl-q` により計算される。

$$
q_{k,n} = \frac{T_{k-1,n} - T_{k,n}}{R_{\mathrm{h},k}}
$$ (fl-q)

ステップ $n$ の、つながり $k$ の水蒸気の流れ $j_{\mathrm{v},k,n}$ は、式 {eq}`fl-jv` により計算される。水蒸気圧 $p_{\mathrm{v},i,n}$ は、 $T_{i,n}$ と $\varphi_{i,n}$ から、「{doc}`湿り空気と水の性質 <03_02_eval_air_water>`」の式 (4) で求める。

$$
j_{\mathrm{v},k,n} = \frac{p_{\mathrm{v},k-1,n} - p_{\mathrm{v},k,n}}{R_{\mathrm{m},k,n}}
$$ (fl-jv)

ステップ $n$ の、つながり $k$ の液水の流れ $j_{\mathrm{l},k,n}$ は、式 {eq}`fl-jl` により計算される。液水の抵抗が無限大のつながりでは、液水の流れは 0 である。

$$
j_{\mathrm{l},k,n} = \frac{\mu_{k-1,n} - \mu_{k,n}}{R_{\mathrm{l},k,n}}
$$ (fl-jl)
