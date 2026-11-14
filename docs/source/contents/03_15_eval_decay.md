# 木材腐朽

## 1 はじめに

### 1.1 概要

本節では、木材腐朽を求める層の材料の質点の、ステップ $n$ から $n + 1$ までの木材腐朽の状態（相対湿度が高い状態が続いた時間、発芽したか、質量減少率）の進め方と、腐朽で生まれる水の求め方を定める。腐朽で生まれる水は、次のステップの水分の収支に加える（「{doc}`熱・水分の収支 <03_11_eval_balance>`」の 5 節）。

木材腐朽は、腐朽が始まるまで（発芽期）と、腐朽が進むあいだ（成長期）に分けて求める。腐朽が始まるのは層の表面の質点だけとし、層の内部へは、隣の質点で腐朽が進んでいるときに広がるものとする。

### 1.2 入力

:::{container} symbols
$T_{i,n+1}$
: ：ステップ $n + 1$ の質点 $i$ の温度。「{doc}`1 ステップの解き方 <03_12_eval_solve>`」で求める（K）

$\varphi_{i,n+1}$
: ：ステップ $n + 1$ の質点 $i$ の相対湿度。 $T_{i,n+1}$ と、質点 $i$ の水分化学ポテンシャルから、「{doc}`湿り空気と水の性質 <03_02_eval_air_water>`」の式 (3) で求める（%）

$\Delta t$
: ：時間刻み（s）

$w_i$
: ：質点 $i$ のコントロールボリュームの幅。「{doc}`壁の構成と質点 <03_04_eval_wall_nodes>`」で求める（m）

$\rho_{\mathrm{w}}$
: ：水の密度。「{doc}`湿り空気と水の性質 <03_02_eval_air_water>`」の値（kg/m³）
:::

木材腐朽を求めるかは、材料の層ごとに与える。

### 1.3 出力

:::{container} symbols
$L_{i,n+1}$
: ：ステップ $n + 1$ の質点 $i$ の質量減少率（-）

$j_{\mathrm{gen},i,n+2}$
: ：ステップ $n + 2$ の水分の収支に加える、質点 $i$ の腐朽で生まれる水（kg/(m² s)）
:::

## 2 記号及び添え字

### 2.1 記号

:::{container} symbols
$f$
: ：質量減少率の式の、時間に比例する部分の係数（1/月）

$g$
: ：質量減少率の式の、時間によらない部分（-）

$G_{i,n}$
: ：ステップ $n$ に、質点 $i$ が発芽しているか（真偽）

$h$
: ：質量減少率に対する水分生成率。 $h = 0.319$（-）

$j_{\mathrm{gen},i,n}$
: ：ステップ $n$ の水分の収支に加える、質点 $i$ の腐朽で生まれる水（kg/(m² s)）

$k_{\mathrm{m}}$
: ：腐朽の反応速度定数（1/s）

$L_{i,n}$
: ：ステップ $n$ の質点 $i$ の質量減少率（-）

$L_{\max}$
: ：質量減少率の上限。 $L_{\max} = 0.6$（-）

$t_{\mathrm{i}}$
: ：発芽までの時間（s）

$t'_{\mathrm{i}}$
: ：質量減少率の式の値が 0 になる時間（月）

$t_{\mathrm{i,min}}$
: ：発芽までの時間の下限。 $t_{\mathrm{i,min}} = 0.5$（月）

$t_{\mathrm{i,neg}}$
: ： $t'_{\mathrm{i}}$ が負のとき、0 とみなす範囲の下端。 $t_{\mathrm{i,neg}} = -0.59$（月）

$t_{\mathrm{i,never}}$
: ：発芽しないとみなすときの、発芽までの時間。 $t_{\mathrm{i,never}} = 100$（月）

$t_{\mathrm{mon}}$
: ：1 か月の長さ。 $t_{\mathrm{mon}} = 30 \times 86400$（s/月）

$WS$
: ：樹種と表面の状態を表す係数。 $WS = -1$（-）

$\Delta t$
: ：時間刻み（s）

$\theta$
: ：温度（deg.C）

$\theta_{i,n}$
: ：ステップ $n$ の質点 $i$ の温度。 $\theta_{i,n} = T_{i,n} - 273.15$（deg.C）

$\theta_{\mathrm{g,min}}$、 $\theta_{\mathrm{g,max}}$
: ：腐朽が進む温度の下限と上限。 $\theta_{\mathrm{g,min}} = 0$、 $\theta_{\mathrm{g,max}} = 40$（deg.C）

$\rho_{\mathrm{w}}$
: ：水の密度（kg/m³）

$\tau_{i,n}$
: ：ステップ $n$ までに、質点 $i$ の相対湿度が $\varphi_{\mathrm{c}}$ 以上の状態が続いた時間（s）

$\varphi$
: ：相対湿度（%）

$\varphi_{i,n}$
: ：ステップ $n$ の質点 $i$ の相対湿度（%）

$\varphi_{\mathrm{c}}$
: ：腐朽が始まる相対湿度の下限（%）

$\varphi_{\mathrm{g}}$
: ：腐朽が進む相対湿度の下限。 $\varphi_{\mathrm{g}} = 98$（%）

$w_i$
: ：質点 $i$ のコントロールボリュームの幅（m）
:::

### 2.2 添え字

:::{container} symbols
$i$
: ：質点の番号

$n$
: ：時刻のステップ
:::

## 3 木材腐朽を求める質点

木材腐朽を求める層の、すべての材料の質点で求める。そのうち、層の両端の質点（外気側の端と室内側の端）を、層の表面の質点とする。質点 $i$ の隣の質点は、同じ層の質点 $i - 1$ と $i + 1$ とする（層の端の質点は、層の内側の 1 つだけ）。

はじめ（ステップ 0）は、どの質点も腐朽していないものとする（ $\tau_{i,0} = 0$、 $G_{i,0}$ は偽、 $L_{i,0} = 0$）。

## 4 発芽期

発芽期は、層の表面の質点だけで求める。表面の質点 $i$ が発芽していないとき、ステップ $n + 1$ の温度 $\theta_{i,n+1}$ と相対湿度 $\varphi_{i,n+1}$ から、次のように進める。

腐朽が始まる相対湿度の下限 $\varphi_{\mathrm{c}}$ は、式 {eq}`dc-rhc` により計算される。

$$
\varphi_{\mathrm{c}} =
\begin{cases}
\min \left( -0.5 \thinspace \theta + 100, \thinspace 100 \right) & \left( \theta \le 15 \right) \\
92.5 & \left( \theta > 15 \right)
\end{cases}
$$ (dc-rhc)

発芽までの時間を求める式の係数 $f$ と $g$ は、式 {eq}`dc-f` と式 {eq}`dc-g` により計算される。

$$
f = 0.1384 \thinspace \theta + 0.4370 \thinspace \varphi - 42.9450 + WS \left( 0.0340 \thinspace \theta - 0.0210 \thinspace \varphi + 1.7210 \right)
$$ (dc-f)

$$
g = -2.2270 \thinspace \theta - 0.0347 \thinspace \varphi + 0.0244 \thinspace \theta \thinspace \varphi + WS \left( -0.504 \thinspace \theta + 0.0096 \thinspace \varphi + 0.0047 \thinspace \theta \thinspace \varphi \right)
$$ (dc-g)

質量減少率の式の値が 0 になる時間 $t'_{\mathrm{i}}$ は、式 {eq}`dc-ti-dash` により計算される。 $f = 0$ のときは、 $t'_{\mathrm{i}} = -\infty$ とする。

$$
t'_{\mathrm{i}} = - \frac{g}{f}
$$ (dc-ti-dash)

発芽までの時間 $t_{\mathrm{i}}$ は、式 {eq}`dc-ti` により計算される。

$$
t_{\mathrm{i}} = t_{\mathrm{mon}} \max \left( t_{\mathrm{i,min}}, \thinspace
\begin{cases}
t'_{\mathrm{i}} & \left( t'_{\mathrm{i}} > 0 \right) \\
0 & \left( t_{\mathrm{i,neg}} < t'_{\mathrm{i}} \le 0 \right) \\
t_{\mathrm{i,never}} & \left( t'_{\mathrm{i}} \le t_{\mathrm{i,neg}} \right)
\end{cases}
\right)
$$ (dc-ti)

ステップ $n + 1$ までに、相対湿度が $\varphi_{\mathrm{c}}$ 以上の状態が続いた時間 $\tau_{i,n+1}$ は、 $\theta = \theta_{i,n+1}$ の $\varphi_{\mathrm{c}}$ を使って、式 {eq}`dc-tau` により計算される。0 deg.C 以下のときと、相対湿度が $\varphi_{\mathrm{c}}$ より低いときは、0 に戻す。

$$
\tau_{i,n+1} =
\begin{cases}
0 & \left( \theta_{i,n+1} \le 0 \quad \text{または} \quad \varphi_{i,n+1} \lt \varphi_{\mathrm{c}} \right) \\
\tau_{i,n} + \Delta t & (\text{それ以外})
\end{cases}
$$ (dc-tau)

ステップ $n + 1$ に発芽しているか $G_{i,n+1}$ は、 $\theta = \theta_{i,n+1}$、 $\varphi = \varphi_{i,n+1}$ の $t_{\mathrm{i}}$ を使って、式 {eq}`dc-germ` により決める。一度発芽した質点は、そのあとも発芽しているものとし、 $\tau_{i,n}$ は求めない。

$$
G_{i,n+1} = \left( \tau_{i,n+1} > t_{\mathrm{i}} \right)
$$ (dc-germ)

## 5 成長期

腐朽の反応速度定数 $k_{\mathrm{m}}$ は、式 {eq}`dc-km` により計算される。式の値が負になる温度では、0 とする。

$$
k_{\mathrm{m}} = \max \left( \left( 2.77 - 3.23 \thinspace \theta + 0.865 \thinspace \theta^2 - 0.0189 \thinspace \theta^3 \right) \times 10^{-10}, \thinspace 0 \right)
$$ (dc-km)

質点 $i$ で、次の 2 つをともに満たすとき、腐朽が進む。

- 質点 $i$ が層の表面の質点で発芽している（ $G_{i,n+1}$ が真）か、隣の質点のステップ $n$ の質量減少率が 0 より大きい（ $L_{i-1,n} > 0$ または $L_{i+1,n} > 0$）。
- 温度が $\theta_{\mathrm{g,min}} < \theta_{i,n+1} \le \theta_{\mathrm{g,max}}$ で、相対湿度が $\varphi_{i,n+1} \ge \varphi_{\mathrm{g}}$ である。

腐朽が進むとき、ステップ $n + 1$ の質量減少率 $L_{i,n+1}$ は、 $\theta = \theta_{i,n+1}$ の $k_{\mathrm{m}}$ を使って、式 {eq}`dc-l` により計算される。腐朽が進まないときは、 $L_{i,n+1} = L_{i,n}$ とする。

$$
L_{i,n+1} = \min \left( L_{i,n} + k_{\mathrm{m}} \thinspace \Delta t, \thinspace L_{\max} \right)
$$ (dc-l)

## 6 腐朽で生まれる水

ステップ $n + 1$ の水分の収支に加える、質点 $i$ の腐朽で生まれる水 $j_{\mathrm{gen},i,n+1}$ は、ステップ $n - 1$ から $n$ までの質量減少率の増えから、式 {eq}`dc-gen` により計算される。ステップ $n + 1$ の未知数（温度と水分化学ポテンシャル）によらない。最初のステップ（ $n = 0$）では、0 とする。木材腐朽を求めない質点では、0 とする。

$$
j_{\mathrm{gen},i,n+1} = h \thinspace \rho_{\mathrm{w}} \frac{L_{i,n} - L_{i,n-1}}{\Delta t} \thinspace w_i
$$ (dc-gen)

## 7 根拠

発芽期と成長期の式、腐朽が表面から始まり内部へ広がるとする理由、時間刻みごとに判定する理由、腐朽で生まれる水を 1 つ前のステップの値で求める理由、反応速度定数が負になる温度で 0 とする理由は、「{doc}`木材腐朽のモデル <04_08_basis_decay>`」に示す。
