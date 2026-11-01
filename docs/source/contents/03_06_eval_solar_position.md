# 太陽位置

## 1 はじめに

### 1.1 概要

本節では、各時刻の太陽の位置（太陽高度と太陽方位角）の計算方法を定める。太陽の位置は、外表面に当たる日射の計算（「{doc}`外表面の日射・夜間放射 <03_09_eval_surface_solar>`」）に使う。

### 1.2 入力

:::{container} symbols
$\phi_{\mathrm{lat}}$
: ：緯度。北緯を正とする（deg）

$\lambda_{\mathrm{lon}}$
: ：経度。東経を正とする（deg）

$t_n$
: ：ステップ $n$ の時刻。1 月 1 日 0 時からの経過時間（s）
:::

### 1.3 出力

:::{container} symbols
$\sin h_{\mathrm{sun},n}$, $\cos h_{\mathrm{sun},n}$
: ：ステップ $n$ の太陽高度の正弦と余弦（-）

$\sin a_{\mathrm{sun},n}$, $\cos a_{\mathrm{sun},n}$
: ：ステップ $n$ の太陽方位角の正弦と余弦（-）
:::

太陽高度と太陽方位角は、角度そのものではなく、正弦と余弦の組で求める。外表面に当たる日射の計算では、正弦と余弦だけを使うためである。

## 2 記号及び添え字

### 2.1 記号

:::{container} symbols
$a_{\mathrm{sun},n}$
: ：ステップ $n$ の太陽方位角。南を 0、西回りを正とする（rad）

$d_n$
: ：ステップ $n$ の通算日。1 月 1 日を 1 とする（d）

$e_n$
: ：ステップ $n$ の均時差（h）

$h_{\mathrm{sun},n}$
: ：ステップ $n$ の太陽高度。地平線より下では負（rad）

$t_n$
: ：ステップ $n$ の時刻。1 月 1 日 0 時からの経過時間（s）

$t_{\mathrm{m},n}$
: ：ステップ $n$ の、その日の標準時（日本の標準時）（h）

$\delta_n$
: ：ステップ $n$ の赤緯（rad）

$\lambda_{\mathrm{lon}}$
: ：経度。東経を正とする（deg）

$\tau_n$
: ：ステップ $n$ の時角（rad）

$\phi_{\mathrm{lat}}$
: ：緯度。北緯を正とする（deg）

$\omega_n$
: ：ステップ $n$ の日角（rad）
:::

### 2.2 添え字

:::{container} symbols
$n$
: ：時刻のステップ
:::

## 3 通算日と標準時

1 年は 365 日とし、2 年目以降も同じ暦を繰り返す。ステップ $n$ の通算日 $d_n$ と、その日の標準時 $t_{\mathrm{m},n}$ は、式 {eq}`sp-day` と式 {eq}`sp-time` により計算される。

$$
d_n = \left\lfloor \frac{t_n \bmod (365 \times 86400)}{86400} \right\rfloor + 1
$$ (sp-day)

$$
t_{\mathrm{m},n} = \frac{t_n \bmod 86400}{3600}
$$ (sp-time)

ここで、 $\lfloor x \rfloor$ は $x$ を超えない最大の整数、 $x \bmod y$ は $x$ を $y$ で割った余りを表す。

## 4 赤緯と均時差

ステップ $n$ の日角 $\omega_n$ は、式 {eq}`sp-omega` により計算される。

$$
\omega_n = \frac{2 \pi d_n}{366}
$$ (sp-omega)

ステップ $n$ の赤緯 $\delta_n$ は、式 {eq}`sp-delta` により計算される。

$$
\begin{aligned}
\delta_n = \ & 0.006322 - 0.405748 \cos(\omega_n + 0.153231) \\
& - 0.005880 \cos(2 \omega_n + 0.207099) - 0.003233 \cos(3 \omega_n + 0.620129)
\end{aligned}
$$ (sp-delta)

ステップ $n$ の均時差 $e_n$ は、式 {eq}`sp-eot` により計算される。

$$
\begin{aligned}
e_n = \ & -0.000279 + 0.122772 \cos(\omega_n + 1.498311) \\
& - 0.165458 \cos(2 \omega_n - 1.261546) - 0.005354 \cos(3 \omega_n - 1.157100)
\end{aligned}
$$ (sp-eot)

## 5 時角

ステップ $n$ の時角 $\tau_n$ は、式 {eq}`sp-tau` により計算される。経度 135 deg（日本の標準時の子午線）との経度の差の分だけ、太陽の南中の時刻がずれることを考える。

$$
\tau_n = \frac{\pi}{12} \left\lbrace \left( t_{\mathrm{m},n} + e_n - 12 \right) + \frac{\lambda_{\mathrm{lon}} - 135}{15} \right\rbrace
$$ (sp-tau)

## 6 太陽高度と太陽方位角

ステップ $n$ の太陽高度の正弦 $\sin h_{\mathrm{sun},n}$ と余弦 $\cos h_{\mathrm{sun},n}$ は、式 {eq}`sp-sinh` と式 {eq}`sp-cosh` により計算される。

$$
\sin h_{\mathrm{sun},n} = \sin \phi_{\mathrm{lat}} \sin \delta_n + \cos \phi_{\mathrm{lat}} \cos \delta_n \cos \tau_n
$$ (sp-sinh)

$$
\cos h_{\mathrm{sun},n} = \sqrt{1 - \left( \sin h_{\mathrm{sun},n} \right)^2}
$$ (sp-cosh)

太陽高度は $-\pi/2$ 以上 $\pi/2$ 以下とするので、余弦は 0 以上である。

ステップ $n$ の太陽方位角の正弦 $\sin a_{\mathrm{sun},n}$ と余弦 $\cos a_{\mathrm{sun},n}$ は、式 {eq}`sp-sina` と式 {eq}`sp-cosa` により計算される。

$$
\sin a_{\mathrm{sun},n} = \frac{\cos \delta_n \sin \tau_n}{\cos h_{\mathrm{sun},n}}
$$ (sp-sina)

$$
\cos a_{\mathrm{sun},n} = \frac{\cos \delta_n \sin \phi_{\mathrm{lat}} \cos \tau_n - \sin \delta_n \cos \phi_{\mathrm{lat}}}{\cos h_{\mathrm{sun},n}}
$$ (sp-cosa)

なお、太陽が天頂にある場合（ $\cos h_{\mathrm{sun},n} = 0$ ）は、太陽方位角は定義されない。日本の緯度では、太陽が天頂に来ることはない。

式 {eq}`sp-sinh`〜{eq}`sp-cosa` の三角関数に入れる緯度は、rad に換算した値とする。
