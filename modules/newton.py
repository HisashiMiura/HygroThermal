"""汎用のニュートン・ラフソン法（設計の文書 2.8 節）。

物理を知らない部品。未知数の並び x_p と、残差の並び F_p（p は未知数・式の番号）を扱う。
ヤコビ行列 J_pq = ∂F_p / ∂x_q は帯の形とし、帯の幅をタプル band = (l, u)（対角の下に l 本、上に u 本）で表す。
帯の幅は解く問題の形で決まるので、物理を知っている呼ぶ側が決めて渡す。
scipy.linalg.solve_banded の詰めた形 ab（J の p 行 q 列は ab[u + p - q, q]）で持つ。

計算方法の文書の評価方法との対応：「1 ステップの解き方」の 4.3〜4.5 節。コメントの「式 (k)」はその節の式の番号
（docs/source/contents/03_12_eval_solve.md）。
"""

from collections.abc import Callable

import numpy as np
from scipy.linalg import LinAlgError, solve_banded

# 残差の関数：未知数の並び x_ps → 残差の並び F_ps
Residual = Callable[[np.ndarray], np.ndarray]

# 修正の幅を半分にする回数の上限（1/2^10 ≈ 1/1000 まで小さくする）
MAX_HALVING = 10

# 残差を計算するときに、未知数が計算できない範囲に入ったときに出るエラー
# （0 で割る、exp があふれる、負の数の対数など。数値の誤り ArithmeticError と、値の誤り ValueError。
# 負の数の小数乗は複素数になり、それを math の関数に渡すと型の誤り TypeError になる）
_RESIDUAL_ERRORS = (ArithmeticError, ValueError, TypeError)


class NewtonNotConverged(RuntimeError):
    """ニュートン法が、決めた回数のうちに収束しなかった（またはヤコビ行列が解けなかった、残差が計算できなかった）。"""


def _try_residual(f: Residual, x_ps: np.ndarray) -> np.ndarray | None:
    """残差を計算する。計算できない（エラーが出る、無限大・NaN になる）ときは None を返す。

    Args:
        f: 残差の関数
        x_ps: 未知数の並び, [p]

    Returns:
        残差の並び F_ps, [p]（計算できなければ None）
    """

    try:
        f_ps = np.asarray(f(x_ps), dtype=float)
    except _RESIDUAL_ERRORS:
        return None

    # 無限大・NaN が入っていたら、計算できなかったとみなす
    if not np.all(np.isfinite(f_ps)):
        return None

    return f_ps


def _banded_jacobian(
        f: Residual,
        x_ps: np.ndarray,
        f_ps: np.ndarray,
        band: tuple[int, int],
        x_max_ps: np.ndarray | None = None
    ) -> np.ndarray:
    """数値微分で、帯のヤコビ行列を詰めた形で作る。

    未知数 q を h_q だけずらしたときの残差の変化から、J の q 列を求める：J の q 列 ≈ (F(x + h_q e_q) - F(x)) / h_q。
    帯の幅 w = l + u + 1 ごとに離れた未知数（q, q + w, q + 2w, …）は、影響する式が重ならないので同時にずらす。
    そのため、残差を計算する回数は w 回（未知数の数によらない）。
    ずらすと上限を超える未知数は、反対の向きにずらす（上限の外の、残差が平らなところで傾きを測らないため）。

    Args:
        f: 残差の関数
        x_ps: 未知数の並び（ここでのヤコビ行列を求める）, [p]
        f_ps: x_ps での残差の並び, [p]（呼ぶ側で計算済みのものを使い回す）
        band: ヤコビ行列の帯の幅 (対角の下の本数 l, 上の本数 u)
        x_max_ps: 未知数の上限, [p]（省略すると上限なし）

    Returns:
        帯のヤコビ行列の詰めた形 ab（形は (l + u + 1, 未知数の数)）

    Raises:
        NewtonNotConverged: 未知数をずらした残差が計算できなかったとき

    Notes:
        計算方法の文書の「1 ステップの解き方」の式 (6)・(7) に対応する。
    """

    x_ps = np.asarray(x_ps, dtype=float)
    n = len(x_ps)

    # 帯の幅：対角の下の本数 l、上の本数 u
    l, u = band

    # 帯の幅
    w = l + u + 1

    # ずらす量 h_q, [q]（値の大きさに比例させる。値が 0 に近いときは 1 を基準にする）
    # ずらすと上限を超える未知数は、反対の向きにずらす
    # 式 (7)
    h_ps = np.sqrt(np.finfo(float).eps) * np.maximum(np.abs(x_ps), 1.0)
    if x_max_ps is not None:
        h_ps = np.where(x_ps + h_ps > x_max_ps, -h_ps, h_ps)

    ab = np.zeros((w, n))

    # 帯の幅 w ごとに離れた未知数の組を、1 組ずつまとめてずらす
    for start in range(w):

        # この組の未知数の番号（start, start + w, start + 2w, …）
        q_qs = np.arange(start, n, w)

        # この組の未知数だけをずらした残差の変化
        x_shifted_ps = x_ps.copy()
        x_shifted_ps[q_qs] += h_ps[q_qs]
        f_shifted_ps = _try_residual(f, x_shifted_ps)
        if f_shifted_ps is None:
            raise NewtonNotConverged('ヤコビ行列を作るときに、残差が計算できませんでした（未知数が計算できない範囲に入った）。')
        df_ps = f_shifted_ps - f_ps

        # 未知数 q が影響する式は p = q - u 〜 q + l。その変化を h_q で割ったものが J の q 列
        # 式 (6)
        for q in q_qs:
            for p in range(max(0, q - u), min(n, q + l + 1)):
                ab[u + p - q, q] = df_ps[p] / h_ps[q]

    return ab


def newton_solve(
        f: Residual,
        x0_ps: np.ndarray,
        band: tuple[int, int],
        tol: float,
        max_iter: int = 20,
        min_iter: int = 0,
        x_max_ps: np.ndarray | None = None
    ) -> tuple[np.ndarray, int]:
    """ニュートン・ラフソン法で、残差 F(x) = 0 となる x を求める。

    m 回目：残差 F と帯のヤコビ行列 J を作り、min_iter 回以上修正していて、残差のいちばん大きいものが tol 以下なら終わる。
    そうでなければ、J Δx = -F を帯行列の直接法で解いて、x を x + Δx に更新する。
    x + Δx で残差が計算できないか、残差の二乗の和が前より小さくならない（いちばん大きいものも許容値以下にならない）ときは、
    Δx を半分にして試す
    （直線探索。最大 MAX_HALVING 回）。答えの近くでは Δx がそのまま通るので、ふつうのニュートン法と同じになる。
    上限 x_max_ps を渡すと、修正した x が上限を超えるときは上限に置き換える（上限の外で残差が平らになり、
    ヤコビ行列が特異に近くなって、次の修正が途方もなく大きくなるのを防ぐ）。

    Args:
        f: 残差の関数
        x0_ps: 最初の仮の値, [p]
        band: ヤコビ行列の帯の幅 (対角の下の本数 l, 上の本数 u)。熱だけなら (1, 1)
        tol: 残差の許容値（残差と同じ単位。熱なら W/m2）
        max_iter: 修正の回数の上限
        min_iter: 修正の回数の下限（最初の仮の値の残差が tol 以下でも、この回数は修正する）
        x_max_ps: 未知数の上限, [p]（上限のない未知数は inf。省略すると上限なし）

    Returns:
        答え x_ps, [p]、修正した回数

    Raises:
        NewtonNotConverged: max_iter 回修正しても残差が tol 以下にならなかったとき、ヤコビ行列が解けなかったとき、
            修正の幅を小さくしても残差が小さくならなかったとき、または最初の仮の値で残差が計算できなかったとき

    Notes:
        計算方法の文書の「1 ステップの解き方」の 4.3〜4.5 節（式 (6)〜(10)）に対応する。
    """

    x_ps = np.array(x0_ps, dtype=float)

    # 未知数の上限, [p]（省略したら、どの未知数も上限なし）
    x_max_ps = np.full(len(x_ps), np.inf) if x_max_ps is None else np.asarray(x_max_ps, dtype=float)

    # 最初の仮の値の残差（計算できなければ、始められないので止める）
    f_ps = _try_residual(f, x_ps)
    if f_ps is None:
        raise NewtonNotConverged('最初の仮の値で、残差が計算できませんでした。')

    for m in range(max_iter + 1):

        # 下限の回数まで修正していて、残差のいちばん大きいものが許容値以下なら、答え（ここまでに m 回修正した）
        # 4.5 節
        if m >= min_iter and np.max(np.abs(f_ps)) <= tol:
            return x_ps, m

        # 修正の回数の上限に達していたら、ここで止める（下で例外を出す）
        if m == max_iter:
            break

        # 帯のヤコビ行列
        ab = _banded_jacobian(f, x_ps, f_ps, band, x_max_ps)

        # J Δx = -F を解く。ヤコビ行列が解けない（特異な）ときは、収束しなかったとして止める
        # 式 (8)
        try:
            dx_ps = solve_banded(band, ab, -f_ps)
        except LinAlgError:
            raise NewtonNotConverged(f'ニュートン法の {m + 1} 回目で、ヤコビ行列が解けませんでした（特異）。') from None

        # 修正の幅を 1, 1/2, 1/4, … と小さくしながら、残差が計算できて、二乗の和が前より小さくなる
        # （または、いちばん大きいものが許容値以下になる）x を探す
        # 式 (9)・(10)
        for s in range(MAX_HALVING + 1):
            x_new_ps = np.minimum(x_ps + dx_ps / 2 ** s, x_max_ps)
            f_new_ps = _try_residual(f, x_new_ps)
            if f_new_ps is not None and (
                    np.sum(f_new_ps ** 2) < np.sum(f_ps ** 2) or np.max(np.abs(f_new_ps)) <= tol):
                break
        else:
            raise NewtonNotConverged(
                f'ニュートン法の {m + 1} 回目で、修正の幅を 1/{2 ** MAX_HALVING} まで小さくしても、残差が小さくなりませんでした'
                f'（残差の最大 {np.max(np.abs(f_ps)):.3g}、許容値 {tol:.3g}）。')

        # x を更新する
        x_ps, f_ps = x_new_ps, f_new_ps

    raise NewtonNotConverged(f'ニュートン法が {max_iter} 回の修正で収束しませんでした（残差の最大 {np.max(np.abs(f_ps)):.3g}、許容値 {tol:.3g}）。')
