"""時間を 1 ステップ進める：熱と水分を一緒に解く（設計の文書 2.8 節）。

熱と水分の収支の残差（balance.py）を、汎用のニュートン法（newton.py）で 0 にして、次の時刻の状態を作る。

計算方法の文書の評価方法との対応：「1 ステップの解き方」。コメントの「式 (k)」は、とくに断りがなければその節の式の番号
（docs/source/contents/03_12_eval_solve.md）。
"""

import numpy as np

from modules.balance import heat_residual_is, moisture_residual_is, ventilation_rate_is
from modules.boundary import IndoorBoundary, OutdoorBoundary
from modules.newton import NewtonNotConverged, newton_solve
from modules.rain import get_rain_on_surface
from modules.surface_radiation import get_solar_heat_gain
from modules.air_water_properties import L_V
from modules.wall_geometry import FilmNode, WallGeometry
from modules.wall_state import WallState

# ヤコビ行列の帯の幅（下 3 本・上 3 本）：未知数を (T_0, μ_0, T_1, μ_1, …) と並べると、
# 質点 i の熱と水分の残差は、質点 i - 1・i・i + 1 の T と μ（左右に 3 つずつ離れた列まで）にしか関係しない
HEAT_MOISTURE_BAND = (3, 3)

# 熱の残差の許容値 ε_h, W/m2（温度にすると、およそ 10^-7 K の精度）
HEAT_RESIDUAL_TOL = 1.0e-6

# 水分の残差の許容値 ε_m, kg/(m2 s)（潜熱にすると熱の許容値と同じ大きさ）
# 式 (3)
MOISTURE_RESIDUAL_TOL = HEAT_RESIDUAL_TOL / L_V

# 収束しなかったとき、前の時刻の相対湿度がこの値以上の質点があれば、結露したとみなす, %
CONDENSATION_RH = 99.9


class Condensation(NewtonNotConverged):
    """結露した：水分が入りすぎて、相対湿度 100% 以下（μ ≤ 0）では収支が合わない。

    結露した水（水膜・結露水）を持つモデルは段階 7 で作る。それまでは、水を消さずにここで止める。
    """


def step(
        geometry: WallGeometry,
        state_n: WallState,
        dt: float,
        outdoor_n: OutdoorBoundary,
        outdoor_npls: OutdoorBoundary,
        indoor_npls: IndoorBoundary
    ) -> tuple[WallState, int]:
    """熱と水分を一緒に 1 ステップ進める。

    未知数を (T_0, μ_0, T_1, μ_1, …) と並べ、熱と水分の残差をそれぞれの許容値で割って（許容値の何倍かにして）、
    交互に並べた残差が 0 になるように、ニュートン法で求める（許容値 1）。最初の仮の値は前の時刻の状態。
    最初の仮の値の残差が許容値より小さくても、1 回は修正する（平衡の近くで、1 ステップの変化が小さいと進まなくなるのを防ぐ）。
    収束せず、前の時刻に相対湿度が CONDENSATION_RH 以上の質点があれば、結露したとして Condensation を出す。

    Args:
        geometry: 壁の形
        state_n: 前の時刻（ステップ n）の状態
        dt: 時間刻み, s
        outdoor_n: ステップ n の外気の境界条件（通気層の換気量だけに使う）
        outdoor_npls: ステップ n + 1 の外気の境界条件
        indoor_npls: ステップ n + 1 の室内の境界条件

    Returns:
        次の時刻の状態、ニュートン法の修正の回数

    Raises:
        Condensation: 結露して、相対湿度 100% 以下では収支が合わないとき
        NewtonNotConverged: それ以外の理由で収束しなかったとき

    Notes:
        計算方法の文書の「1 ステップの解き方」に対応する（ニュートン法の中身は newton.newton_solve）。
    """

    # ステップ n + 1 の、外表面に外から直接入る熱（日射の吸収 - 夜間放射）, W/m2（質点の状態によらないので、1 ステップの中では一定）
    # 評価方法の「外表面の日射・夜間放射」のステップ n に、ここではステップ n + 1 の外気の値を入れる（後退差分）
    q_sol_npls = get_solar_heat_gain(geometry.surface, outdoor_npls.solar)

    # ステップ n + 1 の、外表面に当たる雨水 r_bv,n+1, kg/(m2 s)（質点の状態によらないので、1 ステップの中では一定）
    # 「雨水の浸入」の式 (2)
    r_bv_npls = get_rain_on_surface(geometry.surface, outdoor_npls)

    # ステップ n の、通気層の壁 1 m2 あたりの換気量 q_v,i,n, m3/(m2 s), [i]（ステップ n の質点と外気の温度で求め、1 ステップの中では一定にする。
    # 換気量は √|ρ_a,ex - ρ_a,i| に比例し、外気と同じ温度の所で傾きが無限大になるので、ニュートン法の未知数に含めない）
    # 「通気層の換気」の式 (2)
    q_v_is_n = ventilation_rate_is(geometry, state_n.t_is, outdoor_n)

    # 未知数の並び x_ps を、温度と水分化学ポテンシャルに分けて、許容値で割った残差を交互に並べて返す関数（この 1 ステップ専用）
    def residual(x_ps: np.ndarray) -> np.ndarray:

        # ステップ n + 1 の温度は 0 番目から 1 つおき、水分化学ポテンシャルは 1 番目から 1 つおき
        # 式 (1)
        t_is_npls, mu_is_npls = x_ps[0::2], x_ps[1::2]

        # 熱と水分の残差を、それぞれの許容値で割る
        # 式 (2)
        res_h_is = heat_residual_is(
            geometry, t_is_npls, mu_is_npls, state_n.t_is, dt, outdoor_npls, indoor_npls, q_sol_npls, q_v_is_n, r_bv_npls) / HEAT_RESIDUAL_TOL
        res_m_is = moisture_residual_is(
            geometry, t_is_npls, mu_is_npls, state_n.t_is, state_n.mu_is, dt, outdoor_npls, indoor_npls, q_v_is_n, r_bv_npls) / MOISTURE_RESIDUAL_TOL

        # 熱と水分の残差を交互に並べる
        # 式 (2)
        f_ps = np.empty(len(x_ps))
        f_ps[0::2], f_ps[1::2] = res_h_is, res_m_is
        return f_ps

    # 前の時刻の状態を、(T_0, μ_0, T_1, μ_1, …) の順に並べて、最初の仮の値にする
    # 式 (4)
    x0_ps = np.empty(2 * len(state_n))
    x0_ps[0::2], x0_ps[1::2] = state_n.t_is, state_n.mu_is

    # 未知数の上限：温度は上限なし、水分化学ポテンシャルは 0（相対湿度 100%）、水膜の状態量は上限なし
    # （μ > 0 では相対湿度が 100% で頭打ちになり、残差が μ によらなくなるので、ニュートン法の修正がそこへ行き過ぎないようにする。
    #   水膜の状態量 s ≥ 0 は水膜があることを表し、残差は水膜の水分量で変わるので、上限を設けない）
    # 式 (5)
    is_film_is = np.array([isinstance(node, FilmNode) for node in geometry.nodes])
    x_max_ps = np.empty(len(x0_ps))
    x_max_ps[0::2], x_max_ps[1::2] = np.inf, np.where(is_film_is, np.inf, 0.0)

    # 残差が 0（許容値の 1 倍以下）になる未知数を探す（1 回は必ず修正する）
    try:
        x_ps, m = newton_solve(residual, x0_ps, band=HEAT_MOISTURE_BAND, tol=1.0, min_iter=1, x_max_ps=x_max_ps)
    except NewtonNotConverged as e:

        # 収束しないとき（5 節）

        # 前の時刻に、相対湿度が 100% に近い質点の番号（水膜の質点は、水がたまっても収支が合うので除く）
        rh_is_n = state_n.rh_is
        i_wet_is = np.flatnonzero((rh_is_n >= CONDENSATION_RH) & ~is_film_is)

        # 100% に近い質点がなければ、結露ではない理由で収束しなかったので、そのまま例外を出す
        if len(i_wet_is) == 0:
            raise

        # 100% に近い質点があれば、結露したとして、質点の番号・層の番号・相対湿度を伝える
        where = '、'.join(
            f'質点 {i}（層 {geometry.nodes[i].layer_index}、{rh_is_n[i]:.3f}%）' for i in i_wet_is
        )
        raise Condensation(f'結露しました（相対湿度 100% 以下では水分の収支が合いません）：{where}。') from e

    # 答えを、温度と水分化学ポテンシャルに分ける
    # 式 (11)
    return WallState(t_is=x_ps[0::2], mu_is=x_ps[1::2]), m
