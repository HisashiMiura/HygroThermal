"""収支の計算：つながりごとの流れと、質点ごとの残差（設計の文書 2.8 節）。

物理（流れと収支）はここに書く。解く部分（ニュートン法）は物理を知らない汎用の部品にする。

計算方法の文書の評価方法との対応：流れ（*_flow_ks）は「熱・水分の流れ」、換気量（ventilation_rate_is）は「通気層の換気」、
水分量・熱容量と残差（water_amount_is・heat_capacity_is・*_residual_is）は「熱・水分の収支」。
流れの関数のコメントの「式 (k)」は「熱・水分の流れ」（docs/source/contents/03_05_eval_flows.md）の式の番号。
"""

import numpy as np

from modules.boundary import IndoorBoundary, OutdoorBoundary
from modules.state import State
from modules.air_water_properties import C_A, C_W, L_V, R_V, get_rho_a
from modules.wall_geometry import AirNode, FilmNode, MaterialNode, WallGeometry


def heat_flow_ks(
        geometry: WallGeometry,
        t_is_n: np.ndarray,
        outdoor_n: OutdoorBoundary,
        indoor_n: IndoorBoundary
    ) -> np.ndarray:
    """ステップ n の、つながりごとの熱の流れ q_k,n, W/m2, [k]（外気側から室内側へ流れる向きを正）。

    つながり k は質点 k - 1 と質点 k の間（0 番目は外気と質点 0、最後は最後の質点と室内）。
    残差（熱・水分の収支）からは、ステップ n + 1 の値を渡す。

    Args:
        geometry: 壁の形
        t_is_n: ステップ n の各質点の温度 T_i,n, K, [i]
        outdoor_n: ステップ n の外気の境界条件（外気温 T_ex,n を使う）
        indoor_n: ステップ n の室内の境界条件（室温 T_in,n を使う）

    Returns:
        ステップ n の、つながりごとの熱の流れ q_k,n, W/m2, [k]

    Notes:
        計算方法の文書の「熱・水分の流れ」の式 (9) に対応する。
    """

    # 各つながりの熱抵抗 R_h,k, (m2 K)/W, [k]
    # 式 (6)（Link.r_h）
    r_h_ks = np.array([link.r_h for link in geometry.links])

    # つながり k の左の温度 T_k-1,n と右の温度 T_k,n, K, [k]
    # T_-1,n = T_ex,n、T_N,n = T_in,n とする
    t_left_ks_n = np.concatenate(([outdoor_n.state.t], t_is_n))
    t_right_ks_n = np.concatenate((t_is_n, [indoor_n.state.t]))

    # 熱の流れ q_k,n, W/m2, [k]
    # 式 (9)
    q_ks_n = (t_left_ks_n - t_right_ks_n) / r_h_ks

    return q_ks_n


def _mu_is(geometry: WallGeometry, x_is_n: np.ndarray) -> np.ndarray:
    """ステップ n の各質点の水分化学ポテンシャル μ_i,n, J/kg, [i]。

    状態の水分化学ポテンシャルの欄は、水膜の質点と材料の質点では状態量 s_i,n なので、水分化学ポテンシャルに直す
    （空気の質点はそのまま。解き方で 0 以下に抑えている）。

    Args:
        geometry: 壁の形
        x_is_n: ステップ n の状態の水分化学ポテンシャルの欄, [i]（水膜の質点と材料の質点では s_i,n）

    Returns:
        ステップ n の各質点の水分化学ポテンシャル μ_i,n, J/kg, [i]

    Notes:
        計算方法の文書の「熱・水分の収支」の式 (10)（水膜の質点）・(16)（材料の質点）に対応する。
    """
    return np.array([
        node.mu(x) if isinstance(node, (FilmNode, MaterialNode)) else x
        for node, x in zip(geometry.nodes, x_is_n)
    ])


def _state_left_right_ks(
        t_is_n: np.ndarray,
        mu_is_n: np.ndarray,
        outdoor_n: OutdoorBoundary,
        indoor_n: IndoorBoundary
    ) -> tuple[list[State], list[State]]:
    """ステップ n の、つながり k の左の状態と右の状態の並び, [k]（左の 0 番目は外気、右の最後は室内）。

    Args:
        t_is_n: ステップ n の各質点の温度 T_i,n, K, [i]
        mu_is_n: ステップ n の各質点の水分化学ポテンシャル μ_i,n, J/kg, [i]
        outdoor_n: ステップ n の外気の境界条件
        indoor_n: ステップ n の室内の境界条件

    Returns:
        つながり k の左の状態（質点 k - 1。0 番目は外気）の並び、右の状態（質点 k。最後は室内）の並び
    """

    # 各質点の 1 点の状態, [i]
    state_is_n = [State(t=t_i_n, mu=mu_i_n) for t_i_n, mu_i_n in zip(t_is_n, mu_is_n)]

    # つながり k の左は質点 k - 1（0 番目は外気）、右は質点 k（最後は室内）
    return [outdoor_n.state] + state_is_n, state_is_n + [indoor_n.state]


def vapour_flow_ks(
        geometry: WallGeometry,
        t_is_n: np.ndarray,
        mu_is_n: np.ndarray,
        outdoor_n: OutdoorBoundary,
        indoor_n: IndoorBoundary
    ) -> np.ndarray:
    """ステップ n の、つながりごとの水蒸気の流れ j_v,k,n, kg/(m2 s), [k]（外気側から室内側へ流れる向きを正）。

    水蒸気圧は、質点の温度と水分化学ポテンシャルから求める（偏微分に分けない）。
    透湿抵抗は、左右それぞれの側の状態（相対湿度）で計算する。
    残差（熱・水分の収支）からは、ステップ n + 1 の値を渡す。

    Args:
        geometry: 壁の形
        t_is_n: ステップ n の各質点の温度 T_i,n, K, [i]
        mu_is_n: ステップ n の各質点の水分化学ポテンシャル μ_i,n, J/kg, [i]
        outdoor_n: ステップ n の外気の境界条件
        indoor_n: ステップ n の室内の境界条件

    Returns:
        ステップ n の、つながりごとの水蒸気の流れ j_v,k,n, kg/(m2 s), [k]

    Notes:
        計算方法の文書の「熱・水分の流れ」の式 (10) に対応する。
    """

    # つながり k の左の状態と右の状態, [k]（左の 0 番目は外気、右の最後は室内。水膜の質点は水分化学ポテンシャルに直す）
    state_left_ks_n, state_right_ks_n = _state_left_right_ks(t_is_n, _mu_is(geometry, mu_is_n), outdoor_n, indoor_n)

    # つながり k の左の水蒸気圧 p_v,k-1,n と右の水蒸気圧 p_v,k,n, Pa, [k]
    # 「湿り空気と水の性質」の式 (4)
    p_v_left_ks_n = np.array([s.p_v for s in state_left_ks_n])
    p_v_right_ks_n = np.array([s.p_v for s in state_right_ks_n])

    # 各つながりの透湿抵抗 R_m,k,n, (m2 s Pa)/kg, [k]
    # 式 (7)（Link.r_m。左右それぞれの側の状態で計算する）
    r_m_ks_n = np.array([
        link.r_m(s_left, s_right)
        for link, s_left, s_right in zip(geometry.links, state_left_ks_n, state_right_ks_n)
    ])

    # 水蒸気の流れ j_v,k,n, kg/(m2 s), [k]
    # 式 (10)
    j_v_ks_n = (p_v_left_ks_n - p_v_right_ks_n) / r_m_ks_n

    return j_v_ks_n


def liquid_flow_ks(
        geometry: WallGeometry,
        t_is_n: np.ndarray,
        mu_is_n: np.ndarray,
        outdoor_n: OutdoorBoundary,
        indoor_n: IndoorBoundary
    ) -> np.ndarray:
    """ステップ n の、つながりごとの液水の流れ j_l,k,n, kg/(m2 s), [k]（外気側から室内側へ流れる向きを正）。

    境目（表面・層と層・材料と通気層）と、液相水分伝導率が 0 の所では、液水の抵抗が無限大なので流れは 0。
    残差（熱・水分の収支）からは、ステップ n + 1 の値を渡す。

    Args:
        geometry: 壁の形
        t_is_n: ステップ n の各質点の温度 T_i,n, K, [i]
        mu_is_n: ステップ n の各質点の水分化学ポテンシャル μ_i,n, J/kg, [i]
        outdoor_n: ステップ n の外気の境界条件
        indoor_n: ステップ n の室内の境界条件

    Returns:
        ステップ n の、つながりごとの液水の流れ j_l,k,n, kg/(m2 s), [k]

    Notes:
        計算方法の文書の「熱・水分の流れ」の式 (11) に対応する。
    """

    # つながり k の左の状態と右の状態, [k]（水膜の質点は水分化学ポテンシャルに直す）
    state_left_ks_n, state_right_ks_n = _state_left_right_ks(t_is_n, _mu_is(geometry, mu_is_n), outdoor_n, indoor_n)

    # つながり k の左の水分化学ポテンシャル μ_k-1,n と右の μ_k,n, J/kg, [k]
    mu_left_ks_n = np.array([s.mu for s in state_left_ks_n])
    mu_right_ks_n = np.array([s.mu for s in state_right_ks_n])

    # 各つながりの液水の抵抗 R_l,k,n, (m2 s (J/kg))/kg, [k]
    # 式 (8)（Link.r_l。境目と、液相水分伝導率が 0 の所では無限大）
    r_l_ks_n = np.array([
        link.r_l(s_left, s_right)
        for link, s_left, s_right in zip(geometry.links, state_left_ks_n, state_right_ks_n)
    ])

    # 液水の流れ j_l,k,n, kg/(m2 s), [k]
    # 式 (11)。抵抗が無限大のつながりは、差 / 無限大 = 0
    j_l_ks_n = (mu_left_ks_n - mu_right_ks_n) / r_l_ks_n

    return j_l_ks_n


def water_amount_is(geometry: WallGeometry, t_is_n: np.ndarray, mu_is_n: np.ndarray) -> np.ndarray:
    """ステップ n の、各質点の水分量 W_i,n, kg/m2, [i]。

    Args:
        geometry: 壁の形
        t_is_n: ステップ n の各質点の温度 T_i,n, K, [i]
        mu_is_n: ステップ n の各質点の水分化学ポテンシャル μ_i,n, J/kg, [i]

    Returns:
        ステップ n の、各質点の水分量 W_i,n, kg/m2, [i]

    Notes:
        計算方法の文書の「熱・水分の収支」の式 (1)・(2)・(11)（MaterialNode・AirNode・FilmNode の water_amount）を、質点ごとに並べる。
        水膜の質点では、状態の水分化学ポテンシャルの欄の水膜の状態量 s_i,n から、水膜の水分量を求める。
    """
    w_is_n = np.array([
        node.water_amount(State(t=t_i_n, mu=mu_i_n))
        for node, t_i_n, mu_i_n in zip(geometry.nodes, t_is_n, mu_is_n)
    ])
    return w_is_n


def heat_capacity_is(geometry: WallGeometry, t_is_n: np.ndarray, mu_is_n: np.ndarray) -> np.ndarray:
    """ステップ n の、各質点の熱容量 C_i,n, J/(m2 K), [i]。

    Args:
        geometry: 壁の形
        t_is_n: ステップ n の各質点の温度 T_i,n, K, [i]
        mu_is_n: ステップ n の各質点の水分化学ポテンシャル μ_i,n, J/kg, [i]

    Returns:
        ステップ n の、各質点の熱容量 C_i,n, J/(m2 K), [i]

    Notes:
        計算方法の文書の「熱・水分の収支」の式 (3)・(4)（MaterialNode・AirNode の heat_capacity）を、質点ごとに並べる。水膜の質点は 0（6 節）。
    """
    c_is_n = np.array([
        node.heat_capacity(State(t=t_i_n, mu=mu_i_n))
        for node, t_i_n, mu_i_n in zip(geometry.nodes, t_is_n, mu_is_n)
    ])
    return c_is_n


def ventilation_rate_is(geometry: WallGeometry, t_is_n: np.ndarray, outdoor_n: OutdoorBoundary) -> np.ndarray:
    """ステップ n の、各質点の壁 1 m2 あたりの換気量 q_v,i,n, m3/(m2 s), [i]（通気層の空気の質点だけ。材料の質点は 0）。

    Args:
        geometry: 壁の形
        t_is_n: ステップ n の各質点の温度 T_i,n, K, [i]
        outdoor_n: ステップ n の外気の境界条件（外気温 T_ex,n を使う）

    Returns:
        ステップ n の、各質点の壁 1 m2 あたりの換気量 q_v,i,n, m3/(m2 s), [i]

    Notes:
        計算方法の文書の「通気層の換気」の式 (2)（AirNode.ventilation_rate）を、質点ごとに並べる。
    """
    q_v_is_n = np.array([
        node.ventilation_rate(t_i_n, outdoor_n.state.t, geometry.height, geometry.len_vertical) if isinstance(node, AirNode) else 0.0
        for node, t_i_n in zip(geometry.nodes, t_is_n)
    ])
    return q_v_is_n


def rain_is(geometry: WallGeometry, r_bv_npls: float) -> np.ndarray:
    """ステップ n + 1 の、各質点に浸入する雨水 j_rain,i,n+1, kg/(m2 s), [i]（水膜の質点だけ。ほかの質点は 0）。

    Args:
        geometry: 壁の形
        r_bv_npls: ステップ n + 1 の、外表面に当たる雨水 r_bv,n+1, kg/(m2 s)（rain.get_rain_on_surface）

    Returns:
        ステップ n + 1 の、各質点に浸入する雨水 j_rain,i,n+1, kg/(m2 s), [i]

    Notes:
        計算方法の文書の「熱・水分の収支」の式 (14) に対応する。
    """

    # 浸入する雨水 j_rain,i,n+1, kg/(m2 s), [i]
    # 式 (14)（「熱・水分の収支」）
    j_rain_is_npls = np.array([node.rain_fraction * r_bv_npls if isinstance(node, FilmNode) else 0.0 for node in geometry.nodes])
    return j_rain_is_npls


def drainage_is(geometry: WallGeometry, t_is_n: np.ndarray, mu_is_n: np.ndarray) -> np.ndarray:
    """ステップ n の、各質点から流れ落ちる水 j_drain,i,n, kg/(m2 s), [i]（水膜の質点だけ。ほかの質点は 0）。

    Args:
        geometry: 壁の形
        t_is_n: ステップ n の各質点の温度 T_i,n, K, [i]
        mu_is_n: ステップ n の各質点の状態の水分化学ポテンシャルの欄, [i]（水膜の質点では水膜の状態量 s_i,n）

    Returns:
        ステップ n の、各質点から流れ落ちる水 j_drain,i,n, kg/(m2 s), [i]

    Notes:
        計算方法の文書の「熱・水分の収支」の式 (13)（FilmNode.drainage）を、質点ごとに並べる。
    """
    j_drain_is_n = np.array([
        node.drainage(State(t=t_i_n, mu=mu_i_n)) if isinstance(node, FilmNode) else 0.0
        for node, t_i_n, mu_i_n in zip(geometry.nodes, t_is_n, mu_is_n)
    ])
    return j_drain_is_n


def moisture_residual_is(
        geometry: WallGeometry,
        t_is_npls: np.ndarray,
        mu_is_npls: np.ndarray,
        t_is_n: np.ndarray,
        mu_is_n: np.ndarray,
        dt: float,
        outdoor_npls: OutdoorBoundary,
        indoor_npls: IndoorBoundary,
        q_v_is_n: np.ndarray | None = None,
        r_bv_npls: float = 0.0,
        j_gen_is_npls: np.ndarray | None = None
    ) -> np.ndarray:
    """質点ごとの水分の収支の残差 res_m,i, kg/(m2 s), [i]。

    たまる量 -（左のつながりの流れ - 右のつながりの流れ + 換気で入る水分 + 浸入する雨水 - 流れ落ちる水 + 腐朽で生まれる水）。
    正しいステップ n + 1 の温度・水分化学ポテンシャルを入れると 0 になる。

    Args:
        geometry: 壁の形
        t_is_npls: ステップ n + 1 の各質点の温度（仮の値） T_i,n+1, K, [i]
        mu_is_npls: ステップ n + 1 の各質点の水分化学ポテンシャル（仮の値） μ_i,n+1, J/kg, [i]
        t_is_n: ステップ n の各質点の温度 T_i,n, K, [i]
        mu_is_n: ステップ n の各質点の水分化学ポテンシャル μ_i,n, J/kg, [i]
        dt: 時間刻み Δt, s
        outdoor_npls: ステップ n + 1 の外気の境界条件
        indoor_npls: ステップ n + 1 の室内の境界条件
        q_v_is_n: ステップ n の、各質点の壁 1 m2 あたりの換気量 q_v,i,n, m3/(m2 s), [i]（ventilation_rate_is。省くと換気なし）
        r_bv_npls: ステップ n + 1 の、外表面に当たる雨水 r_bv,n+1, kg/(m2 s)（rain.get_rain_on_surface。省くと雨なし）
        j_gen_is_npls: ステップ n + 1 の、各質点の腐朽で生まれる水 j_gen,i,n+1, kg/(m2 s), [i]
            （decay.water_generation_is。省くと 0）

    Returns:
        質点ごとの水分の収支の残差 res_m,i, kg/(m2 s), [i]

    Notes:
        計算方法の文書の「熱・水分の収支」の式 (8) に対応する（左辺 - 右辺）。換気で入る水分は式 (5)、
        浸入する雨水は式 (14)（「熱・水分の収支」）、流れ落ちる水は式 (13)（「熱・水分の収支」）。
        腐朽で生まれる水は「木材腐朽」の式 (10)。
    """

    # 換気量を省いたときは、換気なし（どの質点も 0）
    if q_v_is_n is None:
        q_v_is_n = np.zeros(len(t_is_npls))

    # 腐朽で生まれる水を省いたときは、0
    if j_gen_is_npls is None:
        j_gen_is_npls = np.zeros(len(t_is_npls))

    # ステップ n + 1 とステップ n の水分量 W_i,n+1・W_i,n, kg/m2, [i]
    # 式 (1)・(2)（「熱・水分の収支」）
    w_is_npls = water_amount_is(geometry, t_is_npls, mu_is_npls)
    w_is_n = water_amount_is(geometry, t_is_n, mu_is_n)

    # ステップ n + 1 の、各つながりの水分の流れ（水蒸気 + 液水） j_v,k,n+1 + j_l,k,n+1, kg/(m2 s), [k]
    # 「熱・水分の流れ」の式 (10)・(11)
    j_ks_npls = (vapour_flow_ks(geometry, t_is_npls, mu_is_npls, outdoor_npls, indoor_npls)
                 + liquid_flow_ks(geometry, t_is_npls, mu_is_npls, outdoor_npls, indoor_npls))

    # 質点 i の左のつながり（つながり i）と右のつながり（つながり i + 1）の水分の流れ, kg/(m2 s), [i]
    j_left_is_npls = j_ks_npls[:-1]
    j_right_is_npls = j_ks_npls[1:]

    # ステップ n + 1 の外気と各質点の水蒸気の密度 ρ_v,ex,n+1・ρ_v,i,n+1, kg/m3
    # 「湿り空気と水の性質」の式 (5)
    rho_v_ex_npls = outdoor_npls.state.p_v / (R_V * outdoor_npls.state.t)
    rho_v_is_npls = np.array([State(t=t_i, mu=mu_i).p_v for t_i, mu_i in zip(t_is_npls, _mu_is(geometry, mu_is_npls))]) / (R_V * t_is_npls)

    # 換気で入る水分 j_vent,i,n+1, kg/(m2 s), [i]（換気量が 0 の質点は 0）
    # 式 (5)（「熱・水分の収支」）
    j_vent_is_npls = q_v_is_n * (rho_v_ex_npls - rho_v_is_npls)

    # 水分の収支の残差 res_m,i, kg/(m2 s), [i]
    # 式 (8) の左辺 - 右辺（「熱・水分の収支」）
    # 浸入する雨水 j_rain,i,n+1 と流れ落ちる水 j_drain,i,n+1, kg/(m2 s), [i]（水膜の質点だけ）
    # 式 (14)・(13)（「熱・水分の収支」）
    j_rain_is_npls = rain_is(geometry, r_bv_npls)
    j_drain_is_npls = drainage_is(geometry, t_is_npls, mu_is_npls)

    res_m_is = (w_is_npls - w_is_n) / dt - (
        j_left_is_npls - j_right_is_npls + j_vent_is_npls + j_rain_is_npls - j_drain_is_npls + j_gen_is_npls)
    return res_m_is


def heat_residual_is(
        geometry: WallGeometry,
        t_is_npls: np.ndarray,
        mu_is_npls: np.ndarray,
        t_is_n: np.ndarray,
        dt: float,
        outdoor_npls: OutdoorBoundary,
        indoor_npls: IndoorBoundary,
        q_sol_npls: float = 0.0,
        q_v_is_n: np.ndarray | None = None,
        r_bv_npls: float = 0.0
    ) -> np.ndarray:
    """質点ごとの熱の収支の残差 res_h,i, W/m2, [i]。

    たまる量 -（左のつながりの熱の流れ - 右のつながりの熱の流れ + 凝縮熱 + 外から直接入る熱 + 換気で入る熱 + 雨水が持ち込む熱）。
    正しいステップ n + 1 の温度・水分化学ポテンシャルを入れると 0 になる。

    Args:
        geometry: 壁の形
        t_is_npls: ステップ n + 1 の各質点の温度（仮の値） T_i,n+1, K, [i]
        mu_is_npls: ステップ n + 1 の各質点の水分化学ポテンシャル（仮の値） μ_i,n+1, J/kg, [i]
        t_is_n: ステップ n の各質点の温度 T_i,n, K, [i]
        dt: 時間刻み Δt, s
        outdoor_npls: ステップ n + 1 の外気の境界条件
        indoor_npls: ステップ n + 1 の室内の境界条件
        q_sol_npls: ステップ n + 1 の、外表面に外から直接入る熱 q_sol,n+1, W/m2（surface_radiation.get_solar_heat_gain）
        q_v_is_n: ステップ n の、各質点の壁 1 m2 あたりの換気量 q_v,i,n, m3/(m2 s), [i]（ventilation_rate_is。省くと換気なし）
        r_bv_npls: ステップ n + 1 の、外表面に当たる雨水 r_bv,n+1, kg/(m2 s)（rain.get_rain_on_surface。省くと雨なし）

    Returns:
        質点ごとの熱の収支の残差 res_h,i, W/m2, [i]

    Notes:
        計算方法の文書の「熱・水分の収支」の式 (9) に対応する（左辺 - 右辺）。換気で入る熱は式 (6)、凝縮熱は式 (7)、
        雨水が持ち込む熱は式 (15)（「熱・水分の収支」）。
    """

    # 換気量を省いたときは、換気なし（どの質点も 0）
    if q_v_is_n is None:
        q_v_is_n = np.zeros(len(t_is_npls))

    # ステップ n + 1 の熱容量 C_i,n+1, J/(m2 K), [i]
    # 式 (3)・(4)（「熱・水分の収支」）
    c_is_npls = heat_capacity_is(geometry, t_is_npls, mu_is_npls)

    # ステップ n + 1 の、各つながりの熱の流れ q_k,n+1 と水蒸気の流れ j_v,k,n+1, W/m2・kg/(m2 s), [k]
    # 「熱・水分の流れ」の式 (9)・(10)
    q_ks_npls = heat_flow_ks(geometry, t_is_npls, outdoor_npls, indoor_npls)
    j_v_ks_npls = vapour_flow_ks(geometry, t_is_npls, mu_is_npls, outdoor_npls, indoor_npls)

    # 質点 i の左のつながり（つながり i）と右のつながり（つながり i + 1）の熱の流れ, W/m2, [i]
    q_left_is_npls = q_ks_npls[:-1]
    q_right_is_npls = q_ks_npls[1:]

    # 凝縮熱 q_cond,i,n+1, W/m2, [i]（材料の質点と水膜の質点は δ_i = 1。通気層の空気の質点は相変化しないので 0）
    # 式 (7)（「熱・水分の収支」）
    is_material_is = np.array([not isinstance(node, AirNode) for node in geometry.nodes])
    q_cond_is_npls = np.where(is_material_is, L_V * (j_v_ks_npls[:-1] - j_v_ks_npls[1:]), 0.0)

    # 外から直接入る熱 q_sol,i,n+1, W/m2, [i]（外表面の質点 0 だけ）
    q_sol_is_npls = np.zeros(len(t_is_npls))
    q_sol_is_npls[0] = q_sol_npls

    # ステップ n + 1 の外気の密度 ρ_a,ex,n+1, kg/m3（入ってくる外気の温度で求める）
    # 「湿り空気と水の性質」の式 (6)
    rho_a_ex_npls = get_rho_a(outdoor_npls.state.t)

    # 換気で入る熱 q_vent,i,n+1, W/m2, [i]（換気量が 0 の質点は 0）
    # 式 (6)（「熱・水分の収支」）
    q_vent_is_npls = C_A * rho_a_ex_npls * q_v_is_n * (outdoor_npls.state.t - t_is_npls)

    # 雨水が持ち込む熱 q_rain,i,n+1, W/m2, [i]（雨水は外気の温度で入ってくる。浸入する雨水が 0 の質点は 0）
    # 式 (15)（「熱・水分の収支」）
    q_rain_is_npls = C_W * rain_is(geometry, r_bv_npls) * (outdoor_npls.state.t - t_is_npls)

    # 熱の収支の残差 res_h,i, W/m2, [i]
    # 式 (9) の左辺 - 右辺（「熱・水分の収支」）
    res_h_is = c_is_npls * (t_is_npls - t_is_n) / dt - (
        q_left_is_npls - q_right_is_npls + q_cond_is_npls + q_sol_is_npls + q_vent_is_npls + q_rain_is_npls)
    return res_h_is
