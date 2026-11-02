"""収支の計算：つながりごとの流れと、質点ごとの残差（設計の文書 2.8 節）。

物理（流れと収支）はここに書く。解く部分（ニュートン法）は物理を知らない汎用の部品にする。
"""

import numpy as np

from modules.boundary import IndoorBoundary, OutdoorBoundary
from modules.state import State
from modules.air_water_properties import C_A, L_V, R_V, get_rho_a
from modules.wall_geometry import AirNode, WallGeometry


def heat_flow_ks(
        geometry: WallGeometry,
        t_is: np.ndarray,
        outdoor: OutdoorBoundary,
        indoor: IndoorBoundary
    ) -> np.ndarray:
    """つながりごとの熱の流れ q_k, W/m2, [k]（外気側から室内側へ流れる向きを正）。

    q_k = (T_k^左 - T_k^右) / R_k。つながり k は質点 k - 1 と質点 k の間（0 番目は外気と質点 0、最後は最後の質点と室内）。
    水蒸気が運ぶ潜熱と液水が運ぶ熱は、熱と水分を一緒に解くとき（段階 5）に足す。

    Args:
        geometry: 壁の形
        t_is: 各質点の温度, K, [i]
        outdoor: 外気の境界条件（外気温 outdoor.state.t を使う）
        indoor: 室内の境界条件（室温 indoor.state.t を使う）
    """

    # 各つながりの熱抵抗 R_k, (m2 K)/W, [k]
    r_h_ks = np.array([link.r_h for link in geometry.links])

    # つながり k の左の温度, K, [k]（0 番目は外気、それ以外は質点 k - 1）
    t_left_ks = np.concatenate(([outdoor.state.t], t_is))

    # つながり k の右の温度, K, [k]（最後は室内、それ以外は質点 k）
    t_right_ks = np.concatenate((t_is, [indoor.state.t]))

    return (t_left_ks - t_right_ks) / r_h_ks


def _state_left_right_ks(
        t_is: np.ndarray,
        mu_is: np.ndarray,
        outdoor: OutdoorBoundary,
        indoor: IndoorBoundary
    ) -> tuple[list[State], list[State]]:
    """つながり k の左の状態と右の状態の並び, [k]（左の 0 番目は外気、右の最後は室内）。

    Args:
        t_is: 各質点の温度, K, [i]
        mu_is: 各質点の水分化学ポテンシャル, J/kg, [i]
        outdoor: 外気の境界条件
        indoor: 室内の境界条件
    """

    # 各質点の 1 点の状態, [i]
    state_is = [State(t=t_i, mu=mu_i) for t_i, mu_i in zip(t_is, mu_is)]

    # つながり k の左は質点 k - 1（0 番目は外気）、右は質点 k（最後は室内）
    return [outdoor.state] + state_is, state_is + [indoor.state]


def vapour_flow_ks(
        geometry: WallGeometry,
        t_is: np.ndarray,
        mu_is: np.ndarray,
        outdoor: OutdoorBoundary,
        indoor: IndoorBoundary
    ) -> np.ndarray:
    """つながりごとの水蒸気の流れ j_v,k, kg/(m2 s), [k]（外気側から室内側へ流れる向きを正）。

    j_v,k = (p_v,左 - p_v,右) / R_m,k（docs/governing_equations.md の式 6.1）。
    水蒸気圧は、質点の温度と水分化学ポテンシャルから p_v = φ p_sat(T)、φ = exp(μ / (R_v T)) で求める（偏微分に分けない）。
    透湿抵抗は、左右それぞれの側の状態（相対湿度）で計算する。

    Args:
        geometry: 壁の形
        t_is: 各質点の温度, K, [i]
        mu_is: 各質点の水分化学ポテンシャル, J/kg, [i]
        outdoor: 外気の境界条件
        indoor: 室内の境界条件
    """

    # つながり k の左の状態と右の状態, [k]（左の 0 番目は外気、右の最後は室内）
    state_left_ks, state_right_ks = _state_left_right_ks(t_is, mu_is, outdoor, indoor)

    # つながり k の左右の水蒸気圧, Pa, [k]
    p_v_left_ks = np.array([s.p_v for s in state_left_ks])
    p_v_right_ks = np.array([s.p_v for s in state_right_ks])

    # 各つながりの透湿抵抗 R_m,k, (m2 s Pa)/kg, [k]（左右それぞれの側の状態で計算する）
    r_m_ks = np.array([
        link.r_m(s_left, s_right)
        for link, s_left, s_right in zip(geometry.links, state_left_ks, state_right_ks)
    ])

    return (p_v_left_ks - p_v_right_ks) / r_m_ks


def liquid_flow_ks(
        geometry: WallGeometry,
        t_is: np.ndarray,
        mu_is: np.ndarray,
        outdoor: OutdoorBoundary,
        indoor: IndoorBoundary
    ) -> np.ndarray:
    """つながりごとの液水の流れ j_l,k, kg/(m2 s), [k]（外気側から室内側へ流れる向きを正）。

    j_l,k = (μ_左 - μ_右) / R_l,k（docs/governing_equations.md の式 6.1）。
    境目（表面・層と層・材料と通気層）と、液相水分伝導率が 0 の所では、液水の抵抗が無限大なので流れは 0。

    Args:
        geometry: 壁の形
        t_is: 各質点の温度, K, [i]
        mu_is: 各質点の水分化学ポテンシャル, J/kg, [i]
        outdoor: 外気の境界条件
        indoor: 室内の境界条件
    """

    # つながり k の左の状態と右の状態, [k]
    state_left_ks, state_right_ks = _state_left_right_ks(t_is, mu_is, outdoor, indoor)

    # つながり k の左右の水分化学ポテンシャル, J/kg, [k]
    mu_left_ks = np.array([s.mu for s in state_left_ks])
    mu_right_ks = np.array([s.mu for s in state_right_ks])

    # 各つながりの液水の抵抗 R_l,k, (m2 s (J/kg))/kg, [k]（境目と、液相水分伝導率が 0 の所では無限大）
    r_l_ks = np.array([
        link.r_l(s_left, s_right)
        for link, s_left, s_right in zip(geometry.links, state_left_ks, state_right_ks)
    ])

    # 抵抗が無限大のつながりは、差 / 無限大 = 0
    return (mu_left_ks - mu_right_ks) / r_l_ks


def water_amount_is(geometry: WallGeometry, t_is: np.ndarray, mu_is: np.ndarray) -> np.ndarray:
    """各質点の水分量 W_i, kg/m2, [i]（docs/governing_equations.md 6.2 節）。

    Args:
        geometry: 壁の形
        t_is: 各質点の温度, K, [i]
        mu_is: 各質点の水分化学ポテンシャル, J/kg, [i]
    """
    return np.array([
        node.water_amount(State(t=t_i, mu=mu_i))
        for node, t_i, mu_i in zip(geometry.nodes, t_is, mu_is)
    ])


def heat_capacity_is(geometry: WallGeometry, t_is: np.ndarray, mu_is: np.ndarray) -> np.ndarray:
    """各質点の熱容量 C_i, J/(m2 K), [i]（乾燥した材料 ＋ 含まれる水。docs/governing_equations.md 6.2 節）。

    Args:
        geometry: 壁の形
        t_is: 各質点の温度, K, [i]
        mu_is: 各質点の水分化学ポテンシャル, J/kg, [i]
    """
    return np.array([
        node.heat_capacity(State(t=t_i, mu=mu_i))
        for node, t_i, mu_i in zip(geometry.nodes, t_is, mu_is)
    ])


def ventilation_rate_is(geometry: WallGeometry, t_is: np.ndarray, outdoor: OutdoorBoundary) -> np.ndarray:
    """各質点の、壁 1 m2 あたりの換気量 q_v,i, m3/(m2 s), [i]（通気層の空気の質点だけ。材料の質点は 0）。

    上流の空気は外気とする（壁が 1 枚なので）。

    Args:
        geometry: 壁の形
        t_is: 各質点の温度, K, [i]（solver.step では、前の時刻の温度を渡す）
        outdoor: 外気の境界条件
    """
    return np.array([
        node.ventilation_rate(t_i, outdoor.state.t, geometry.height, geometry.len_vertical) if isinstance(node, AirNode) else 0.0
        for node, t_i in zip(geometry.nodes, t_is)
    ])


def moisture_residual_is(
        geometry: WallGeometry,
        t_is_npls: np.ndarray,
        mu_is_npls: np.ndarray,
        t_is_n: np.ndarray,
        mu_is_n: np.ndarray,
        dt: float,
        outdoor: OutdoorBoundary,
        indoor: IndoorBoundary,
        q_v_is: np.ndarray | None = None
    ) -> np.ndarray:
    """質点ごとの水分の収支の残差 res_m,i, kg/(m2 s), [i]。

    res_m,i = (W_i,n+1 - W_i,n) / Δt - (j_i - j_(i+1) + j_vent,i)（たまる量 -（左のつながりの流れ - 右のつながりの流れ + 換気で入る水分））。
    水分の流れ j = 水蒸気の流れ j_v + 液水の流れ j_l。流れはすべてステップ n + 1 の値で計算する（後退差分）。
    たまる量を水分量の変化で書くので、水分量が厳密に保存される。
    換気で入る水分 j_vent,i = q_v,i (ρ_v,ex - ρ_v,i)（水蒸気の密度 ρ_v = p_v / (R_v T) の差。通気層の空気の質点だけ）。
    正しいステップ n + 1 の温度・水分化学ポテンシャルを入れると 0 になる（docs/governing_equations.md の式 6.2）。

    Args:
        geometry: 壁の形
        t_is_npls: 各質点のステップ n + 1 の温度（仮の値）, K, [i]
        mu_is_npls: 各質点のステップ n + 1 の水分化学ポテンシャル（仮の値）, J/kg, [i]
        t_is_n: 各質点のステップ n の温度, K, [i]
        mu_is_n: 各質点のステップ n の水分化学ポテンシャル, J/kg, [i]
        dt: 時間刻み, s
        outdoor: ステップ n + 1 の外気の境界条件
        indoor: ステップ n + 1 の室内の境界条件
        q_v_is: 各質点の、壁 1 m2 あたりの換気量, m3/(m2 s), [i]（ventilation_rate_is。省くと換気なし）
    """

    # 換気量を省いたときは、換気なし（どの質点も 0）
    if q_v_is is None:
        q_v_is = np.zeros(len(t_is_npls))

    # ステップ n + 1 とステップ n の水分量, kg/m2, [i]
    w_is_npls = water_amount_is(geometry, t_is_npls, mu_is_npls)
    w_is_n = water_amount_is(geometry, t_is_n, mu_is_n)

    # 各つながりの水分の流れ（水蒸気 + 液水）, kg/(m2 s), [k]（ステップ n + 1 の値）
    j_ks = (vapour_flow_ks(geometry, t_is_npls, mu_is_npls, outdoor, indoor)
            + liquid_flow_ks(geometry, t_is_npls, mu_is_npls, outdoor, indoor))

    # 質点 i の左のつながり（つながり i）の流れ, kg/(m2 s), [i]（j_ks の最後を除く。正なら質点 i に入ってくる）
    j_left_is = j_ks[:-1]

    # 質点 i の右のつながり（つながり i + 1）の流れ, kg/(m2 s), [i]（j_ks の先頭を除く。正なら質点 i から出ていく）
    j_right_is = j_ks[1:]

    # 外気と各質点の水蒸気の密度, kg/m3（ステップ n + 1 の値。「湿り空気と水の性質」の式 (5)）
    rho_v_ex = outdoor.state.p_v / (R_V * outdoor.state.t)
    rho_v_is = np.array([State(t=t_i, mu=mu_i).p_v for t_i, mu_i in zip(t_is_npls, mu_is_npls)]) / (R_V * t_is_npls)

    # 換気で入る水分, kg/(m2 s), [i]（換気量が 0 の質点は 0）
    j_vent_is = q_v_is * (rho_v_ex - rho_v_is)

    # たまる量 -（左のつながりの流れ - 右のつながりの流れ + 換気で入る水分）
    return (w_is_npls - w_is_n) / dt - (j_left_is - j_right_is + j_vent_is)


def heat_residual_is(
        geometry: WallGeometry,
        t_is_npls: np.ndarray,
        mu_is_npls: np.ndarray,
        t_is_n: np.ndarray,
        dt: float,
        outdoor: OutdoorBoundary,
        indoor: IndoorBoundary,
        q_sol: float = 0.0,
        q_v_is: np.ndarray | None = None
    ) -> np.ndarray:
    """質点ごとの熱の収支の残差 res_h,i, W/m2, [i]（docs/governing_equations.md の式 6.3）。

    res_h,i = C_i,n+1 (T_i,n+1 - T_i,n) / Δt - (q_i - q_(i+1) + δ_i L_v (j_v,i - j_v,(i+1)) + q_sol,i + q_vent,i)
    （たまる量 -（左のつながりの熱伝導 - 右のつながりの熱伝導 + 凝縮熱 + 外から直接入る熱 + 換気で入る熱））。
    流れはすべてステップ n + 1 の値で計算する（後退差分）。
    凝縮熱 L_v (j_v,i - j_v,(i+1)) は、その質点の凝縮量（蒸発量）を水蒸気の収支から見積もったもの（材料では、水蒸気のまま
    たまる量はわずかなので、凝縮量 ≈ 左から入ってくる水蒸気 - 右へ出ていく水蒸気）。L_v は蒸発潜熱。
    δ_i は材料の質点で 1、通気層の空気の質点で 0。空気の質点では、たまる水分はすべて水蒸気のままで相変化は起きないので、
    凝縮熱を数えない（換気で入った水蒸気がつながりから出ていくとき、通気層で蒸発したと見積もらないため）。
    C_i は水の分を足した熱容量で、ステップ n + 1 の仮の値で計算する。移動する水が運ぶ顕熱は入れない（量が非常に小さい）。
    外表面の質点（0 番目）には、外から直接入る熱 q_sol（日射の吸収 - 夜間放射）も入る。
    通気層の空気の質点には、換気で入る熱 q_vent,i = cρ q_v,i (T_ex - T_i) も入る（上流の空気は外気）。
    正しいステップ n + 1 の温度・水分化学ポテンシャルを入れると 0 になる。

    Args:
        geometry: 壁の形
        t_is_npls: 各質点のステップ n + 1 の温度（仮の値）, K, [i]
        mu_is_npls: 各質点のステップ n + 1 の水分化学ポテンシャル（仮の値）, J/kg, [i]
        t_is_n: 各質点のステップ n の温度, K, [i]
        dt: 時間刻み, s
        outdoor: ステップ n + 1 の外気の境界条件
        indoor: ステップ n + 1 の室内の境界条件
        q_sol: 外表面に外から直接入る熱, W/m2（surface_radiation.get_solar_heat_gain）
        q_v_is: 各質点の、壁 1 m2 あたりの換気量, m3/(m2 s), [i]（ventilation_rate_is。省くと換気なし）
    """

    # 換気量を省いたときは、換気なし（どの質点も 0）
    if q_v_is is None:
        q_v_is = np.zeros(len(t_is_npls))

    # 各質点のステップ n + 1 の熱容量 C_i, J/(m2 K), [i]（乾燥した材料 + 含まれる水）
    c_is_npls = heat_capacity_is(geometry, t_is_npls, mu_is_npls)

    # 各つながりの熱伝導の流れと水蒸気の流れ, W/m2・kg/(m2 s), [k]（ステップ n + 1 の値）
    q_ks = heat_flow_ks(geometry, t_is_npls, outdoor, indoor)
    j_v_ks = vapour_flow_ks(geometry, t_is_npls, mu_is_npls, outdoor, indoor)

    # 質点 i の左のつながり（つながり i）の熱伝導の流れ, W/m2, [i]（q_ks の最後を除く。正なら質点 i に入ってくる）
    q_left_is = q_ks[:-1]

    # 質点 i の右のつながり（つながり i + 1）の熱伝導の流れ, W/m2, [i]（q_ks の先頭を除く。正なら質点 i から出ていく）
    q_right_is = q_ks[1:]

    # 凝縮熱, W/m2, [i]（凝縮量 ≈ 左から入ってくる水蒸気 - 右へ出ていく水蒸気。通気層の空気の質点は相変化しないので 0）
    is_material_is = np.array([not isinstance(node, AirNode) for node in geometry.nodes])
    q_cond_is = np.where(is_material_is, L_V * (j_v_ks[:-1] - j_v_ks[1:]), 0.0)

    # 外から直接入る熱（日射の吸収 - 夜間放射）, W/m2, [i]（外表面の質点 0 だけ）
    q_sol_is = np.zeros(len(t_is_npls))
    q_sol_is[0] = q_sol

    # 換気で入る熱, W/m2, [i]（換気量が 0 の質点は 0。空気の密度は、入ってくる外気の温度で求める）
    q_vent_is = C_A * get_rho_a(outdoor.state.t) * q_v_is * (outdoor.state.t - t_is_npls)

    # たまる量 -（左のつながりの熱伝導 - 右のつながりの熱伝導 + 凝縮熱 + 外から直接入る熱 + 換気で入る熱）
    return c_is_npls * (t_is_npls - t_is_n) / dt - (q_left_is - q_right_is + q_cond_is + q_sol_is + q_vent_is)
