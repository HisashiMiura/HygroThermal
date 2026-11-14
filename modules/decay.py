"""木材腐朽：木材の質量減少率と、腐朽で生まれる水（計算方法の文書の評価方法「木材腐朽」）。

各値のコメントは、1 行目に「名前 記号, 単位」（無次元は -）、2 行目に「式 (k)」と補足を書く。
コメントの「式 (k)」は、計算方法の文書の「木材腐朽」（docs/source/contents/03_15_eval_decay.md）の式の番号。

木材腐朽を求める層（入力の wood_decay）の材料の質点で、時間刻みごとに、ステップの終わりの温度・相対湿度から腐朽の状態を進める。
腐朽が始まるまで（発芽期）は Nofal・Kumaran（2011）、腐朽が進む速さ（成長期）と腐朽で生まれる水は齋藤ら（2008）の式による。
腐朽が始まるのは層の表面の質点（層の両端の質点）だけとし、内部へは、同じ層の隣の質点の質量減少率が 0 より大きいときに広がる。
腐朽で生まれる水は、次のステップの水分の収支に足す（「熱・水分の収支」の式 (8)）。
"""

from dataclasses import dataclass

import numpy as np

from modules.air_water_properties import RHO_W
from modules.wall_geometry import MaterialNode, WallGeometry

# 1 か月（発芽までの時間の単位）, s
_SECONDS_PER_MONTH = 30.0 * 86400.0

# 樹種と表面の状態を表す係数 WS, -（式 (2)・(3)。値の根拠は作成者に確かめる）
WS = -1.0

# 発芽までの時間の下限 t_i,min, 月（式 (5)）
_T_INIT_MIN = 0.5

# 発芽までの時間の式の値が負のとき、0 とみなす範囲の下端 t_i,neg, 月（式 (5)。これより小さければ発芽しないとみなす）
_T_INIT_NEGATIVE_LIMIT = -0.59

# 発芽しないとみなすときの、発芽までの時間 t_i,never, 月（式 (5)）
_T_INIT_NEVER = 100.0

# 腐朽が進む相対湿度の下限 φ_g, %（5 節）
RH_GROWTH = 98.0

# 腐朽が進む温度の範囲 θ_g,min・θ_g,max, deg.C（5 節。θ_g,min より高く θ_g,max 以下）
THETA_GROWTH_MIN = 0.0
THETA_GROWTH_MAX = 40.0

# 質量減少率の上限 L_max, -（式 (9)）
MASS_LOSS_MAX = 0.6

# 質量減少率に対する水分生成率 h, (m3/m3)/-（式 (10)）
H_GENERATION = 0.319


def rh_critical(theta: float) -> float:
    """腐朽が始まる相対湿度の下限 φ_c, %（式 (1)）。

    Args:
        theta: 温度 θ, deg.C
    """

    # 腐朽が始まる相対湿度の下限 φ_c, %
    # 式 (1)
    return min(-0.5 * theta + 100.0, 100.0) if theta <= 15.0 else 92.5


def initial_response_time(theta: float, rh: float) -> float:
    """発芽までの時間 t_i, s（式 (2)〜(5)）。

    Args:
        theta: 温度 θ, deg.C
        rh: 相対湿度 φ, %
    """

    # 質量減少率の式の係数 f・g（f は 1/月、g は -）
    # 式 (2)・(3)
    f = 0.1384 * theta + 0.4370 * rh - 42.9450 + WS * (0.0340 * theta - 0.0210 * rh + 1.7210)
    g = -2.2270 * theta - 0.0347 * rh + 0.0244 * theta * rh + WS * (-0.504 * theta + 0.0096 * rh + 0.0047 * theta * rh)

    # 質量減少率が 0 になる時間 t'_i, 月（f = 0 なら、発芽しないとみなせるように負の無限大）
    # 式 (4)
    t_i_dash = -g / f if f != 0.0 else -np.inf

    # 発芽までの時間 t_i, 月：0 以下なら、t_i,neg より大きければ 0、それより小さければ発芽しないとみなす。下限は t_i,min
    # 式 (5)
    if t_i_dash > 0.0:
        t_i = t_i_dash
    elif t_i_dash > _T_INIT_NEGATIVE_LIMIT:
        t_i = 0.0
    else:
        t_i = _T_INIT_NEVER
    t_i = max(t_i, _T_INIT_MIN)

    return t_i * _SECONDS_PER_MONTH


def reaction_rate(theta: float) -> float:
    """腐朽の反応速度定数 k_m, 1/s（式 (8)。負になる温度では 0）。

    Args:
        theta: 温度 θ, deg.C
    """

    # 反応速度定数 k_m, 1/s
    # 式 (8)
    k_m = max((2.77 - 3.23 * theta + 0.865 * theta ** 2 - 0.0189 * theta ** 3) * 1.0e-10, 0.0)
    return k_m


@dataclass(frozen=True, eq=False)
class DecayNodes:
    """木材腐朽を求める質点と、その表面・隣の関係（「木材腐朽」の 3 節）。"""

    # 木材腐朽を求める質点か, [i]
    is_target_is: np.ndarray

    # 層の表面の質点（腐朽が始まりうる質点）か, [i]
    is_surface_is: np.ndarray

    # 同じ層の、左・右の隣の質点の番号, [i]（隣がない・腐朽を求めない質点は -1）
    i_left_is: np.ndarray
    i_right_is: np.ndarray

    # 質点のコントロールボリュームの幅 w_i, m, [i]（腐朽を求めない質点は 0）
    width_is: np.ndarray

    @classmethod
    def from_geometry(cls, geometry: WallGeometry):
        """壁の形から作る（木材腐朽を求める層の材料の質点）。"""

        nodes = geometry.nodes
        n = len(nodes)

        # 木材腐朽を求める質点
        is_target_is = np.array([isinstance(node, MaterialNode) and node.wood_decay for node in nodes])

        # 同じ層の隣の質点（材料の質点は層の中で続いて並び、層と層のあいだには水膜の質点がある）
        def same_layer(i: int, j: int) -> bool:
            return 0 <= j < n and is_target_is[j] and nodes[j].layer_index == nodes[i].layer_index
        i_left_is = np.array([i - 1 if is_target_is[i] and same_layer(i, i - 1) else -1 for i in range(n)])
        i_right_is = np.array([i + 1 if is_target_is[i] and same_layer(i, i + 1) else -1 for i in range(n)])

        # 層の表面の質点：層の端の質点（左か右に、同じ層の隣がない）
        is_surface_is = is_target_is & ((i_left_is < 0) | (i_right_is < 0))

        width_is = np.array([node.width if t else 0.0 for node, t in zip(nodes, is_target_is)])

        return cls(is_target_is=is_target_is, is_surface_is=is_surface_is, i_left_is=i_left_is, i_right_is=i_right_is,
                   width_is=width_is)

    @property
    def any(self) -> bool:
        """木材腐朽を求める質点があるか。"""
        return bool(self.is_target_is.any())


@dataclass(frozen=True, eq=False)
class DecayState:
    """ステップ n の木材腐朽の状態, [i]（腐朽を求めない質点は 0・False のまま）。"""

    # 相対湿度が φ_c 以上の状態が続いた時間 τ_i,n, s（表面の質点だけ）
    time_wet_is: np.ndarray

    # 発芽したか G_i,n（表面の質点だけ。一度発芽したら戻らない）
    germinated_is: np.ndarray

    # 質量減少率 L_i,n, -
    mass_loss_is: np.ndarray

    @classmethod
    def initial(cls, n: int):
        """腐朽していない状態（質点の数 n）。"""
        return cls(time_wet_is=np.zeros(n), germinated_is=np.zeros(n, dtype=bool), mass_loss_is=np.zeros(n))


def advance(
        state_n: DecayState,
        nodes: DecayNodes,
        theta_is_npls: np.ndarray,
        rh_is_npls: np.ndarray,
        dt: float
    ) -> DecayState:
    """木材腐朽の状態を、ステップ n から n + 1 へ進める。

    Args:
        state_n: ステップ n の木材腐朽の状態
        nodes: 木材腐朽を求める質点
        theta_is_npls: ステップ n + 1 の各質点の温度 θ_i,n+1, deg.C, [i]
        rh_is_npls: ステップ n + 1 の各質点の相対湿度 φ_i,n+1, %, [i]
        dt: 時間刻み Δt, s

    Returns:
        ステップ n + 1 の木材腐朽の状態

    Notes:
        計算方法の文書の「木材腐朽」の式 (6)・(7)・(9) と 5 節の条件に対応する。
    """

    time_wet_is = state_n.time_wet_is.copy()
    germinated_is = state_n.germinated_is.copy()
    mass_loss_is = state_n.mass_loss_is.copy()

    for i in np.flatnonzero(nodes.is_target_is):
        theta, rh = theta_is_npls[i], rh_is_npls[i]

        # 発芽期（表面の質点だけ）
        if nodes.is_surface_is[i] and not germinated_is[i]:

            # 相対湿度が φ_c 以上の状態が続いた時間 τ_i,n+1, s（0 deg.C 以下か、φ_c より低ければ 0 に戻す）
            # 式 (6)
            time_wet_is[i] = 0.0 if theta <= 0.0 or rh < rh_critical(theta) else time_wet_is[i] + dt

            # 発芽したか G_i,n+1（続いた時間が、発芽までの時間を超えたら発芽する）
            # 式 (7)
            germinated_is[i] = time_wet_is[i] > initial_response_time(theta, rh)

        # 成長期：表面の質点で発芽しているか、同じ層の隣の質点のステップ n の質量減少率が 0 より大きく、
        # 温度が θ_g,min より高く θ_g,max 以下、相対湿度が φ_g 以上なら、腐朽が進む（5 節）
        started = nodes.is_surface_is[i] and germinated_is[i]
        spread = any(j >= 0 and state_n.mass_loss_is[j] > 0.0 for j in (nodes.i_left_is[i], nodes.i_right_is[i]))
        if (started or spread) and THETA_GROWTH_MIN < theta <= THETA_GROWTH_MAX and rh >= RH_GROWTH:

            # 質量減少率 L_i,n+1, -
            # 式 (9)
            mass_loss_is[i] = min(mass_loss_is[i] + reaction_rate(theta) * dt, MASS_LOSS_MAX)

    return DecayState(time_wet_is=time_wet_is, germinated_is=germinated_is, mass_loss_is=mass_loss_is)


def water_generation_is(
        nodes: DecayNodes,
        state_n: DecayState,
        state_nmns: DecayState,
        dt: float
    ) -> np.ndarray:
    """ステップ n + 1 の水分の収支に足す、腐朽で生まれる水 j_gen,i,n+1, kg/(m2 s), [i]。

    ステップ n - 1 から n までの質量減少率の増えから求める（ステップ n + 1 の未知数によらない）。

    Args:
        nodes: 木材腐朽を求める質点
        state_n: ステップ n の木材腐朽の状態
        state_nmns: ステップ n - 1 の木材腐朽の状態
        dt: 時間刻み Δt, s

    Returns:
        各質点の、腐朽で生まれる水 j_gen,i,n+1, kg/(m2 s), [i]

    Notes:
        計算方法の文書の「木材腐朽」の式 (10) に対応する。
    """

    # 腐朽で生まれる水 j_gen,i,n+1, kg/(m2 s), [i]
    # 式 (10)
    j_gen_is_npls = H_GENERATION * RHO_W * (state_n.mass_loss_is - state_nmns.mass_loss_is) / dt * nodes.width_is
    return j_gen_is_npls
