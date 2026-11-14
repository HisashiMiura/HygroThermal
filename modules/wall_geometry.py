"""壁の形（時間で変わらないもの）：質点とつながり。

壁の入力（wall_input.WallInput）から作る。

計算方法の文書の評価方法との対応：質点の分け方とつながり（_make_nodes・_make_links）は「壁の構成と質点」、
つながりの抵抗（Link・Interface、質点の r_h_left など）は「熱・水分の流れ」、水分量・熱容量（water_amount・heat_capacity）と水膜の質点（FilmNode）は「熱・水分の収支」、
換気量は「通気層の換気」。コメントの「式 (k)」は、とくに断りがなければ「壁の構成と質点」
（docs/source/contents/03_04_eval_wall_nodes.md）の式の番号。

質点の置き方（設計の文書 2.3 節）：層の両端（表面・異種材料の境目）に必ず質点を置く。
層の中の質点はコントロールボリュームの中心にあり、層の端の質点は表面・境目の上にある。

つながりは、質点に「質点からコントロールボリュームの端までの抵抗」を聞いて、境目の抵抗と足す。
材料の質点・空気の質点・水膜の質点は、同じ名前（width・heat_capacity・water_amount・r_h_left・r_h_right・layer_index）で使える。
層と層の境目には、いつも水膜の質点（厚さのない質点）を置く。
材料の質点と水膜の質点の状態の水分化学ポテンシャルの欄には状態量 s を入れ、0 を超えた分で自由水・水膜の水を表す（mu で水分化学ポテンシャルに直す）。
"""

from collections.abc import Callable
from dataclasses import dataclass
from functools import cached_property
from typing import Self

import numpy as np

from modules.direction import Direction
from modules.material_loaders import load_materials
from modules.materials import Material, Materials
from modules.state import State
from modules.air_water_properties import C_A, C_W, R_V, RHO_W, get_p_vs, get_rho_a
from modules.wall_input import AirLayerInput, FilmInterfaceInput, InterfaceInput, LayerInput, WallInput


# 重力加速度, m/s2
_G = 9.8

# 水膜の係数 k_i を求める基準の温度 T_ref, K（「熱・水分の収支」の式 (12)）
_FILM_T_REF = 293.15

# 水膜の係数 k_i を求める基準の時間 τ, s（「熱・水分の収支」の式 (12)）
_FILM_TAU = 3600.0

# 水膜が流れ落ちずに留まれる水分量 W_f,max, kg/m2（「熱・水分の収支」の式 (13)）
FILM_W_MAX = 0.005

# 水膜の留まれる量を超えた水が流れ落ちる時定数 τ_d, s（「熱・水分の収支」の式 (13)）
FILM_TAU_DRAIN = 600.0

# 材料の自由水の係数 k_i を求める基準の温度 T_ref, K（「熱・水分の収支」の式 (19)）
_FREE_WATER_T_REF = 293.15

# 平衡含水率の、相対湿度 100% の手前の傾きを求める相対湿度の比の差 Δr, -（「熱・水分の収支」の式 (18)）
_FREE_WATER_DR = 1.0e-4

# 材料の自由水の係数 k_i を求めるときの、乾燥密度 × 平衡含水率の傾きの下限 (ρ_s u'_1)_min, kg/m3
# （吸湿しない材料（平衡含水率が 0）でも、自由水を表せるようにする。「熱・水分の収支」の式 (19)）
_FREE_WATER_RHO_DU_MIN = 1.0


@dataclass(frozen=True, kw_only=True)
class MaterialNode:
    """材料の質点。"""

    # 質点からコントロールボリュームの左（外気側）の端までの距離 d_L,i, m（層の左の端の質点は 0。式 (2)）
    d_left: float

    # 質点からコントロールボリュームの右（室内側）の端までの距離 d_R,i, m（層の右の端の質点は 0。式 (3)）
    d_right: float

    # 材料
    material: Material

    # どの層の質点か（外気側から 0, 1, 2, …）
    layer_index: int

    # 空隙率 ψ, -（None なら、自由水の上限を確かめない。「熱・水分の収支」の式 (20)）
    porosity: float | None = None

    # 木材腐朽を求めるか（「木材腐朽」）
    wood_decay: bool = False

    def __post_init__(self):
        """作るときに値を確かめる。"""

        # 左右の端までの距離が 0 以上でなければエラー（層の端の質点は、表面・境目の側が 0）
        if not (self.d_left >= 0.0 and self.d_right >= 0.0):
            raise ValueError(f'左右の端までの距離 d_left・d_right は 0 以上にしてください（{self.d_left}, {self.d_right}）。')

        # 幅（左右の距離の和）が 0 より大きくなければエラー
        if not self.width > 0.0:
            raise ValueError(f'質点の幅（d_left + d_right）は 0 より大きくしてください（{self.width}）。')

    @property
    def width(self) -> float:
        """コントロールボリュームの幅 w_i, m（左右の端までの距離の和。式 (4)）"""
        return self.d_left + self.d_right

    @property
    def heat_capacity_dry(self) -> float:
        """乾燥した材料の熱容量, J/(m2 K)（幅 × 容積比熱 c_s ρ_s。「熱・水分の収支」の式 (3) の第 1 項）"""
        return self.width * self.material.c_rho

    @staticmethod
    def mu(s: float) -> float:
        """状態量 s から、水分化学ポテンシャル μ, J/kg（自由水があれば 0。「熱・水分の収支」の式 (16)）"""

        # 水分化学ポテンシャル μ_i,n, J/kg
        # 式 (16)（「熱・水分の収支」）
        mu_i_n = min(s, 0.0)
        return mu_i_n

    @cached_property
    def k(self) -> float:
        """状態量を自由水の量に直す係数 k_i, (kg/m2)/(J/kg)（「熱・水分の収支」の式 (18)・(19)）

        目盛りの換算で、計算の結果には影響しない。相対湿度 100% の前後で、水分量の傾きがつながるように決める。
        材料と幅で決まるので、1 回だけ求めて覚えておく。
        """

        # 平衡含水率の、相対湿度 100% の手前の傾き u'_1, kg/kg
        # 式 (18)（「熱・水分の収支」）
        du_1 = (self.material.get_u(100.0) - self.material.get_u(100.0 * (1.0 - _FREE_WATER_DR))) / _FREE_WATER_DR

        # 係数 k_i, (kg/m2)/(J/kg)
        # 式 (19)（「熱・水分の収支」）
        k_i = max(self.material.rho_s * du_1, _FREE_WATER_RHO_DU_MIN) * self.width / (R_V * _FREE_WATER_T_REF)
        return k_i

    def free_water(self, state: State) -> float:
        """自由水の量 W_fw,i,n, kg/m2（平衡含水率の相対湿度 100% の値を超えてたまる水。「熱・水分の収支」の式 (17)）

        Args:
            state: この質点のステップ n の状態（mu の欄に状態量 s_i,n が入っている）
        """

        # 自由水の量 W_fw,i,n, kg/m2（自由水がない s ≤ 0 では 0）
        # 式 (17)（「熱・水分の収支」）
        w_fw_i_n = self.k * max(state.mu, 0.0)
        return w_fw_i_n

    @property
    def water_limit(self) -> float | None:
        """空隙を水で満たす水分量 W_max,i, kg/m2（空隙率 × 水の密度 × 幅。空隙率がなければ None。「熱・水分の収支」の式 (20)）"""

        # 空隙を水で満たす水分量 W_max,i, kg/m2
        # 式 (20)（「熱・水分の収支」）
        return None if self.porosity is None else self.porosity * RHO_W * self.width

    def water_amount(self, state: State) -> float:
        """水分量 W_i,n, kg/m2（乾燥密度 × 質量含水率 × 幅 ＋ 自由水の量。「熱・水分の収支」の式 (1)）

        Args:
            state: この質点のステップ n の状態（mu の欄に状態量 s_i,n が入っている）
        """

        # 相対湿度 φ_i,n, %（自由水がある s > 0 では、水分化学ポテンシャルが 0 なので 100%）
        # 式 (16)（「熱・水分の収支」）
        rh_i_n = 100.0 if state.mu > 0.0 else state.rh

        # 水分量 W_i,n, kg/m2
        # 式 (1)（「熱・水分の収支」）
        w_i_n = self.material.rho_s * self.material.get_u(rh_i_n) * self.width + self.free_water(state)
        return w_i_n

    def heat_capacity(self, state: State) -> float:
        """熱容量 C_i,n, J/(m2 K)（乾燥した材料 ＋ 含まれる水。「熱・水分の収支」の式 (3)）

        Args:
            state: この質点のステップ n の状態
        """

        # 熱容量 C_i,n, J/(m2 K)
        # 式 (3)（「熱・水分の収支」）
        c_i_n = self.heat_capacity_dry + C_W * self.water_amount(state)
        return c_i_n

    @property
    def r_h_left(self) -> float:
        """質点から左（外気側）の端までの熱抵抗 r_h,L,i, (m2 K)/W（「熱・水分の流れ」の式 (3)）"""
        return self.d_left / self.material.lambda_h

    @property
    def r_h_right(self) -> float:
        """質点から右（室内側）の端までの熱抵抗 r_h,R,i, (m2 K)/W（「熱・水分の流れ」の式 (3)）"""
        return self.d_right / self.material.lambda_h

    def r_m_left(self, state: State) -> float:
        """質点から左（外気側）の端までの透湿抵抗 r_m,L,i,n, (m2 s Pa)/kg（「熱・水分の流れ」の式 (4)。その状態の相対湿度での湿気伝導率を使う）"""
        return self.d_left / self.material.get_lambda_dsh_m(rh=state.rh)

    def r_m_right(self, state: State) -> float:
        """質点から右（室内側）の端までの透湿抵抗 r_m,R,i,n, (m2 s Pa)/kg（「熱・水分の流れ」の式 (4)。その状態の相対湿度での湿気伝導率を使う）"""
        return self.d_right / self.material.get_lambda_dsh_m(rh=state.rh)

    def r_l_left(self, state: State) -> float:
        """質点から左（外気側）の端までの液水の抵抗 r_l,L,i,n, (m2 s (J/kg))/kg（「熱・水分の流れ」の式 (5)）"""
        return self._r_l(self.d_left, state)

    def r_l_right(self, state: State) -> float:
        """質点から右（室内側）の端までの液水の抵抗 r_l,R,i,n, (m2 s (J/kg))/kg（「熱・水分の流れ」の式 (5)）"""
        return self._r_l(self.d_right, state)

    def _r_l(self, d: float, state: State) -> float:
        """距離 d の液水の抵抗, (m2 s (J/kg))/kg（「熱・水分の流れ」の式 (5)）。

        液水伝導率が 0（液水が移動しない材料、または液水拡散係数の範囲の外）なら、0 で割らずに無限大にする。
        """

        lambda_l = self.material.get_lambda_dsh_mu_l(rh=state.rh, t=state.t)

        # 液水伝導率が 0 なら無限大（液水は移動しない）
        return float('inf') if lambda_l == 0.0 else d / lambda_l


@dataclass(frozen=True, kw_only=True)
class AirNode:
    """通気層の空気の質点（よく混ざった空気の 1 つの質点）。

    空気はよく混ざっているとみなし、質点から端までの抵抗は 0。
    材料と通気層の間の伝達は、つながりの境目の抵抗が持つ。
    """

    # 通気層の厚さ, m
    thickness: float

    # 流量係数, -（相当開口面積 = α × 厚さ × 壁の横の長さ。壁 1 m2 あたりの換気量では横の長さが消えるので、α のまま持つ）
    alpha: float

    # どの層の質点か（外気側から 0, 1, 2, …）
    layer_index: int

    def __post_init__(self):
        """作るときに値を確かめる。"""

        # 厚さが 0 より大きくなければエラー
        if not self.thickness > 0.0:
            raise ValueError(f'通気層の厚さ thickness は 0 より大きくしてください（{self.thickness}）。')

        # 流量係数が 0 以上でなければエラー（0 は換気のない空気層）
        if not self.alpha >= 0.0:
            raise ValueError(f'流量係数 alpha は 0 以上にしてください（{self.alpha}）。')

    @property
    def width(self) -> float:
        """コントロールボリュームの幅 w_i, m（通気層の厚さ d_l。「壁の構成と質点」の 4 節）"""
        return self.thickness

    def water_amount(self, state: State) -> float:
        """水分量 W_i,n, kg/m2（空気の中の水蒸気：水蒸気の密度 × 厚さ。「熱・水分の収支」の式 (2)）

        Args:
            state: この質点のステップ n の状態
        """

        # 水蒸気の密度 ρ_v,i,n, kg/m3
        # 「湿り空気と水の性質」の式 (5)
        rho_v_i_n = state.p_v / (R_V * state.t)

        # 水分量 W_i,n, kg/m2
        # 式 (2)（「熱・水分の収支」）
        w_i_n = rho_v_i_n * self.thickness
        return w_i_n

    def ventilation_rate(
            self,
            t_i_n: float,
            t_ex_n: float,
            height: float,
            len_vertical: float
        ) -> float:
        """ステップ n の、壁 1 m2 あたりの換気量 q_v,i,n, m3/(m2 s)（温度差による浮力換気）。

        上流の空気は外気とする（壁が 1 枚なので）。

        Args:
            t_i_n: ステップ n の通気層の空気の温度 T_i,n, K
            t_ex_n: ステップ n の外気温 T_ex,n, K
            height: 壁の下端と上端の高さの差 H, m（WallGeometry.height）
            len_vertical: 壁の縦の長さ L_vert, m

        Returns:
            ステップ n の、壁 1 m2 あたりの換気量 q_v,i,n, m3/(m2 s)

        Notes:
            計算方法の文書の「通気層の換気」の式 (2) に対応する。
            壁 1 m2 あたりにすると壁の横の長さが消えることは、計算方法の文書の根拠「通気層の換気のモデル」の 5 節に書いた。
            Fortran の Whamcavity の式には外気の密度で割る所がない（作成者に確かめること）。
        """

        # ステップ n の外気と通気層の空気の密度 ρ_a,ex,n・ρ_a,i,n, kg/m3
        # 「湿り空気と水の性質」の式 (6)
        rho_a_ex_n = get_rho_a(t=t_ex_n)
        rho_a_i_n = get_rho_a(t=t_i_n)

        # ステップ n の、壁 1 m2 あたりの換気量 q_v,i,n, m3/(m2 s)
        # 式 (2)（「通気層の換気」）
        q_v_i_n = self.alpha * self.thickness / len_vertical * (2.0 * _G * height * abs(rho_a_ex_n - rho_a_i_n) / rho_a_ex_n) ** 0.5
        return q_v_i_n

    def heat_capacity(self, state: State) -> float:
        """熱容量 C_i,n, J/(m2 K)（空気の比熱 × 空気の密度 × 厚さ。「熱・水分の収支」の式 (4)）

        水蒸気の分は小さいので足さない。

        Args:
            state: この質点のステップ n の状態
        """

        # 空気の密度 ρ_a,i,n, kg/m3（この質点の温度で求める）
        # 「湿り空気と水の性質」の式 (6)
        rho_a_i_n = get_rho_a(state.t)

        # 熱容量 C_i,n, J/(m2 K)
        # 式 (4)（「熱・水分の収支」）
        c_i_n = C_A * rho_a_i_n * self.thickness
        return c_i_n

    @property
    def r_h_left(self) -> float:
        """質点から左の端までの熱抵抗, (m2 K)/W（よく混ざった空気なので 0）"""
        return 0.0

    @property
    def r_h_right(self) -> float:
        """質点から右の端までの熱抵抗, (m2 K)/W（よく混ざった空気なので 0）"""
        return 0.0

    def r_m_left(self, state: State) -> float:
        """質点から左の端までの透湿抵抗, (m2 s Pa)/kg（よく混ざった空気なので 0）"""
        return 0.0

    def r_m_right(self, state: State) -> float:
        """質点から右の端までの透湿抵抗, (m2 s Pa)/kg（よく混ざった空気なので 0）"""
        return 0.0


@dataclass(frozen=True, kw_only=True)
class FilmNode:
    """水膜の質点：層と層の境目の上に置く、厚さのない質点。境目にたまる水（結露した水・浸入した雨水）を受け持つ。

    状態の水分化学ポテンシャルの欄には、水膜の状態量 s_i を入れる（負なら水膜なし、0 以上なら水膜あり）。
    水分化学ポテンシャルは min(s, 0)、水膜の水分量は k_i × max(s, 0)。熱容量は 0（熱の収支は定常）。
    厚さがないので、質点から端までの抵抗は 0。左右の境目のコンダクタンスは、つながりの境目が持つ。
    """

    # どの層の質点か（境目の左の層の番号。「壁の構成と質点」の 4 節）
    layer_index: int

    # 水膜の左（外気側）のつながりの、境目の湿気コンダクタンス G_m,L,i, kg/(m2 s Pa)
    g_m_left: float

    # 水膜の右（室内側）のつながりの、境目の湿気コンダクタンス G_m,R,i, kg/(m2 s Pa)
    g_m_right: float

    # 外表面に当たる雨水のうち、この水膜に浸入する割合 f_p,i, -
    rain_fraction: float = 0.0

    @property
    def width(self) -> float:
        """コントロールボリュームの幅 w_i, m（厚さがないので 0。「壁の構成と質点」の 4 節）"""
        return 0.0

    @property
    def k(self) -> float:
        """水膜の状態量を水膜の水分量に直す係数 k_i, (kg/m2)/(J/kg)（「熱・水分の収支」の式 (12)）

        目盛りの換算で、計算の結果には影響しない。水膜なしと水膜ありで、水分の収支の傾きの大きさがそろうように決める。
        """

        # 係数 k_i, (kg/m2)/(J/kg)
        # 式 (12)（「熱・水分の収支」）
        k_i = _FILM_TAU * (self.g_m_left + self.g_m_right) * get_p_vs(_FILM_T_REF) / (R_V * _FILM_T_REF)
        return k_i

    def water_amount(self, state: State) -> float:
        """水分量（水膜の水分量） W_i,n, kg/m2（「熱・水分の収支」の式 (11)）

        Args:
            state: この質点のステップ n の状態（mu の欄に水膜の状態量 s_i,n が入っている）
        """

        # 水膜の水分量 W_i,n, kg/m2（水膜がない s < 0 では 0）
        # 式 (11)（「熱・水分の収支」）
        w_i_n = self.k * max(state.mu, 0.0)
        return w_i_n

    def heat_capacity(self, state: State) -> float:
        """熱容量 C_i,n, J/(m2 K)（厚さがないので 0。「熱・水分の収支」の 6 節）"""
        return 0.0

    def drainage(self, state: State) -> float:
        """流れ落ちる水 j_drain,i,n, kg/(m2 s)（留まれる量を超えた分が、時定数 τ_d で壁の外へ流れ出す。「熱・水分の収支」の式 (13)）

        Args:
            state: この質点のステップ n の状態（mu の欄に水膜の状態量 s_i,n が入っている）
        """

        # 流れ落ちる水 j_drain,i,n, kg/(m2 s)
        # 式 (13)（「熱・水分の収支」）
        j_drain_i_n = max(self.water_amount(state) - FILM_W_MAX, 0.0) / FILM_TAU_DRAIN
        return j_drain_i_n

    @staticmethod
    def mu(s: float) -> float:
        """水膜の状態量 s から、水分化学ポテンシャル μ, J/kg（水膜があれば 0。「熱・水分の収支」の式 (10)）"""

        # 水分化学ポテンシャル μ_i,n, J/kg
        # 式 (10)（「熱・水分の収支」）
        mu_i_n = min(s, 0.0)
        return mu_i_n

    @property
    def r_h_left(self) -> float:
        """質点から左の端までの熱抵抗, (m2 K)/W（厚さがないので 0）"""
        return 0.0

    @property
    def r_h_right(self) -> float:
        """質点から右の端までの熱抵抗, (m2 K)/W（厚さがないので 0）"""
        return 0.0

    def r_m_left(self, state: State) -> float:
        """質点から左の端までの透湿抵抗, (m2 s Pa)/kg（厚さがないので 0）"""
        return 0.0

    def r_m_right(self, state: State) -> float:
        """質点から右の端までの透湿抵抗, (m2 s Pa)/kg（厚さがないので 0）"""
        return 0.0


# 質点：材料の質点・空気の質点・水膜の質点のどれか
Node = MaterialNode | AirNode | FilmNode


@dataclass(frozen=True)
class Outdoor:
    """外気。つながりの左（外気側）の端になる。

    質点ではないので、端までの抵抗は 0（外表面の伝達は、つながりの境目が持つ）。
    外気の温度・湿度は、時刻ごとの境界条件から取り出す（段階 6）。
    """

    @property
    def r_h_right(self) -> float:
        """右の端までの熱抵抗, (m2 K)/W（質点ではないので 0）"""
        return 0.0

    def r_m_right(self, state: State) -> float:
        """右の端までの透湿抵抗, (m2 s Pa)/kg（質点ではないので 0）"""
        return 0.0


@dataclass(frozen=True)
class Indoor:
    """室内。つながりの右（室内側）の端になる。

    質点ではないので、端までの抵抗は 0（内表面の伝達は、つながりの境目が持つ）。
    室内の温度・湿度は、時刻ごとの境界条件から取り出す（段階 6）。
    """

    @property
    def r_h_left(self) -> float:
        """左の端までの熱抵抗, (m2 K)/W（質点ではないので 0）"""
        return 0.0

    def r_m_left(self, state: State) -> float:
        """左の端までの透湿抵抗, (m2 s Pa)/kg（質点ではないので 0）"""
        return 0.0


# 外気と室内（値を持たないので、1 つずつ作って使い回す）
OUTDOOR = Outdoor()
INDOOR = Indoor()


@dataclass(frozen=True, kw_only=True)
class Interface:
    """計算の側の境目：表面の伝達、層と層の接触、材料と通気層の間の伝達、塗膜・防湿層などの抵抗。

    入力の境目（wall_input.InterfaceInput、コンダクタンス）から、_make_links で作る。
    境目では液水は移動しない。
    """

    # 境目の熱抵抗 R_h,int,k, (m2 K)/W（「熱・水分の流れ」の式 (1)）
    r_h: float

    # 境目の透湿抵抗 R_m,int,k, (m2 s Pa)/kg（「熱・水分の流れ」の式 (2)）
    r_m: float

    def __post_init__(self):
        """作るときに値を確かめる。"""

        # 熱抵抗・透湿抵抗が 0 より大きくなければエラー
        for key in ('r_h', 'r_m'):
            value = getattr(self, key)
            if not value > 0.0:
                raise ValueError(f'境目の {key} は 0 より大きくしてください（{value}）。')


@dataclass(frozen=True, kw_only=True)
class Link:
    """つながり：隣り合う 2 つの質点、または端の質点と外気・室内の間。

    抵抗 = 左から右の端までの抵抗 + 境目の抵抗 + 右から左の端までの抵抗。
    層の中のつながりは境目を持たない。層と層・材料と通気層・表面のつながりは、境目（Interface）を持つ。
    """

    # 左（外気側）：質点か外気
    left: Node | Outdoor

    # 右（室内側）：質点か室内
    right: Node | Indoor

    # 境目。層の中のつながりは None（境目がない）
    interface: Interface | None

    def __post_init__(self):
        """作るときに、つながりの形を確かめる。"""

        # 左が質点か外気でなければエラー（室内を左に置くことはできない）
        if not isinstance(self.left, (MaterialNode, AirNode, FilmNode, Outdoor)):
            raise ValueError(f'つながりの左 left は、質点か外気にしてください（{type(self.left).__name__}）。')

        # 右が質点か室内でなければエラー（外気を右に置くことはできない）
        if not isinstance(self.right, (MaterialNode, AirNode, FilmNode, Indoor)):
            raise ValueError(f'つながりの右 right は、質点か室内にしてください（{type(self.right).__name__}）。')

        # 外気・室内に面するつながりに境目がなければエラー（表面の伝達を境目で表すため）
        if isinstance(self.left, Outdoor) or isinstance(self.right, Indoor):
            if self.interface is None:
                raise ValueError('外気・室内に面するつながりには、境目 interface（表面の伝達）が必要です。')
            return

        # 水膜の質点に面するつながりに境目がなければエラー（水膜の左右のコンダクタンスを境目で表すため）
        if isinstance(self.left, FilmNode) or isinstance(self.right, FilmNode):
            if self.interface is None:
                raise ValueError('水膜の質点に面するつながりには、境目 interface が必要です。')
            return

        # 左右がどちらも質点のとき、同じ層なら境目なし、違う層なら境目ありでなければエラー
        same_layer = self.left.layer_index == self.right.layer_index
        if same_layer and self.interface is not None:
            raise ValueError(f'同じ層（{self.left.layer_index} 番目）の中のつながりに、境目 interface は指定できません。')
        if not same_layer and self.interface is None:
            raise ValueError(f'層 {self.left.layer_index} 番目と {self.right.layer_index} 番目の間のつながりには、境目 interface が必要です。')

    @property
    def r_h(self) -> float:
        """熱抵抗 R_h,k, (m2 K)/W（「熱・水分の流れ」の式 (6)）"""

        # 境目の熱抵抗（層の中なら 0）
        r_interface = self.interface.r_h if self.interface is not None else 0.0

        # 左から右の端まで ＋ 境目 ＋ 右から左の端まで（外気・室内・空気の質点の分は 0）
        return self.left.r_h_right + r_interface + self.right.r_h_left

    def r_m(self, state_left: State, state_right: State) -> float:
        """透湿抵抗 R_m,k,n, (m2 s Pa)/kg（「熱・水分の流れ」の式 (7)）

        左右の半分は、それぞれの側の状態（相対湿度）で計算する。外気・室内の側には、外気・室内の状態を渡す（使わないが、渡し方をそろえる）。

        Args:
            state_left: 左（質点か外気）の状態
            state_right: 右（質点か室内）の状態
        """

        # 境目の透湿抵抗（層の中なら 0）
        r_interface = self.interface.r_m if self.interface is not None else 0.0

        # 左から右の端まで ＋ 境目 ＋ 右から左の端まで（外気・室内・空気の質点の分は 0）
        return self.left.r_m_right(state_left) + r_interface + self.right.r_m_left(state_right)

    def r_l(self, state_left: State, state_right: State) -> float:
        """液水の抵抗 R_l,k,n, (m2 s (J/kg))/kg（「熱・水分の流れ」の式 (8)）

        境目（層と層、材料と通気層、表面）では液水は移動しないので無限大。層の中だけ、左右の材料の質点の液水の抵抗を足す。

        Args:
            state_left: 左の状態
            state_right: 右の状態
        """

        # 境目のあるつながりなら無限大（境目では液水は移動しない）
        if self.interface is not None:
            return float('inf')

        # 層の中（左右とも同じ層の材料の質点）
        return self.left.r_l_right(state_left) + self.right.r_l_left(state_right)


@dataclass(frozen=True, kw_only=True)
class Surface:
    """計算の側の外表面：向きと、日射・放射・風に関わる性質。

    入力の外表面（wall_input.SurfaceInput）から、WallGeometry.from_input で作る。
    値は入力で確かめ済みなので、ここでは確かめない。
    外表面の日射の熱は surface_radiation.get_solar_heat_gain が、外表面に当たる雨水は rain.get_rain_on_surface が、この値を使って求める。
    """

    # 方位（真上・真下を向いた面は TOP・BOTTOM）
    direction: Direction

    # 傾斜角, 度（0：上向きの水平、90：垂直、180：下向きの水平）
    angle: float

    # 地面の日射反射率（アルベド）, -
    albedo: float

    # 外表面の日射吸収率, -
    absorptance: float

    # 外表面の長波長放射率, -
    emissivity: float

    # 風速を求める評価高さ, m
    eva_height: float

    # 雨水暴露係数 F_E, -
    exposure_factor: float

    # 雨水付着係数 F_D, -
    deposition_factor: float



# 材料を探す関数：(出どころ, 名称) → 材料
FindMaterial = Callable[[str, str], Material]


def material_finder() -> FindMaterial:
    """Excel から材料を探す関数を作る。出どころごとに、最初に使うときに 1 回だけ読み込む。

    例:
        find = material_finder()
        find('whamcavity', '構造用合板1')
    """

    # 出どころごとの材料の一覧（読み込んだものを覚えておく）
    materials_by_source: dict[str, Materials] = {}

    def find(source: str, name: str) -> Material:

        # まだ読み込んでいない出どころなら、読み込んで覚えておく
        if source not in materials_by_source:
            materials_by_source[source] = Materials(load_materials(source))

        return materials_by_source[source].get_material(name=name)

    return find


def _make_nodes(
        layers: tuple[LayerInput, ...],
        interfaces: tuple[InterfaceInput | FilmInterfaceInput, ...],
        find_material: FindMaterial
    ) -> tuple[Node, ...]:
    """層の並び（外気側から順）から、質点の並び（外気側から順）を作る。

    材料の層（厚さ d_l、分割数 N_div,l）は、間隔 Δx_l で N_div,l 個の質点にする。層の両端の質点は表面・境目の上にあり、
    その側までの距離は 0、内側の端までの距離は Δx_l / 2。層の中の質点は、左右とも Δx_l / 2。
    通気層は、空気の質点 1 つにする。
    層と層の境目には、水膜の質点を 1 つ置く（層の番号は左の層）。

    内部用の関数。ふだんは WallGeometry.from_input から、WallInput の layers と interfaces を渡して呼ぶ。
    いちばん外側・いちばん内側が材料の層であること、通気層が 2 枚続かないこと、層と層の境目が水膜の左右 2 組であることなどは
    WallInput で確かめ済みなので、ここでは確かめない。

    Args:
        layers: 層の入力の並び（WallInput.layers）
        interfaces: 境目の入力の並び（WallInput.interfaces。層の数 + 1 個）。水膜の質点の湿気コンダクタンスに使う
        find_material: 材料を探す関数（ふだんは material_finder() で作る）

    Returns:
        質点の並び（外気側から順）

    Notes:
        計算方法の文書の「壁の構成と質点」の 3 節（式 (1)〜(3)）・4 節に対応する。
    """

    nodes: list[Node] = []

    for i, layer in enumerate(layers):

        # 層と層の境目（層 i - 1 と層 i の間）に、水膜の質点を置く（層の番号は左の層 i - 1）
        if i > 0:
            film = interfaces[i]
            nodes.append(FilmNode(layer_index=i - 1, g_m_left=film.left.g_m, g_m_right=film.right.g_m,
                                  rain_fraction=film.rain_fraction))

        # 通気層なら、空気の質点を 1 つ作る
        if isinstance(layer, AirLayerInput):
            nodes.append(AirNode(thickness=layer.thickness, alpha=layer.alpha, layer_index=i))
            continue

        # 材料の層なら、材料を探す。見つからなければ（KeyError）、どの層かを付けてエラーにする
        try:
            material = find_material(layer.source, layer.name)
        except KeyError:
            raise ValueError(f'層 {i} 番目の材料が見つかりません（出どころ {layer.source!r}、名称 {layer.name!r}）。') from None

        # 質点の間隔 Δx_l, m
        # 式 (1)。分割数は入力の確認で 2 以上なので、0 で割ることはない
        dx = layer.thickness / (layer.n_div - 1)

        # 層の左の端・中・右の端の質点を作る（端の質点は、表面・境目の側までの距離が 0）
        for j in range(layer.n_div):

            # 質点からコントロールボリュームの左の端までの距離 d_L, m
            # 式 (2)
            d_left = 0.0 if j == 0 else dx / 2

            # 質点からコントロールボリュームの右の端までの距離 d_R, m
            # 式 (3)
            d_right = 0.0 if j == layer.n_div - 1 else dx / 2

            nodes.append(MaterialNode(d_left=d_left, d_right=d_right, material=material, layer_index=i,
                                      porosity=layer.porosity, wood_decay=layer.wood_decay))

    return tuple(nodes)


def _interface(x: InterfaceInput) -> Interface:
    """入力の境目のコンダクタンス（G_h・G_m）を、計算の側の境目の抵抗（R_h,int・R_m,int）に変える（「熱・水分の流れ」の式 (1)・(2)）。"""

    # 境目の熱抵抗 R_h,int,k, (m2 K)/W と透湿抵抗 R_m,int,k, (m2 s Pa)/kg
    # 「熱・水分の流れ」の式 (1)・(2)
    return Interface(r_h=1.0 / x.g_h, r_m=1.0 / x.g_m)


def _make_links(nodes: tuple[Node, ...], interfaces: tuple[InterfaceInput | FilmInterfaceInput, ...]) -> tuple[Link, ...]:
    """質点の並びと境目の並びから、つながりの並び（外気側から順、質点の数 + 1 個）を作る。

    k 番目のつながりは、質点 k - 1 と質点 k の間（0 番目は外気と質点 0、最後は最後の質点と室内）。
    j 番目の境目は、層 j - 1 と層 j の間（0 番目は外表面、最後は内表面）。
    層と層の境目には水膜の質点があり、水膜の左のつながりに境目の左の分、右のつながりに右の分のコンダクタンスが入る。

    内部用の関数。ふだんは WallGeometry.from_input から、_make_nodes で作った質点と WallInput の interfaces を渡して呼ぶ。
    正しく作れているかは、Link を作るときの確認（境目の有無と層の関係）でも確かめられる。

    Args:
        nodes: 質点の並び（_make_nodes で作ったもの）
        interfaces: 入力の境目の並び（WallInput.interfaces。層の数 + 1 個）。計算の側の境目（Interface）に変えて使う

    Returns:
        つながりの並び（外気側から順、質点の数 + 1 個）

    Notes:
        計算方法の文書の「壁の構成と質点」の 5 節に対応する。
    """

    # 外気と最初の質点（外表面）
    links = [Link(left=OUTDOOR, right=nodes[0], interface=_interface(interfaces[0]))]

    # 隣り合う質点どうし（0 番目と 1 番目、1 番目と 2 番目、…）
    for left, right in zip(nodes[:-1], nodes[1:]):

        # 水膜の質点の左のつながり：境目（層 l と層 l + 1 の間は interfaces[l + 1]）の左の分
        if isinstance(right, FilmNode):
            interface = _interface(interfaces[right.layer_index + 1].left)

        # 水膜の質点の右のつながり：境目の右の分
        elif isinstance(left, FilmNode):
            interface = _interface(interfaces[left.layer_index + 1].right)

        # 同じ層の中のつながり：境目なし
        else:
            interface = None

        links.append(Link(left=left, right=right, interface=interface))

    # 最後の質点と室内（内表面）
    links.append(Link(left=nodes[-1], right=INDOOR, interface=_interface(interfaces[-1])))

    return tuple(links)


@dataclass(frozen=True, kw_only=True)
class WallGeometry:
    """壁の形（時間で変わらないもの）：質点・つながり・外表面・壁の大きさ。

    入口は from_input（壁の入力から作る）。壁の入力そのものは持たない（入力と計算を切り離すため）。
    質点 k は、つながり k（左）とつながり k + 1（右）の間にある。
    各層の初期の温湿度は持たない（状態の側で扱う）。
    """

    # 質点の並び（外気側から順、添え字 i）
    nodes: tuple[Node, ...]

    # つながりの並び（外気側から順、質点の数 + 1 個、添え字 k。つながり k は質点 k - 1 と質点 k の間）
    links: tuple[Link, ...]

    # 外表面
    surface: Surface

    # 縦の長さ, m
    len_vertical: float

    # 横の長さ, m
    len_horizontal: float

    def __post_init__(self):
        """作るときに、質点とつながりの並びが合っているか確かめる（from_input で作れば、いつも合っている）。

        壁の大きさは入力で確かめ済みなので、ここでは確かめない。
        """

        # つながりの数が「質点の数 + 1」でなければエラー
        if len(self.links) != len(self.nodes) + 1:
            raise ValueError(f'つながりの数は、質点の数 + 1（{len(self.nodes) + 1}）にしてください（{len(self.links)}）。')

        # つながり k の右と、つながり k + 1 の左が、どちらも質点 k でなければエラー
        for k, node in enumerate(self.nodes):
            if self.links[k].right is not node or self.links[k + 1].left is not node:
                raise ValueError(f'質点 {k} 番目の左右のつながりが、その質点につながっていません。')

    @classmethod
    def from_input(cls, wall_input: WallInput, find_material: FindMaterial | None = None) -> Self:
        """壁の入力から壁の形を作る（入力から計算への橋渡し）。

        Args:
            wall_input: 壁の入力（確かめ済み）
            find_material: 材料を探す関数。渡さなければ、Excel から探す（material_finder()）
        """

        # 材料を探す関数がなければ、Excel から探す関数を作る
        if find_material is None:
            find_material = material_finder()

        # 質点の並びと、つながりの並び
        nodes = _make_nodes(wall_input.layers, wall_input.interfaces, find_material)
        links = _make_links(nodes, wall_input.interfaces)

        # 入力の外表面を、計算の側の外表面に写す
        s = wall_input.surface
        surface = Surface(direction=s.direction, angle=s.angle, albedo=s.albedo,
                          absorptance=s.absorptance, emissivity=s.emissivity, eva_height=s.eva_height,
                          exposure_factor=s.exposure_factor, deposition_factor=s.deposition_factor)

        return cls(nodes=nodes, links=links, surface=surface,
                   len_vertical=wall_input.len_vertical, len_horizontal=wall_input.len_horizontal)

    @property
    def height(self) -> float:
        """壁の下端と上端の高さの差 H, m（通気層の浮力換気に使う）。

        通気層は、壁の下端から上端まで通っているとする。

        Notes:
            計算方法の文書の「通気層の換気」の式 (1) に対応する。
        """

        # 壁の下端と上端の高さの差 H, m
        # 式 (1)（「通気層の換気」）
        h = self.len_vertical * np.sin(np.radians(self.surface.angle))
        return h

    @property
    def area(self) -> float:
        """壁の面積, m2（縦の長さ × 横の長さ）"""
        return self.len_vertical * self.len_horizontal
