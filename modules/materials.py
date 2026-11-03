"""材料の性質：材料ごとの熱・水分の性質（計算方法の文書の評価方法「材料の性質」）。

各値のコメントは、1 行目に「名前 記号, 単位」（無次元は -）、2 行目に「式 (k)」と補足を書く。
コメントの「式 (k)」は、計算方法の文書の「材料の性質」（docs/source/contents/03_03_eval_materials.md）の式の番号。
平衡含水率・湿気伝導率・液水拡散係数の式（式 (1)〜(4)）は、material_props.py の部品に書く。
"""

from dataclasses import dataclass
from modules import thermo_dynamics
from modules.air_water_properties import get_mu
from modules.material_props import (
    FormulaIsotherm,
    PowerExpIsotherm,
    ConstantVaporConductivity,
    PowerExpVaporConductivity,
    NoLiquidDiffusivity,
    PiecewiseExpPolyDiffusivity,
)


@dataclass
class Material:
    """材料 1 つの熱・水分の性質（計算方法の文書の評価方法「材料の性質」）。

    平衡含水率・湿気伝導率・液水拡散係数の式の形は、部品（material_props.py）が持つ。
    """

    # 名前
    name: str

    # 熱伝導率 λ_h, W/(m K)
    lambda_h: float

    # 比熱 c_s, J/(kg K)
    c_s: float

    # 乾燥密度 ρ_s, kg/m3
    rho_s: float

    # 平衡含水率の部品（4 節）
    isotherm: FormulaIsotherm | PowerExpIsotherm

    # 湿気伝導率の部品（5 節）
    vapor_conductivity: ConstantVaporConductivity | PowerExpVaporConductivity

    # 液水拡散係数の部品（6 節）
    liquid_diffusivity: NoLiquidDiffusivity | PiecewiseExpPolyDiffusivity

    def get_u(self, rh: float) -> float:
        """平衡含水率 u を求める。

        Args:
            rh: 相対湿度 φ, %

        Returns:
            平衡含水率（質量含水率） u, kg/kg

        Notes:
            計算方法の文書の「材料の性質」の 4 節に対応する。
        """

        # 相対湿度の比 r, -
        # r = φ / 100（rh * 0.01 と書く。rh / 100 にすると、95 % ちょうどが 0.95 になり、境目 0.95 の式の切り替わりが今と変わる）
        r = rh * 0.01

        return self.isotherm.u(r)

    def get_lambda_dsh_m(self, rh: float) -> float:
        """湿気伝導率 λ'_m を求める。

        Args:
            rh: 相対湿度 φ, %

        Returns:
            湿気伝導率 λ'_m, kg/(m s Pa)

        Notes:
            計算方法の文書の「材料の性質」の 5 節に対応する。
        """

        # 相対湿度の比 r, -
        # r = φ / 100（rh * 0.01 と書く。rh / 100 にすると、95 % ちょうどが 0.95 になり、境目 0.95 の式の切り替わりが今と変わる）
        r = rh * 0.01

        return self.vapor_conductivity.lambda_dsh_m(r)

    def get_lambda_dsh_mu_l(self, rh: float, t: float) -> float:
        """液相水分伝導率 λ'_μl を求める。

        Args:
            rh: 相対湿度 φ, %
            t: 温度 T, K

        Returns:
            液相水分伝導率 λ'_μl, kg/(m s (J/kg))

        Notes:
            計算方法の文書の「材料の性質」の 7 節（式 (5)・(6)）に対応する。
        """

        # 液水拡散係数 D_w, m2/s
        # 6 節（相対湿度 φ の平衡含水率 u で求める）
        d_w = self.liquid_diffusivity.d_w(self.get_u(rh))

        # 液水が移動しない（拡散係数が 0）ときは 0（下の差分を計算しなくてよい）
        if d_w == 0.0:
            return 0.0

        # ずらした相対湿度 φ_+・φ_-, %
        # 式 (5)
        rh_plus, rh_minus = rh + 0.005, rh - 0.005

        # 相対湿度 φ_+・φ_- の平衡含水率 u_+・u_-, kg/kg
        # 4 節
        u_plus, u_minus = self.get_u(rh_plus), self.get_u(rh_minus)

        # 相対湿度 φ_+・φ_- の水分化学ポテンシャル μ_+・μ_-, J/kg
        # 「湿り空気と水の性質」の式 (2)
        mu_plus, mu_minus = get_mu(rh_plus, t), get_mu(rh_minus, t)

        # 液相水分伝導率 λ'_μl, kg/(m s (J/kg))
        # 式 (6)。微分 du/dμ を差分で近似する
        lambda_dsh_mu_l = self.rho_s * d_w * abs((u_plus - u_minus) / (mu_plus - mu_minus))

        return lambda_dsh_mu_l

    def get_dpsi_dmu(self, mu: float, t: float):
        """水分化学ポテンシャル変化に対する含水率変化 (DPDU)（古いコードのため。新しいコードは使わない）。

        Args:
            mu: 水分化学ポテンシャル, J/kg
            t: 絶対温度, K
        Returns:
            水分化学ポテンシャル変化に対する含水率変化, (m3/m3)/(J/kg)
        """

        return thermo_dynamics.get_dpsi_dmu(mu=mu, t=t, gma=self.rho_s, get_u=self.get_u)

    @property
    def c_rho(self):
        """乾燥した材料の容積比熱 c_s ρ_s, J/(m3 K)（3 節）"""
        return self.c_s * self.rho_s


class Materials:
    """材料の一覧。名前で材料を探す。

    材料の一覧は、呼ぶ側がローダー（material_loaders.py）で読み込んで渡す。
    例: Materials(load_whamcavity())
    """

    def __init__(self, ms: list[Material]):

        self.ms: list[Material] = ms

    def get_material(self, name: str) -> Material:
        
        ms = []

        for m in self.ms:

            if m.name == name:

                ms.append(m)
        
        if len(ms) == 0:
            raise KeyError('指定の名前の材料が見つかりませんでした。')
        
        elif len(ms) > 1:
            raise ValueError(f'同じ名前の材料が複数あります: {name}')
        
        else:

            return ms[0]
        
    def get_material_name(self) -> list[str]:
        """登録されている材料名称のリストを取得する。"""

        return [m.name for m in self.ms]
