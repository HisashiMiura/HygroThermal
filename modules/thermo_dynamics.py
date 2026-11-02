"""古い熱・水分の物性の関数（段階 8 で消す）。

湿り空気と水の定数と関数は air_water_properties.py に移した。ここでは、古いコードのためにそれを使う。
含水率の微分（get_dpsi_dmu）と液相水分伝導率（get_lambda_dsh_mu_l）は、計算方法の文書の「材料の性質」を書くときに移す。
"""

import math
from typing import Callable

from modules.air_water_properties import C_W, L_V, R_V, RHO_W, T_0, get_mu, get_p_vs, get_rh, get_rho_a

# 古いコードのための、Fortran から引き継いだ名前
ATP = T_0
RG = R_V
RW = L_V
CPL = C_W
ROW = RHO_W
ROW_CP = 1300.0  # 空気の容積比熱, J/(m3 K)（新しいコードは、空気の比熱 C_A と密度 get_rho_a に分けた）

RMVA = 0.6217
"""乾き空気の平均分子量(約28.966)に対する水蒸気の分子量(約18.015)の比"""

P_ATM = 101325.0
"""標準大気圧, Pa"""

COND_H_O = 22.4
"""室外側表面熱伝達率, W/(m2 K)"""

COND_H_I = 9.2
"""室内側表面熱伝達率, W/(m2 K)"""

COND_H_AIR = 9.2
"""通気層表面熱伝達率, W/(m2 K)"""

COND_M_O = 2.0e-11
"""室外側表面湿気伝達率, (kg/s)/(m2 Pa)"""

COND_M_I = 3.43e-08
"""室内側表面湿気伝達率, (kg/s)/(m2 Pa)"""

COND_M_AIR = 3.43e-08
"""通気層表面湿気伝達率, (kg/s)/(m2 Pa)"""

R_M_BS = 2.4e5
"""バックシーラー透水抵抗(=2.4e+5 m2sPa/kg by 長村)"""

def get_x(p_v: float) -> float:
    """絶対湿度を求める。

    Args:
        p_v: 水蒸気圧, Pa

    Returns:
        絶対湿度, kg/kg(DA)
    """

    return 0.622 * p_v / (P_ATM - p_v)

def get_rho(t: float) -> float:
    """空気の密度, kg/m3（古いコードのため。air_water_properties.get_rho_a を使う）。"""

    return get_rho_a(t)


def get_wp(rh: float, t: float) -> float:
    """水分化学ポテンシャル, J/kg（古いコードのため。air_water_properties.get_mu を使う）。"""

    return get_mu(rh, t)


def GOFF(t: float) -> tuple[float, float]:
    """飽和水蒸気圧（古いコードのため。air_water_properties.get_p_vs を使う）。

    Args:
        t: 絶対温度, K

    Returns:
        FS: mmHg
        VP: Pa
    """

    vp = get_p_vs(t)

    FS = vp / 133.322

    return FS, vp


def FUNCX(rh: float, vp: float) -> float:
    """混合比（乾き空気 1 kg あたりの水蒸気重量 [kg]）を求める。

    Args:
        rh: 相対湿度, %
        vp: 飽和水蒸気圧, Pa

    Returns:
        混合比, kg/kg(DA)
    """

    # 133.322 Pa / mmHg
    # 大気圧 101324.72 Pa = 760.0 * 133.322

    return vp * rh * 0.01 * RMVA / (760.0 * 133.322 - vp * rh * 0.01)
    

def get_dpv_dmu(rh: float, t: float) -> float:
    """含水率変化に対する絶対湿度の変化率(DGDU)を計算する。

    Args:
        rh: 相対湿度, %
        tmp: 絶対温度, K
    
    Returns:
        水分化学ポテンシャルに対する水蒸気圧の微分, Pa / (J / kg)
    """

    # --- DGDU の計算 (湿度差による微分) ---
    rh1 = rh + 0.01
    rh2 = rh - 0.01
    
    # 外部定義されている前提の関数
    # get_wp(相対湿度, 温度, パラメータ)
    wp1 = get_wp(rh1, t)
    wp2 = get_wp(rh2, t)
    
    # GOFF関数(飽和水蒸気圧等)から絶対湿度を計算
    # 以前の回答通り、2つの戻り値を想定
    _, vp = GOFF(t)
    x1 = vp * rh1 * 0.01
    x2 = vp * rh2 * 0.01
    
    # ゼロ除算を回避しつつ微分布を計算
    if abs(wp1 - wp2) > 1e-15:
        dgdu = abs((x1 - x2) / (wp1 - wp2))
    else:
        dgdu = 0.0

    return dgdu


def get_dpv_dt(rh: float, t: float) -> float:
    """温度変化に対する絶対湿度の変化率(DGDT)を計算する。

    Args:
        rh: 相対湿度, %
        tmp: 絶対温度, K
    
    Returns:
        絶対温度に対する水蒸気圧の微分, Pa / K

    """

    # --- DGDT の計算 (温度差による微分) ---
    tp1 = t + 0.1
    tp2 = t - 0.1
    
    # 温度を変えて飽和水蒸気圧を取得
    _, vp1 = GOFF(tp1)
    _, vp2 = GOFF(tp2)
    
    x1_t = vp1 * rh * 0.01
    x2_t = vp2 * rh * 0.01
    
    # 温度勾配による絶対湿度の変化率
    # tp1 - tp2 は常に 0.2 なのでゼロ除算の心配はほぼありません
    dgdt = abs((x1_t - x2_t) / (tp1 - tp2))

    return dgdt


def get_lambda_dsh_mu_l(
        rh: float,
        t: float,
        rho_s: float,
        get_u: Callable[[float], float],
        d_w: Callable[[float], float]
    ) -> float:
    """水分化学ポテンシャル勾配に対する液相水分伝導率 (DIFF) を計算する。

    Fortran（Whamcavity.f95）の SUBROUTINE DIFF に当たる。
    Fortran では構造用合板の係数と範囲が直書きされていたが、ここでは液水拡散係数を求める関数 d_w を受け取る。

    液相水分伝導率 = 材料の乾燥密度 × 液水拡散係数 × |Δu / Δμ|

    液水拡散係数は、出典（RR-110 の式 (5)、ASHRAE RP-1018）で J_l = -ρ0 D_w grad u（ρ0：乾燥密度、u：質量含水率 kg/kg）と
    定義されているので、掛けるのは材料の乾燥密度（docs/governing_equations.md の式 3.4）。
    Fortran の DIFF は水の密度 998 を掛けているが、意図して変えた（2026-10-24、docs/note_fortran_diff_density.md）。

    Args:
        rh: 相対湿度(0.0～100.0), %
        t: 絶対温度, K
        rho_s: 材料の乾燥密度, kg/m3
        get_u: 相対湿度(%)から質量含水率(kg/kg)を求める関数
        d_w: 質量含水率(kg/kg)から液水拡散係数(m2/s)を求める関数
    Returns:
        液相水分伝導率, (kg/s)/(m (J/kg))
    """

    # 液水拡散係数, m2/s
    dw = d_w(get_u(rh))

    # 液水が移動しない（拡散係数が 0）ときは 0
    if dw == 0.0:
        return 0.0

    # 含水率の水分化学ポテンシャルによる微分を、相対湿度 ±0.005% の差分で近似する
    rh1 = rh + 0.005
    rh2 = rh - 0.005

    # 質量含水率, kg/kg
    u1 = get_u(rh1)
    u2 = get_u(rh2)

    # 水分化学ポテンシャル, J/kg
    wp1 = get_wp(rh1, t)
    wp2 = get_wp(rh2, t)

    # ゼロ除算のチェック
    if abs(wp1 - wp2) <= 1e-12:
        return 0.0

    # 拡散係数(m2/s) は質量含水率の勾配(kg/kg)を駆動力とする係数なので、
    # 乾燥密度を掛けて含水量(kg/m3)の勾配の基準に、Δu/Δμ を掛けて水分化学ポテンシャル差(J/kg)基準に直す。
    return rho_s * dw * abs((u1 - u2) / (wp1 - wp2))


def get_dpsi_dmu(mu: float, t: float, gma: float, get_u: callable) -> float:
    """水分化学ポテンシャル変化に対する含水率変化 (DPDU) を計算する。
        (m3/m3)/(J/kg)
        分子：含水率, m3/m3
        分母：水分化学ポテンシャル, J/kg
    
    Args:
        mu: 水分化学ポテンシャル, J/kg
        t: 絶対温度, K
        gma: 材料密度, kg/m3
        get_u: 含水率を計算する関数
    """

    # 差分幅の計算 (d1の1%), J/kg
    dmu = abs(mu * 0.01)
    
    # 微小変化させた含水率, J/kg
    mu1 = mu + dmu
    mu2 = mu - dmu
    
    # 相対湿度, %
    rh1 = get_rh(mu1, t)
    rh2 = get_rh(mu2, t)
    
    # 質量基準の含水率外部定義されている前提の AHGANS, 含水率
    u1 = get_u(rh1)
    u2 = get_u(rh2)
    
    # 体積基準の含水率, m3/m3
    # gma: 材料密度　kg/m3
    # ROW: 水の密度　kg/m3
    psi1 = u1 * gma / ROW
    psi2 = u2 * gma / ROW
    
    # 中央差分による勾配(微分値)の近似
    # 0.5 * (VGT1 - VGT2) / DW
    # dw を2回たしているので2でわっている。
    dpsi_dmu = (psi1 - psi2) / (2 * dmu) 
        
    return dpsi_dmu
