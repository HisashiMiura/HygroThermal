"""外表面の日射・夜間放射：壁の外表面に外から直接入る熱（計算方法の文書の評価方法「外表面の日射・夜間放射」）。

各値のコメントは、1 行目に「名前 記号, 単位」（無次元は -）、2 行目に「式 (k)」と補足を書く。
コメントの「式 (k)」は、計算方法の文書の外表面の日射・夜間放射の節（docs/source/contents/03_09_eval_surface_solar.md）の式の番号。
今の surface_solar.py・wall_surface.py（古いコード）は、段階 8 で消す。
"""

import numpy as np

from modules.boundary import SolarRadiation
from modules.wall_geometry import Surface


def get_solar_heat_gain(surface: Surface, solar: SolarRadiation) -> float:
    """外表面に外から直接入る熱 q_sol（日射の吸収 - 夜間放射）, W/m2。

    下向きの面も特別扱いしない（直達と天空は 0 になり、地面の反射だけが入る）。

    Args:
        surface: 外表面（向き・傾斜角・アルベド・日射吸収率・長波長放射率）
        solar: 太陽の位置と、日射・夜間放射

    Returns:
        外表面に外から直接入る熱 q_sol, W/m2
    """

    # 外表面の傾斜角 β_w, rad
    # 三角関数に入れるので rad に直す
    beta_w = np.radians(surface.angle)

    # 外表面の方位角 α_w, rad
    # 上向き・下向きの面は方位がないので 0 とする（sin β_w = 0 なので答えに効かない）
    alpha_w = 0.0 if surface.direction.alpha is None else np.radians(surface.direction.alpha)

    # 入射角の余弦 cos ψ, -
    # 式 (1)
    cos_psi = (np.cos(beta_w) * solar.sin_h
               + np.sin(beta_w) * solar.cos_h * (solar.cos_a * np.cos(alpha_w) + solar.sin_a * np.sin(alpha_w)))

    # 外表面に当たる直達日射量 I_D, W/m2
    # 式 (2)。太陽が地平線の下か、面の裏側にあれば 0（垂直の面では、地平線の下の太陽でも cos ψ が正になるので、太陽高度も確かめる）
    i_d = solar.i_dn * cos_psi if solar.sin_h > 0.0 and cos_psi > 0.0 else 0.0

    # 天空に対する形態係数 F_sky, -
    # 式 (3)
    f_sky = (1.0 + np.cos(beta_w)) / 2.0

    # 外表面に当たる天空日射量 I_S, W/m2
    # 式 (4)
    i_s = f_sky * solar.i_sky

    # 水平面全天日射量 I_T, W/m2
    # 式 (5)。直達の分は、太陽が地平線の上にあるときだけ（1 時間の値を一定として使うと、日の出・日の入りをまたぐ 1 時間では、
    # 太陽が地平線の下でも直達日射が残るため）
    i_t = solar.i_dn * max(solar.sin_h, 0.0) + solar.i_sky

    # 外表面に当たる地面の反射日射量 I_R, W/m2
    # 式 (6)
    i_r = surface.albedo * (1.0 - f_sky) * i_t

    # 外表面に外から直接入る熱 q_sol, W/m2
    # 式 (7)
    return float(surface.absorptance * (i_d + i_s + i_r) - surface.emissivity * f_sky * solar.r_eff)
