"""外表面の日射・夜間放射：壁の外表面に外から直接入る熱（計算方法の文書の評価方法「外表面の日射・夜間放射」）。

各値のコメントは、1 行目に「名前 記号, 単位」（無次元は -）、2 行目に「式 (k)」と補足を書く。
コメントの「式 (k)」は、計算方法の文書の外表面の日射・夜間放射の節（docs/source/contents/03_09_eval_surface_solar.md）の式の番号。
変数の名前は、その節の記号に合わせる。ステップ n の値には、記号と同じく _n を付ける（例：cos ψ_n は cos_psi_n）。時刻によらない値（β_w・α_w・F_sky）には付けない。
今の surface_solar.py・wall_surface.py（古いコード）は、段階 8 で消す。
"""

import numpy as np

from modules.boundary import SolarRadiation
from modules.wall_geometry import Surface


def get_solar_heat_gain(surface: Surface, solar_n: SolarRadiation) -> float:
    """ステップ n の、外表面に外から直接入る熱 q_sol,n（日射の吸収 - 夜間放射）, W/m2。

    下向きの面も特別扱いしない（直達と天空は 0 になり、地面の反射だけが入る）。

    Args:
        surface: 外表面（方位・傾斜角・アルベド・日射吸収率・長波長放射率。時刻によらない）
        solar_n: ステップ n の太陽の位置と、日射・夜間放射

    Returns:
        ステップ n の、外表面に外から直接入る熱 q_sol,n, W/m2
    """

    # 外表面の傾斜角 β_w, rad
    # 三角関数に入れるので rad に直す
    beta_w = np.radians(surface.angle)

    # 外表面の方位角 α_w, rad
    # 上向き・下向きの面は方位がないので 0 とする（sin β_w = 0 なので答えに効かない）
    alpha_w = 0.0 if surface.direction.alpha is None else np.radians(surface.direction.alpha)

    # 入射角の余弦 cos ψ_n, -
    # 式 (1)
    cos_psi_n = (np.cos(beta_w) * solar_n.sin_h
                 + np.sin(beta_w) * solar_n.cos_h * (solar_n.cos_a * np.cos(alpha_w) + solar_n.sin_a * np.sin(alpha_w)))

    # 外表面に当たる直達日射量 I_D,n, W/m2
    # 式 (2)。太陽が地平線の下か、面の裏側にあれば 0（垂直の面では、地平線の下の太陽でも cos ψ_n が正になるので、太陽高度も確かめる）
    i_d_n = solar_n.i_dn * cos_psi_n if solar_n.sin_h > 0.0 and cos_psi_n > 0.0 else 0.0

    # 天空に対する形態係数 F_sky, -
    # 式 (3)
    f_sky = (1.0 + np.cos(beta_w)) / 2.0

    # 外表面に当たる天空日射量 I_S,n, W/m2
    # 式 (4)
    i_s_n = f_sky * solar_n.i_sky

    # 水平面全天日射量 I_T,n, W/m2
    # 式 (5)。直達の分は、太陽が地平線の上にあるときだけ（1 時間の値を一定として使うと、日の出・日の入りをまたぐ 1 時間では、
    # 太陽が地平線の下でも直達日射が残るため）
    i_t_n = solar_n.i_dn * max(solar_n.sin_h, 0.0) + solar_n.i_sky

    # 外表面に当たる地面の反射日射量 I_R,n, W/m2
    # 式 (6)
    i_r_n = surface.albedo * (1.0 - f_sky) * i_t_n

    # 外表面に外から直接入る熱 q_sol,n, W/m2
    # 式 (7)
    q_sol_n = float(surface.absorptance * (i_d_n + i_s_n + i_r_n) - surface.emissivity * f_sky * solar_n.r_eff)

    return q_sol_n
