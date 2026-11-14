"""雨水の浸入：外表面に当たる雨水（風で運ばれる雨）（計算方法の文書の評価方法「雨水の浸入」）。

各値のコメントは、1 行目に「名前 記号, 単位」（無次元は -）、2 行目に「式 (k)」と補足を書く。
コメントの「式 (k)」は、計算方法の文書の「雨水の浸入」（docs/source/contents/03_14_eval_rain.md）の式の番号。
外表面に当たった雨水のうち、層と層の境目の水膜に浸入する割合（FilmNode.rain_fraction）を掛けた量が、その水膜に入る
（「熱・水分の収支」の 6 節。balance.moisture_residual_is・heat_residual_is）。
"""

import numpy as np

from modules.boundary import OutdoorBoundary
from modules.wall_geometry import Surface

# 経験的な定数 F_L, kg s/(m3 mm)（ASHRAE Standard 160 の式 (4-6)）
F_L = 0.2


def get_rain_on_surface(surface: Surface, outdoor_n: OutdoorBoundary) -> float:
    """ステップ n の、外表面に当たる雨水 r_bv,n, kg/(m2 s)。

    鉛直な壁の式（ASHRAE Standard 160 の式 (4-6)）。鉛直でない面で雨水の浸入を使わないことは、壁の入力で確かめる。

    Args:
        surface: 外表面（方位、雨水暴露係数・雨水付着係数）
        outdoor_n: ステップ n の外気の境界条件（降水量、高さ 10 m の風速の北向き・東向きの成分）

    Returns:
        ステップ n の、外表面に当たる雨水 r_bv,n, kg/(m2 s)

    Notes:
        計算方法の文書の「雨水の浸入」の式 (1)・(2) に対応する。
    """

    # 外表面の方位角 α_w, rad
    # 上向き・下向きの面は方位がないので 0 とする（雨水の浸入は鉛直な壁だけで使う）
    alpha_w = 0.0 if surface.direction.alpha is None else np.radians(surface.direction.alpha)

    # 高さ 10 m の風速のうち、外表面に向かって吹く成分 U_10,w,n, m/s（外表面から離れる向きなら 0）
    # 式 (1)
    u_10_w_n = max(outdoor_n.wind_east * np.sin(alpha_w) + outdoor_n.wind_north * np.cos(alpha_w), 0.0)

    # 外表面に当たる雨水 r_bv,n, kg/(m2 s)
    # 式 (2)。降水量 1 mm/h は 1 kg/(m2 h) なので、3600 で割って kg/(m2 s) にする
    r_bv_n = surface.exposure_factor * surface.deposition_factor * F_L * u_10_w_n * outdoor_n.rain / 3600.0

    return float(r_bv_n)
