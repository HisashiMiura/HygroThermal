"""04_06「境目にたまる水（水膜）のモデル」10.3 節：7 層の壁で、雨水を浸入させた 2 年間の計算。

通気層と石膏ボードの境目（境目 2）に、外表面に当たる雨水の 1% を浸入させる。
水膜の右（石膏ボードの側）の湿気コンダクタンスを、デフォルト（1.0e-5）、2e-9、2e-10 の 3 通りにする。
比べる相手は、同じ湿気コンダクタンスで、雨水を浸入させない計算。
"""

import dataclasses

import numpy as np

import _common as c
from modules.rain import get_rain_on_surface
from modules.state import State

# 雨水を浸入させる境目の番号（通気層と石膏ボードの境目）
_INTERFACE = 2

# 水膜の左（通気層の空気の側）のつながりの湿気コンダクタンス, kg/(m2 s Pa)
_G_M_AIR = 3.43e-8


def _wall_patch(g_m_right: float | None, rain_fraction: float) -> tuple:
    """境目 2 の水膜の右の湿気コンダクタンスと、浸入する雨水の割合を変える。"""
    wall_input = c.ex._main_wall_input

    def patched():
        w = wall_input()
        interfaces = list(w.interfaces)
        f = interfaces[_INTERFACE]
        right = f.right if g_m_right is None else dataclasses.replace(f.right, g_m=g_m_right)
        interfaces[_INTERFACE] = dataclasses.replace(f, right=right, rain_fraction=rain_fraction)
        return dataclasses.replace(w, interfaces=tuple(interfaces))

    return ((c.ex, '_main_wall_input', patched),)


def _weather():
    """例の壁の気象データ。"""
    from modules.outdoor_climate import RadiationInterval, RadiationKind, RadiationUnit, WeatherOutdoor
    return WeatherOutdoor.from_csv(
        c.os.path.join(c.ROOT, 'modules', 'weather.csv'), latitude=36.7, longitude=137.21,
        radiation_kind=RadiationKind.INTERVAL, radiation_interval=RadiationInterval.BEFORE,
        radiation_unit=RadiationUnit.MJ_PER_M2_H, rain_interval=RadiationInterval.BEFORE, wind_height=6.5, wind_exponent=0.25)


def main():
    geometry = c.example_geometry()
    i_air = next(i for i, n in enumerate(geometry.nodes) if type(n).__name__ == 'AirNode')
    i_film, i_gyp = i_air + 1, i_air + 2

    # 外表面に当たる雨水（2 年目、1 時間ごと）
    weather = _weather()
    base = c.run_example('base')
    y2 = c.year2(base)
    time_y2 = base['time_ns'][y2]
    r_bv = np.array([get_rain_on_surface(geometry.surface, weather.get_outdoor_boundary(t)) for t in time_y2])
    rain_hours = r_bv > 0.0
    print(f'外表面に当たった雨水 {r_bv.sum() * 3600:.0f} kg/m2（水平面の降水量 {weather.rain_hs.sum():.0f} mm）、'
          f'雨水が当たった時間 {rain_hours.sum()} 時間、水膜に浸入した雨水 {0.01 * r_bv.sum() * 3600:.1f} kg/m2')

    for label, g_m in (('デフォルト（1.0e-5）', None), ('2e-9', 2.0e-9), ('2e-10', 2.0e-10)):
        tag = 'default' if g_m is None else f'{g_m:.0e}'
        dry = base if g_m is None else c.run_example(f'film_dry_{tag}', patches=_wall_patch(g_m, 0.0))
        wet = c.run_example(f'film_rain_{tag}', patches=_wall_patch(g_m, 0.01))
        s = wet['mu_is_ns'][i_film, y2]
        k = geometry.nodes[i_film].k if g_m is None else dataclasses.replace(geometry.nodes[i_film], g_m_right=g_m).k
        w_film = k * np.maximum(s, 0.0)
        print(f'[{label}] 水膜があった時間 {(s > 0).sum()} h、水膜の水の最大 {w_film.max() * 1000:.1f} g/m2')
        for i, name in ((i_gyp, '石膏ボードの通気層側の端'), (i_air, '通気層の空気')):
            d_rh = c.rh_is_ns(wet['t_is_ns'][i, y2], wet['mu_is_ns'][i, y2]) - c.rh_is_ns(dry['t_is_ns'][i, y2], dry['mu_is_ns'][i, y2])
            print(f'  {name}（質点 {i}）：相対湿度の上がり 最大 {d_rh.max():.0f}%、平均 {d_rh.mean():.1f}%')
        print(f'  温度の差の最大 {np.abs(wet["t_is_ns"][:, y2] - dry["t_is_ns"][:, y2]).max():.1f} K、'
              f'修正の回数 平均 {dry["iterations_ns"].mean():.2f} 回 → {wet["iterations_ns"].mean():.2f} 回、最大 {wet["iterations_ns"].max()} 回')

    # 水膜が全面を覆うとした、通気層への蒸発の勢い（雨水を浸入させない計算の通気層の空気で求める）
    t_film, t_air, mu_air = base['t_is_ns'][i_film, y2], base['t_is_ns'][i_air, y2], base['mu_is_ns'][i_air, y2]
    cap = np.array([_G_M_AIR * (State(t=a, mu=0.0).p_v - State(t=b, mu=min(m, 0.0)).p_v) for a, b, m in zip(t_film, t_air, mu_air)])
    j_rain = 0.01 * r_bv
    print(f'雨の時間の、通気層への蒸発の勢い：中央値 {np.median(cap[rain_hours]) * 3.6e6:.0f} g/(m2 h)、'
          f'浸入する雨水：中央値 {np.median(j_rain[rain_hours]) * 3.6e6:.2f} g/(m2 h)、最大 {j_rain.max() * 3.6e6:.0f} g/(m2 h)')


if __name__ == '__main__':
    main()
