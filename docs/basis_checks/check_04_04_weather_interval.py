"""04_04「気象データと時刻」：表 4（気象データの日射の表し方の確かめ）と、6 節の表 3・表 5（区間を取り違えた影響）。

表 3：区間を後として読んだ計算と、正しく（前として）読んだ計算の差（時間刻み 3600 s）
表 5：区間を中央として読んだ計算と、前として読んだ計算の差（時間刻み 900 s。記録は 1 時間ごと）
"""

import math
import os

import numpy as np
import pandas as pd

import _common as c
from modules.outdoor_climate import RadiationInterval, WeatherOutdoor
from modules.solar_position import get_solar_position
from modules.wall_geometry import MaterialNode

_LATITUDE, _LONGITUDE = 36.7, 137.21


def _interval_patch(interval: RadiationInterval) -> tuple:
    """気象データの日射の区間を、interval として読むように差し替える。"""
    from_csv = WeatherOutdoor.from_csv.__func__

    def patched(cls, *a, **k):
        k['radiation_interval'] = interval
        return from_csv(cls, *a, **k)

    return ((WeatherOutdoor, 'from_csv', classmethod(patched)),)


def _solar(t_s: float) -> tuple[float, float]:
    """時刻 t_s（1 月 1 日 0 時からの経過時間）の太陽高度と、方位（北から時計回り）, deg。"""
    sin_h, cos_h, sin_a, cos_a = get_solar_position(t_s, _LATITUDE, _LONGITUDE)
    return math.degrees(math.asin(sin_h)), (180.0 + math.degrees(math.atan2(sin_a, cos_a))) % 360.0


def _sunrise(day: int) -> float:
    """通算日 day の日の出の時刻, 時（1 分刻みで太陽高度の正弦が正になる所を探し、前の 1 分とのあいだを直線で補う）。"""
    sin_h = lambda minute: get_solar_position((day - 1) * 86400.0 + minute * 60.0, _LATITUDE, _LONGITUDE)[0]
    for minute in range(1, 24 * 60):
        if sin_h(minute) > 0.0:
            s0, s1 = sin_h(minute - 1), sin_h(minute)
            return (minute - 1 + s0 / (s0 - s1)) / 60.0
    return float('nan')


def table_4():
    """表 4：気象データの日射の表し方の確かめ。"""
    w = pd.read_csv(os.path.join(c.ROOT, 'modules', 'weather.csv'), encoding='cp932')
    print('表 4　気象データの日射の表し方の確かめ')
    for month, day, hour in ((1, 1, 12), (6, 16, 16)):
        row = w[(w['月'] == month) & (w['日'] == day) & (w['時'] == hour)].iloc[0]
        n_day = int(pd.Timestamp(2001, month, day).dayofyear)
        t_half = (n_day - 1) * 86400.0 + (hour - 0.5) * 3600.0
        t_full = (n_day - 1) * 86400.0 + hour * 3600.0
        print(f'  {month} 月 {day} 日 {hour} 時の欄：高度 {row["太陽高度"]:.1f}°・方位 {row["太陽方位"]:.0f}°。'
              f'式では {hour - 1} 時 30 分が {_solar(t_half)[0]:.1f}°・{_solar(t_half)[1]:.0f}°、'
              f'{hour} 時が {_solar(t_full)[0]:.1f}°・{_solar(t_full)[1]:.0f}°')
    for month, day in ((1, 1), (6, 21)):
        n_day = int(pd.Timestamp(2001, month, day).dayofyear)
        sunrise = _sunrise(n_day)
        rows = w[(w['月'] == month) & (w['日'] == day)]
        first = int(rows[rows['推定日射量'] > 0.0]['時'].min())
        print(f'  {month} 月 {day} 日の日の出 {int(sunrise)} 時 {int(sunrise % 1 * 60):02d} 分ごろ、日射が 0 でなくなるのは {first} 時の欄から')
    for column in ('推定日射量', '法線面直達日射量', '水平面天空日射量'):
        print(f'  {column}で重みを付けた時刻の平均：{(w["時"] * w[column]).sum() / w[column].sum():.2f} 時')
    print(f'  南中時刻（経度による）：{12.0 - (_LONGITUDE - 135.0) / 15.0:.2f} 時')


def _diff_table(title: str, base: dict, other: dict):
    """表 3・表 5：外表面・通気層の空気・グラスウールの室内側の端の差。"""
    geometry = c.example_geometry()
    i_air = next(i for i, n in enumerate(geometry.nodes) if type(n).__name__ == 'AirNode')
    i_gw = max(i for i, n in enumerate(geometry.nodes) if isinstance(n, MaterialNode) and 'グラスウール' in n.material.name)
    y2 = c.year2(base)
    print(title)
    for i, name in ((0, '外表面'), (i_air, '通気層の空気'), (i_gw, 'グラスウールの室内側の端')):
        t_b, t_o = base['t_is_ns'][i, y2], other['t_is_ns'][i, y2]
        rh_b = c.rh_is_ns(t_b, base['mu_is_ns'][i, y2])
        rh_o = c.rh_is_ns(t_o, other['mu_is_ns'][i, y2])
        d_t, d_rh = t_o - t_b, rh_o - rh_b
        d_max = c.daily_max(t_o) - c.daily_max(t_b)
        print(f'  {name}（質点 {i}）：温度の差 {np.abs(d_t).max():.2f} K・{np.abs(d_t).mean():.2f} K、'
              f'1 日の最高温度の差 {np.abs(d_max).max():.2f} K・{d_max.mean():.2f} K、'
              f'相対湿度の差 {np.abs(d_rh).max():.1f}%・{np.abs(d_rh).mean():.1f}%')
    t0_b, t0_o = base['t_is_ns'][0, y2][:8760].reshape(365, 24), other['t_is_ns'][0, y2][:8760].reshape(365, 24)
    print(f'  外表面の 1 日の最高温度になる時刻のずれ（平均）：{(t0_o.argmax(axis=1) - t0_b.argmax(axis=1)).mean():.1f} 時間')


def main():
    table_4()
    base = c.run_example('base')
    after = c.run_example('interval_after', patches=_interval_patch(RadiationInterval.AFTER))
    _diff_table('表 3　区間を後として読んだ場合（時間刻み 3600 s、2 年目）', base, after)
    before_900 = c.run_example('interval_before_900', dt=900.0)
    centered_900 = c.run_example('interval_centered_900', dt=900.0, patches=_interval_patch(RadiationInterval.CENTERED))
    _diff_table('表 5　区間を中央として読んだ場合（時間刻み 900 s、2 年目）', before_900, centered_900)


if __name__ == '__main__':
    main()
