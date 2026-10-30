import math
import os

import numpy as np
import pandas as pd
import pytest

from modules import solar_position as sp
from modules.time_axis import SECONDS_PER_DAY, SECONDS_PER_HOUR, SECONDS_PER_YEAR

_H, _D = SECONDS_PER_HOUR, SECONDS_PER_DAY

# 今の気象データ（緯度 36.7 度、経度 137.21 度）
_CSV = os.path.join(os.path.dirname(__file__), '..', 'modules', 'weather.csv')


def _altitude_deg(t: float, phi_lat: float = 36.7, lambda_lon: float = 137.21) -> float:
    """時刻 t（s）の太陽高度, deg。"""
    sin_h, _, _, _ = sp.get_solar_position(t, phi_lat=phi_lat, lambda_lon=lambda_lon)
    return math.degrees(math.asin(sin_h))


def _noon_altitude_deg(d: int, phi_lat: float) -> float:
    """経度 135 度で、通算日 d の 1 日のうちで最も高い太陽高度, deg（1 分刻みで探す）。"""
    return max(_altitude_deg((d - 1) * _D + m * 60.0, phi_lat=phi_lat, lambda_lon=135.0) for m in range(10 * 60, 14 * 60))


def _noon_time_h(d: int) -> float:
    """経度 135 度で、通算日 d に太陽が真南に来る（方位角の正弦が負から正に変わる）標準時, h（1 秒刻みで探す）。"""
    for s in range(11 * 3600, 13 * 3600):
        _, _, sin_a, _ = sp.get_solar_position((d - 1) * _D + s, phi_lat=35.0, lambda_lon=135.0)
        if sin_a >= 0.0:
            return s / 3600.0
    raise AssertionError('南中が見つからない')


class TestDayAndTime:

    def test_day_boundary(self):
        # 式 (1)・(2)：日付の変わり目の前後 1 秒で、太陽高度がほぼ続いている（通算日と標準時が正しく繰り上がる）。
        # 赤緯と均時差は日ごとの値なので、変わり目で 1 日分（0.4 度以下）だけ段がつく
        for day in (1, 31, 59, 364):
            t = day * _D
            assert _altitude_deg(t - 1.0) == pytest.approx(_altitude_deg(t), abs=0.5)

    def test_new_year(self):
        # 2 年目も同じ暦を繰り返す：365 日後の同じ時刻は、同じ太陽高度
        t = 100 * _D + 10.25 * _H
        assert _altitude_deg(t + SECONDS_PER_YEAR) == pytest.approx(_altitude_deg(t))


class TestDeclinationAndEquationOfTime:

    @pytest.mark.parametrize('d, expected_deg', [(172, 90.0 - 35.0 + 23.44), (355, 90.0 - 35.0 - 23.44), (80, 55.0)])
    def test_declination(self, d, expected_deg):
        # 式 (3)・(4)：北緯 35 度の南中の高度 = 90 - 緯度 + 赤緯。赤緯は夏至に約 +23.4 度、冬至に約 -23.4 度、春分に約 0 度
        assert _noon_altitude_deg(d, phi_lat=35.0) == pytest.approx(expected_deg, abs=0.6)

    @pytest.mark.parametrize('d, expected_min', [(307, 16.4), (42, -14.2)])
    def test_equation_of_time(self, d, expected_min):
        # 式 (5)・(6)：経度 135 度では、標準時 12 時 - 均時差 に南中する。
        # 均時差は 11 月 3 日ごろに約 +16 分（11 時 44 分ごろ南中）、2 月 11 日ごろに約 -14 分（12 時 14 分ごろ南中）
        assert (12.0 - _noon_time_h(d)) * 60.0 == pytest.approx(expected_min, abs=1.0)


class TestSolarPosition:

    def test_noon_at_135(self):
        # 経度 135 度の南中では、方位角は真南（正弦 0、余弦 1）
        t = (172 - 1) * _D + _noon_time_h(172) * _H
        _, _, sin_a, cos_a = sp.get_solar_position(t, phi_lat=35.0, lambda_lon=135.0)
        assert sin_a == pytest.approx(0.0, abs=1e-3)
        assert cos_a == pytest.approx(1.0)

    def test_unit_circle(self):
        # 正弦と余弦の 2 乗の和は 1
        for t in np.linspace(0.0, SECONDS_PER_YEAR, 37):
            sin_h, cos_h, sin_a, cos_a = sp.get_solar_position(t, phi_lat=36.7, lambda_lon=137.21)
            assert sin_h ** 2 + cos_h ** 2 == pytest.approx(1.0)
            assert sin_a ** 2 + cos_a ** 2 == pytest.approx(1.0)

    def test_morning_and_afternoon(self):
        # 午前は東の空（方位角が負）、午後は西の空（正）
        _, _, sin_a_am, _ = sp.get_solar_position(100 * _D + 9 * _H, phi_lat=36.7, lambda_lon=137.21)
        _, _, sin_a_pm, _ = sp.get_solar_position(100 * _D + 15 * _H, phi_lat=36.7, lambda_lon=137.21)
        assert sin_a_am < 0.0 < sin_a_pm

    def test_weather_csv(self):
        # 気象データ weather.csv の太陽高度の列（h 時の欄は h - 0.5 時の位置）と、2 度以内で合う
        d = pd.read_csv(_CSV, encoding='shift_jis')
        month_days = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        for row in (10, 12, 4000, 4005, 8000):
            r = d.iloc[row]
            day = sum(month_days[:int(r['月']) - 1]) + int(r['日'])
            t = (day - 1) * _D + (r['時'] - 0.5) * _H
            sin_h, _, _, _ = sp.get_solar_position(t, phi_lat=36.7, lambda_lon=137.21)
            assert math.degrees(math.asin(sin_h)) == pytest.approx(r['太陽高度'], abs=2.0)
