import os

import numpy as np
import pandas as pd
import pytest

from modules.boundary import OutdoorBoundary
from modules.outdoor_climate import RadiationInterval, RadiationKind, RadiationUnit, WeatherOutdoor
from modules.thermo_dynamics import ATP, GOFF
from modules.time_axis import HOURS_PER_YEAR, SECONDS_PER_HOUR, SECONDS_PER_YEAR

_H = 3600.0

# 今の気象データ（緯度 36.7 度、経度 137.21 度、日射は前の 1 時間の積算値 MJ/(m2 h)）
_CSV = os.path.join(os.path.dirname(__file__), '..', 'modules', 'weather.csv')


def _weather(kind: RadiationKind, interval: RadiationInterval = RadiationInterval.BEFORE) -> WeatherOutdoor:
    """確かめやすい作りもののデータ：行 h の気温は 0.001 h deg.C、直達日射は h W/m2、天空日射は 2h、夜間放射は -h。"""
    h_hs = np.arange(HOURS_PER_YEAR, dtype=float)
    return WeatherOutdoor(
        latitude=36.7,
        longitude=137.21,
        theta_hs=0.001 * h_hs,
        p_v_hs=np.full(HOURS_PER_YEAR, 500.0),
        i_dn_hs=h_hs,
        i_sky_hs=2.0 * h_hs,
        r_eff_hs=-h_hs,
        radiation_kind=kind,
        radiation_interval=interval,
    )


class TestInstant:

    def test_on_the_hour(self):
        # ちょうどの時刻は、その行の値
        w = _weather(RadiationKind.INSTANT)
        assert w.get_outdoor_boundary(12 * _H).state.theta == pytest.approx(0.012)
        assert w.get_outdoor_boundary(12 * _H).solar.i_dn == pytest.approx(12.0)

    def test_between(self):
        # 11 時 30 分は、11 時と 12 時の行の平均（気温も日射も直線で補う）
        w = _weather(RadiationKind.INSTANT)
        assert w.get_outdoor_boundary(11.5 * _H).state.theta == pytest.approx(0.0115)
        assert w.get_outdoor_boundary(11.5 * _H).solar.i_dn == pytest.approx(11.5)

    def test_new_year(self):
        # 12 月 31 日 23 時 30 分は、最後の行と最初の行の平均。365 日後は同じ値（2 年目も同じデータを繰り返す）
        w = _weather(RadiationKind.INSTANT)
        last = HOURS_PER_YEAR - 1
        assert w.get_outdoor_boundary((last + 0.5) * _H).state.theta == pytest.approx(0.5 * 0.001 * last)
        assert w.get_outdoor_boundary(SECONDS_PER_YEAR + 40.25 * _H).state.theta == pytest.approx(w.get_outdoor_boundary(40.25 * _H).state.theta)


class TestInterval:

    @pytest.mark.parametrize('hour', [11.0 + 1 / 3600, 11.5, 12.0])
    def test_before(self, hour):
        # 「h 時」の欄が h - 1 〜 h 時の値：11 時 00 分 01 秒〜12 時 00 分 00 秒は、どれも 12 時の欄の値
        s = _weather(RadiationKind.INTERVAL, RadiationInterval.BEFORE).get_outdoor_boundary(hour * _H).solar
        assert (s.i_dn, s.i_sky, s.r_eff) == (12.0, 24.0, -12.0)

    def test_before_next(self):
        # 12 時 00 分 01 秒からは 13 時の欄の値
        assert _weather(RadiationKind.INTERVAL).get_outdoor_boundary(12 * _H + 1.0).solar.i_dn == 13.0

    @pytest.mark.parametrize('hour, expected', [(11.5 + 1 / 3600, 12.0), (12.0, 12.0), (12.5, 12.0), (12.5 + 1 / 3600, 13.0)])
    def test_centered(self, hour, expected):
        # 「h 時」の欄が h - 0.5 〜 h + 0.5 時の値：11 時 30 分 01 秒〜12 時 30 分 00 秒は 12 時の欄の値
        w = _weather(RadiationKind.INTERVAL, RadiationInterval.CENTERED)
        assert w.get_outdoor_boundary(hour * _H).solar.i_dn == expected

    @pytest.mark.parametrize('hour, expected', [(12.0 + 1 / 3600, 12.0), (12.5, 12.0), (13.0, 12.0), (13.0 + 1 / 3600, 13.0)])
    def test_after(self, hour, expected):
        # 「h 時」の欄が h 〜 h + 1 時の値：12 時 00 分 01 秒〜13 時 00 分 00 秒は 12 時の欄の値
        w = _weather(RadiationKind.INTERVAL, RadiationInterval.AFTER)
        assert w.get_outdoor_boundary(hour * _H).solar.i_dn == expected

    def test_temperature_is_instant(self):
        # 日射が 1 時間の値でも、気温は瞬時値として直線で補う
        assert _weather(RadiationKind.INTERVAL).get_outdoor_boundary(11.5 * _H).state.theta == pytest.approx(0.0115)

    def test_daily_total(self):
        # 1 月 1 日の直達日射を 3 分刻み（ステップの終わりの時刻の値）で足し合わせると、1〜24 時の欄の合計 × 1 時間と一致する
        w = _weather(RadiationKind.INTERVAL)
        dt = 180.0
        total = sum(w.get_outdoor_boundary((n + 1) * dt).solar.i_dn * dt for n in range(int(24 * _H / dt)))
        assert total == pytest.approx(sum(range(1, 25)) * _H)

    def test_new_year(self):
        # BEFORE では、1 月 1 日 0 時の欄は前の年の 12 月 31 日 23〜24 時の値（2 年目の 0 時ちょうども同じ欄）
        w = _weather(RadiationKind.INTERVAL)
        assert w.get_outdoor_boundary(SECONDS_PER_YEAR - 0.5 * _H).solar.i_dn == 0.0
        assert w.get_outdoor_boundary(SECONDS_PER_YEAR).solar.i_dn == 0.0


class TestAt:

    def test_boundary(self):
        # 外気の境界条件：気温と水蒸気圧が作りもののデータのとおりで、日射も持つ
        w = _weather(RadiationKind.INTERVAL)
        outdoor = w.get_outdoor_boundary(12 * _H)
        assert isinstance(outdoor, OutdoorBoundary)
        assert outdoor.state.theta == pytest.approx(0.012)
        assert outdoor.state.p_v == pytest.approx(500.0)
        assert outdoor.solar.i_dn == 12.0

    def test_sun(self):
        # 1 月 1 日の正午ごろは太陽が南の空にあり、真夜中は地平線の下
        w = _weather(RadiationKind.INTERVAL)
        noon, midnight = w.get_outdoor_boundary(12 * _H).solar, w.get_outdoor_boundary(0.0).solar
        assert noon.sin_h > 0.4
        assert abs(noon.sin_a) < 0.2 and noon.cos_a > 0.9
        assert midnight.sin_h < 0.0

    def test_rh_not_over_100(self):
        # 行の間を補った時刻でも、相対湿度は 100% を超えない（どちらの行も飽和でも、間の時刻の飽和水蒸気圧は直線より低い）
        theta_hs = np.where(np.arange(HOURS_PER_YEAR) % 2 == 0, 0.0, 20.0)
        p_v_hs = np.array([GOFF(th + ATP)[1] for th in theta_hs])
        w = WeatherOutdoor(latitude=36.7, longitude=137.21, theta_hs=theta_hs, p_v_hs=p_v_hs,
                           i_dn_hs=np.zeros(HOURS_PER_YEAR), i_sky_hs=np.zeros(HOURS_PER_YEAR),
                           r_eff_hs=np.zeros(HOURS_PER_YEAR), radiation_kind=RadiationKind.INTERVAL)
        assert w.get_outdoor_boundary(0.5 * _H).state.rh == pytest.approx(100.0)


class TestFromCsv:

    @pytest.fixture(scope='class')
    def weather(self):
        return WeatherOutdoor.from_csv(_CSV, latitude=36.7, longitude=137.21, radiation_kind=RadiationKind.INTERVAL,
                                       radiation_interval=RadiationInterval.BEFORE, radiation_unit=RadiationUnit.MJ_PER_M2_H)

    def test_values(self, weather):
        # 1 月 1 日 12 時の気温は CSV の 12 時の行、11 時 30 分の直達日射は 12 時の行（MJ/(m2 h) を W/m2 に直したもの）
        d = pd.read_csv(_CSV, encoding='shift_jis')
        assert weather.get_outdoor_boundary(12 * _H).state.theta == pytest.approx(d['気温'][12])
        assert weather.get_outdoor_boundary(11.5 * _H).solar.i_dn == pytest.approx(d['法線面直達日射量'][12] * 1e6 / 3600)

    def test_no_clamp(self, weather):
        # 今のデータには、相対湿度が 100% を超える行はない
        assert weather.n_rh_clamped == 0

    def test_readonly(self, weather):
        with pytest.raises(ValueError):
            weather.theta_hs[0] = 99.0

    def test_clamp(self, tmp_path):
        # 相対湿度が 100% を超える行は、水蒸気圧を飽和水蒸気圧に丸め、行の数を数える
        n = HOURS_PER_YEAR
        x_hs = np.full(n, 3.0)
        x_hs[[5, 100]] = 30.0  # 0 deg.C で 30 g/kg は飽和（約 3.8 g/kg）を大きく超える
        d = pd.DataFrame({'気温': np.zeros(n), '絶対湿度': x_hs, '気圧': np.full(n, 101325.0),
                          '法線面直達日射量': np.zeros(n), '水平面天空日射量': np.zeros(n), '夜間放射': np.zeros(n)})
        path = tmp_path / 'weather.csv'
        d.to_csv(path, index=False, encoding='shift_jis')
        w = WeatherOutdoor.from_csv(str(path), latitude=36.7, longitude=137.21, radiation_kind=RadiationKind.INTERVAL,
                                    radiation_interval=RadiationInterval.BEFORE, radiation_unit=RadiationUnit.MJ_PER_M2_H)
        assert w.n_rh_clamped == 2
        assert w.get_outdoor_boundary(5 * _H).state.rh == pytest.approx(100.0)


class TestInvalid:

    @pytest.mark.parametrize('kwargs, match', [
        ({'latitude': 137.21}, 'latitude'),
        ({'longitude': 200.0}, 'longitude'),
        ({'radiation_kind': 'interval'}, 'radiation_kind'),
        ({'theta_hs': np.zeros(10)}, 'theta_hs'),
        ({'i_dn_hs': np.full(HOURS_PER_YEAR, -1.0)}, '日射量'),
    ])
    def test_invalid(self, kwargs, match):
        args = {'latitude': 36.7, 'longitude': 137.21, 'theta_hs': np.zeros(HOURS_PER_YEAR),
                'p_v_hs': np.full(HOURS_PER_YEAR, 500.0), 'i_dn_hs': np.zeros(HOURS_PER_YEAR),
                'i_sky_hs': np.zeros(HOURS_PER_YEAR), 'r_eff_hs': np.zeros(HOURS_PER_YEAR),
                'radiation_kind': RadiationKind.INTERVAL} | kwargs
        with pytest.raises(ValueError, match=match):
            WeatherOutdoor(**args)
