import math

import pytest

from modules.boundary import NO_SUN, SolarRadiation
from modules.direction import Direction
from modules.surface_radiation import get_solar_heat_gain
from modules.wall_geometry import Surface


def _surface(direction: Direction, angle: float) -> Surface:
    """アルベド 0.2、日射吸収率 0.8、長波長放射率 0.9 の外表面。"""
    return Surface(direction=direction, angle=angle, albedo=0.2, absorptance=0.8, emissivity=0.9, eva_height=4.2, exposure_factor=1.0, deposition_factor=1.0)


def _sun(altitude: float, azimuth: float, i_dn: float, i_sky: float, r_eff: float = 0.0) -> SolarRadiation:
    """太陽高度 altitude・方位角 azimuth（度、南が 0、西回りが正）の太陽と、日射・夜間放射。"""
    h, a = math.radians(altitude), math.radians(azimuth)
    return SolarRadiation(sin_h=math.sin(h), cos_h=math.cos(h), sin_a=math.sin(a), cos_a=math.cos(a),
                          i_dn=i_dn, i_sky=i_sky, r_eff=r_eff)


class TestSolarHeatGain:

    def test_no_sun(self):
        # 日射も夜間放射もなければ 0
        assert get_solar_heat_gain(_surface(Direction.S, 90.0), NO_SUN) == 0.0

    def test_roof_zenith(self):
        # 水平の屋根、太陽が真上：I = 直達 800 + 天空 100 = 900（地面は見えない）。q = 0.8 × 900 - 0.9 × 1 × 50
        q = get_solar_heat_gain(_surface(Direction.TOP, 0.0), _sun(90.0, 0.0, i_dn=800.0, i_sky=100.0, r_eff=50.0))
        assert q == pytest.approx(0.8 * 900.0 - 0.9 * 50.0)

    def test_south_wall(self):
        # 南向きの垂直の壁、太陽が真南・高度 30 度：
        # I = 直達 × cos 30° + 天空 × 0.5 +（直達 × sin 30° + 天空）× アルベド × 0.5
        q = get_solar_heat_gain(_surface(Direction.S, 90.0), _sun(30.0, 0.0, i_dn=600.0, i_sky=150.0, r_eff=40.0))
        i = 600.0 * math.cos(math.radians(30.0)) + 150.0 * 0.5 + (600.0 * 0.5 + 150.0) * 0.2 * 0.5
        assert q == pytest.approx(0.8 * i - 0.9 * 0.5 * 40.0)

    def test_sun_behind(self):
        # 東向きの垂直の壁、太陽が真西：直達は当たらず、天空と地面の反射だけ
        q = get_solar_heat_gain(_surface(Direction.E, 90.0), _sun(30.0, 90.0, i_dn=600.0, i_sky=150.0))
        i = 150.0 * 0.5 + (600.0 * 0.5 + 150.0) * 0.2 * 0.5
        assert q == pytest.approx(0.8 * i)

    def test_bottom(self):
        # 下向きの面（軒裏など）：直達と天空は 0、地面の反射だけが全部入る。夜間放射も出ていかない
        q = get_solar_heat_gain(_surface(Direction.BOTTOM, 180.0), _sun(30.0, 0.0, i_dn=600.0, i_sky=150.0, r_eff=40.0))
        assert q == pytest.approx(0.8 * (600.0 * 0.5 + 150.0) * 0.2)

    def test_sun_below_horizon(self):
        # 太陽が地平線の下（高度 -5 度）で直達日射が残っていても（1 時間の値を一定に使う日の出の前など）、
        # 直達も、地面の反射の直達の分も 0 で、地面の反射が負にならない。
        # 垂直の面では、地平線の下の太陽でも入射角の余弦が正（cos 5°）になるので、それだけでは 0 にならない
        q = get_solar_heat_gain(_surface(Direction.S, 90.0), _sun(-5.0, 0.0, i_dn=300.0, i_sky=50.0))
        assert q == pytest.approx(0.8 * (50.0 * 0.5 + 50.0 * 0.2 * 0.5))

    def test_night(self):
        # 夜（日射 0）の垂直の壁：夜間放射の半分（天空の形態係数 0.5）× 放射率だけ熱が出ていく
        q = get_solar_heat_gain(_surface(Direction.N, 90.0), _sun(-10.0, 0.0, i_dn=0.0, i_sky=0.0, r_eff=60.0))
        assert q == pytest.approx(-0.9 * 0.5 * 60.0)
