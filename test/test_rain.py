import math
import warnings

import pytest

from modules.boundary import OutdoorBoundary
from modules.direction import Direction
from modules.rain import F_L, get_rain_on_surface
from modules.wall_geometry import Surface


def _surface(direction: Direction, f_e: float = 1.0, f_d: float = 1.0) -> Surface:
    """鉛直な壁の外表面（雨水暴露係数 f_e、雨水付着係数 f_d）。"""
    return Surface(direction=direction, angle=90.0, albedo=0.2, absorptance=0.8, emissivity=0.9, eva_height=4.2,
                   exposure_factor=f_e, deposition_factor=f_d)


def _outdoor(rain: float, north: float, east: float) -> OutdoorBoundary:
    return OutdoorBoundary.from_theta_rh(theta=10.0, rh=90.0, rain=rain, wind_north=north, wind_east=east)


class TestRainOnSurface:

    @pytest.mark.parametrize('direction, north, east', [
        (Direction.S, 3.0, 0.0),    # 南向きの壁に、南から北へ吹く風
        (Direction.W, 0.0, 3.0),    # 西向きの壁に、西から東へ吹く風
        (Direction.N, -3.0, 0.0),   # 北向きの壁に、北から南へ吹く風
        (Direction.E, 0.0, -3.0),   # 東向きの壁に、東から西へ吹く風
    ])
    def test_facing_wind(self, direction, north, east):
        # 壁に真正面から 3 m/s で吹き、降水量 10 mm/h：F_E F_D F_L U r_h / 3600
        r = get_rain_on_surface(_surface(direction), _outdoor(10.0, north, east))
        assert r == pytest.approx(1.0 * 1.0 * F_L * 3.0 * 10.0 / 3600.0)

    def test_leeward(self):
        # 壁から離れる向きの風では、雨水は当たらない
        assert get_rain_on_surface(_surface(Direction.S), _outdoor(10.0, -3.0, 0.0)) == 0.0

    def test_oblique(self):
        # 南西向きの壁に、南から北へ吹く風：壁に向かう成分は 3 cos 45 deg
        r = get_rain_on_surface(_surface(Direction.SW), _outdoor(10.0, 3.0, 0.0))
        assert r == pytest.approx(F_L * 3.0 * math.cos(math.radians(45.0)) * 10.0 / 3600.0)

    def test_factors(self):
        # 雨水暴露係数・雨水付着係数に比例する
        r1 = get_rain_on_surface(_surface(Direction.S), _outdoor(10.0, 3.0, 0.0))
        r2 = get_rain_on_surface(_surface(Direction.S, f_e=1.4, f_d=0.5), _outdoor(10.0, 3.0, 0.0))
        assert r2 == pytest.approx(1.4 * 0.5 * r1)

    def test_no_rain(self):
        assert get_rain_on_surface(_surface(Direction.S), _outdoor(0.0, 3.0, 0.0)) == 0.0
