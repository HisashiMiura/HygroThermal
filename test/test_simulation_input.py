import copy
import json
import os

import pytest

from modules.indoor_climate import CosineIndoor
from modules.outdoor_climate import RadiationInterval, RadiationKind, RadiationUnit, WeatherOutdoor
from modules.simulation_input import OutdoorInput, SimulationInput, SimulationSettings

_ROOT = os.path.join(os.path.dirname(__file__), '..')

# 計算例の壁の入力のファイル
_EXAMPLE_JSON = os.path.join(_ROOT, 'examples', 'main_wall.json')


def _example_dict() -> dict:
    """計算例の壁の入力のファイルの中身。"""
    with open(_EXAMPLE_JSON, encoding='utf-8') as f:
        return json.load(f)


def _from_dict(d: dict) -> SimulationInput:
    """計算例の壁の入力のファイルと同じフォルダを起点にして作る。"""
    return SimulationInput.from_dict(d, base_dir=os.path.dirname(_EXAMPLE_JSON))


class TestExample:
    """計算例の壁の入力のファイル（examples/main_wall.json）を読んだ中身。"""

    def test_wall(self):
        # 7 層（通気層 1 枚）、境目 8 個。構造用合板の 2 層で木材腐朽を求める
        w = SimulationInput.from_json(_EXAMPLE_JSON).wall
        names = [getattr(x, 'name', '通気層') for x in w.layers]
        assert names == ['サイディング', '通気層', '石膏ボード', '構造用合板1', 'グラスウール1', '構造用合板1', '石膏ボード']
        assert [getattr(x, 'wood_decay', False) for x in w.layers] == [False, False, False, True, False, True, False]
        assert len(w.interfaces) == 8 and (w.len_vertical, w.len_horizontal) == (7.0, 0.42)

    def test_outdoor(self):
        o = SimulationInput.from_json(_EXAMPLE_JSON).outdoor
        assert o.weather_file == os.path.abspath(os.path.join(_ROOT, 'modules', 'weather.csv'))
        assert (o.latitude, o.longitude, o.wind_height, o.wind_exponent) == (36.7, 137.21, 6.5, 0.25)
        assert o.radiation_kind is RadiationKind.INTERVAL
        assert o.radiation_interval is RadiationInterval.BEFORE
        assert o.radiation_unit is RadiationUnit.MJ_PER_M2_H
        assert o.rain_interval is RadiationInterval.BEFORE

    def test_weather_same_as_example(self):
        # 読んだ気象データの、外気の境界条件が計算例のスクリプトと同じ
        weather = _example_weather()
        loaded = SimulationInput.from_json(_EXAMPLE_JSON).outdoor.load()
        for t in (0.0, 3600.0 * 12.5, 3600.0 * 5000.0):
            assert loaded.get_outdoor_boundary(t) == weather.get_outdoor_boundary(t)

    def test_indoor(self):
        assert SimulationInput.from_json(_EXAMPLE_JSON).indoor == CosineIndoor(theta_in_ave=22.5, theta_in_amp=4.5, rh_in=60.0)

    def test_settings(self):
        s = SimulationInput.from_json(_EXAMPLE_JSON).simulation
        assert (s.dt, s.years, s.n_steps, s.record_every) == (3600.0, 2, 2 * 8760, 1)


def _example_weather() -> WeatherOutdoor:
    """付属の気象データを、入力のファイルを通さずに、同じ設定で読んだもの。"""
    return WeatherOutdoor.from_csv(
        os.path.join(_ROOT, 'modules', 'weather.csv'), latitude=36.7, longitude=137.21,
        radiation_kind=RadiationKind.INTERVAL, radiation_interval=RadiationInterval.BEFORE,
        radiation_unit=RadiationUnit.MJ_PER_M2_H, rain_interval=RadiationInterval.BEFORE, wind_height=6.5, wind_exponent=0.25)


class TestSimulationInputErrors:

    def test_missing_key(self):
        d = _example_dict()
        del d['indoor']
        with pytest.raises(ValueError, match='indoor'):
            _from_dict(d)

    def test_unknown_key(self):
        d = _example_dict()
        d['output'] = {}
        with pytest.raises(ValueError, match='output'):
            _from_dict(d)

    def test_not_dict(self):
        with pytest.raises(ValueError, match='dict'):
            _from_dict([])

    def test_part_not_dict(self):
        d = _example_dict()
        d['simulation'] = 3600
        with pytest.raises(ValueError, match='simulation'):
            _from_dict(d)

    def test_wall_error(self):
        # 壁の入力のエラーは、壁の入力が出す
        d = _example_dict()
        del d['wall']['surface']
        with pytest.raises(ValueError, match='surface'):
            _from_dict(d)

    def test_bad_json(self, tmp_path):
        p = tmp_path / 'input.json'
        p.write_text('{"wall": ', encoding='utf-8')
        with pytest.raises(ValueError, match='JSON'):
            SimulationInput.from_json(str(p))

    def test_weather_file_relative_to_json(self, tmp_path):
        # 気象データのファイルは、入力のファイルがあるフォルダから探す
        d = _example_dict()
        d['outdoor']['weather_file'] = os.path.abspath(os.path.join(_ROOT, 'modules', 'weather.csv'))
        p = tmp_path / 'input.json'
        p.write_text(json.dumps(d, ensure_ascii=False), encoding='utf-8')
        assert SimulationInput.from_json(str(p)).outdoor.weather_file == d['outdoor']['weather_file']
        d['outdoor']['weather_file'] = 'weather.csv'
        p.write_text(json.dumps(d, ensure_ascii=False), encoding='utf-8')
        with pytest.raises(ValueError, match=str(tmp_path)):
            SimulationInput.from_json(str(p))


class TestOutdoorInput:

    def _d(self) -> dict:
        return copy.deepcopy(_example_dict()['outdoor'])

    def _from(self, d: dict) -> OutdoorInput:
        return OutdoorInput.from_dict(d, base_dir=os.path.dirname(_EXAMPLE_JSON))

    def test_units(self):
        d = self._d()
        d['radiation_unit'] = 'W/m2'
        assert self._from(d).radiation_unit is RadiationUnit.W_PER_M2

    @pytest.mark.parametrize('key, value', [
        ('radiation_kind', 'average'),
        ('radiation_interval', 'after_hour'),
        ('rain_interval', 'BEFORE'),
        ('radiation_unit', 'MJ/m2h'),
    ])
    def test_unknown_word(self, key, value):
        d = self._d()
        d[key] = value
        with pytest.raises(ValueError, match=key):
            self._from(d)

    @pytest.mark.parametrize('value', ['36.7', True, None])
    def test_not_number(self, value):
        d = self._d()
        d['latitude'] = value
        with pytest.raises(ValueError, match='latitude'):
            self._from(d)

    def test_no_file(self):
        d = self._d()
        d['weather_file'] = 'no_such_file.csv'
        with pytest.raises(ValueError, match='no_such_file.csv'):
            self._from(d)

    def test_file_not_string(self):
        d = self._d()
        d['weather_file'] = 1
        with pytest.raises(ValueError, match='weather_file'):
            self._from(d)

    def test_range_checked_when_loaded(self):
        # 値の範囲は、気象データを読むときに確かめる
        d = self._d()
        d['latitude'] = 100.0
        o = self._from(d)
        with pytest.raises(ValueError, match='緯度'):
            o.load()


class TestIndoorInput:

    def _d(self) -> dict:
        return copy.deepcopy(_example_dict())

    def test_d_peak(self):
        d = self._d()
        d['indoor']['d_peak'] = 200
        assert _from_dict(d).indoor.d_peak == 200

    @pytest.mark.parametrize('value', ['sine', None])
    def test_unknown_type(self, value):
        d = self._d()
        d['indoor']['type'] = value
        with pytest.raises(ValueError, match='type'):
            _from_dict(d)

    def test_unknown_key(self):
        d = self._d()
        d['indoor']['rh'] = 50.0
        with pytest.raises(ValueError, match="'rh'"):
            _from_dict(d)

    def test_value_error(self):
        # 値のエラーは室内の条件が出し、どの入力かを頭に付ける
        d = self._d()
        d['indoor']['rh_in'] = 120.0
        with pytest.raises(ValueError, match='室内の入力.*rh_in'):
            _from_dict(d)

    def test_d_peak_float(self):
        d = self._d()
        d['indoor']['d_peak'] = 200.5
        with pytest.raises(ValueError, match='d_peak'):
            _from_dict(d)


class TestSimulationSettings:

    def test_default_record_interval(self):
        s = SimulationSettings.from_dict({'dt': 900, 'years': 1})
        assert (s.record_interval, s.n_steps, s.record_every) == (3600.0, 4 * 8760, 4)

    def test_daily(self):
        s = SimulationSettings.from_dict({'dt': 86400, 'years': 1, 'record_interval': 86400})
        assert (s.n_steps, s.record_every) == (365, 1)

    @pytest.mark.parametrize('d, match', [
        ({'dt': 0, 'years': 1}, 'dt'),
        ({'dt': -3600, 'years': 1}, 'dt'),
        ({'dt': '3600', 'years': 1}, 'dt'),
        ({'dt': 3600, 'years': 0}, 'years'),
        ({'dt': 3600, 'years': 1.5}, 'years'),
        ({'dt': 3600, 'years': True}, 'years'),
        ({'dt': 7000, 'years': 1}, '整数倍'),
        ({'dt': 900, 'years': 1, 'record_interval': 1000}, 'record_interval'),
        ({'dt': 3600, 'years': 1, 'record_interval': 0}, 'record_interval'),
        ({'dt': 3600}, 'years'),
        ({'dt': 3600, 'years': 1, 'step': 1}, 'step'),
    ])
    def test_errors(self, d, match):
        with pytest.raises(ValueError, match=match):
            SimulationSettings.from_dict(d)
