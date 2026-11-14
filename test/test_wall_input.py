import dataclasses
import warnings

import pytest

from modules.direction import Direction
from modules.wall_input import (MaterialLayerInput, AirLayerInput, layer_input_from_dict, InterfaceInput, FilmInterfaceInput,
                                SurfaceInput, WallInput)


def _material_layer() -> MaterialLayerInput:
    """材料の層の正しい入力（構造用合板1、9 mm）。"""
    return MaterialLayerInput(source='whamcavity', name='構造用合板1', thickness=0.009, n_div=6,
                      initial_temperature=10.0, initial_rh=80.0)


class TestMaterialLayerInput:

    def test_values(self):
        layer = _material_layer()
        assert (layer.source, layer.name, layer.thickness, layer.n_div) == ('whamcavity', '構造用合板1', 0.009, 6)
        assert (layer.initial_temperature, layer.initial_rh) == (10.0, 80.0)

    def test_frozen(self):
        # 作ったあとは値を変えられない
        layer = _material_layer()
        with pytest.raises(dataclasses.FrozenInstanceError):
            layer.thickness = 0.1

    def test_keyword_only(self):
        # 名前を付けずに（位置で）値を渡すとエラー（順番の間違いを防ぐ）
        with pytest.raises(TypeError):
            MaterialLayerInput('whamcavity', '構造用合板1', 0.009, 6, 10.0, 80.0)


def _material_layer_dict() -> dict:
    """材料の層の正しい入力の dict（構造用合板1、9 mm）。"""
    return {'source': 'whamcavity', 'name': '構造用合板1', 'thickness': 0.009, 'n_div': 6,
            'initial_temperature': 10.0, 'initial_rh': 80.0}


class TestMaterialLayerInputFromDict:

    def test_read(self):
        # dict から作ったものと、直接作ったものが同じ
        assert MaterialLayerInput.from_dict(_material_layer_dict()) == _material_layer()

    @pytest.mark.parametrize('key', ['source', 'name', 'thickness', 'n_div', 'initial_temperature', 'initial_rh'])
    def test_missing_key(self, key):
        d = _material_layer_dict()
        del d[key]
        with pytest.raises(ValueError, match=key):
            MaterialLayerInput.from_dict(d)

    def test_unknown_key(self):
        d = _material_layer_dict()
        d['thicknes'] = 0.01
        with pytest.raises(ValueError, match='thicknes'):
            MaterialLayerInput.from_dict(d)


class TestMaterialLayerInputValues:

    @pytest.mark.parametrize('key, value', [
        ('source', 'excel'),
        ('name', ''),
        ('thickness', 0.0),
        ('thickness', -0.01),
        ('n_div', 1),
        ('n_div', 2.5),
        ('n_div', True),
        ('initial_rh', 0.0),
        ('initial_rh', 100.0),
    ])
    def test_invalid_value(self, key, value):
        d = _material_layer_dict()
        d[key] = value
        with pytest.raises(ValueError, match=key):
            MaterialLayerInput.from_dict(d)

    def test_n_div_2(self):
        # 分割数 2（両端だけ）は作れる
        d = _material_layer_dict()
        d['n_div'] = 2
        assert MaterialLayerInput.from_dict(d).n_div == 2


def _air_layer_dict() -> dict:
    """通気層の正しい入力の dict（25 mm、α = 0.085）。"""
    return {'air_layer': True, 'thickness': 0.025, 'alpha': 0.085, 'initial_temperature': 15.0, 'initial_rh': 50.0}


class TestAirLayerInput:

    def test_read(self):
        layer = AirLayerInput.from_dict(_air_layer_dict())
        assert (layer.thickness, layer.alpha, layer.n_div) == (0.025, 0.085, 1)
        assert (layer.initial_temperature, layer.initial_rh) == (15.0, 50.0)

    def test_frozen(self):
        layer = AirLayerInput.from_dict(_air_layer_dict())
        with pytest.raises(dataclasses.FrozenInstanceError):
            layer.alpha = 0.1

    def test_alpha_0(self):
        # α = 0（換気のない空気層）は作れる
        d = _air_layer_dict()
        d['alpha'] = 0.0
        assert AirLayerInput.from_dict(d).alpha == 0.0

    @pytest.mark.parametrize('key', ['air_layer', 'thickness', 'alpha', 'initial_temperature', 'initial_rh'])
    def test_missing_key(self, key):
        d = _air_layer_dict()
        del d[key]
        with pytest.raises(ValueError, match=key):
            AirLayerInput.from_dict(d)

    @pytest.mark.parametrize('key, value', [
        ('air_layer', False),
        ('thickness', 0.0),
        ('alpha', -0.1),
        ('initial_rh', 0.0),
        ('initial_rh', 100.0),
    ])
    def test_invalid_value(self, key, value):
        d = _air_layer_dict()
        d[key] = value
        with pytest.raises(ValueError, match=key):
            AirLayerInput.from_dict(d)

    @pytest.mark.parametrize('key', ['name', 'source', 'n_div'])
    def test_not_allowed_key(self, key):
        # 通気層には材料も分割数も指定できない
        d = _air_layer_dict()
        d[key] = 'x'
        with pytest.raises(ValueError, match=key):
            AirLayerInput.from_dict(d)


class TestLayerInputFromDict:

    def test_material_layer(self):
        # 'air_layer' がなければ材料の層
        layer = layer_input_from_dict(_material_layer_dict())
        assert isinstance(layer, MaterialLayerInput)
        assert layer == MaterialLayerInput.from_dict(_material_layer_dict())

    def test_air_layer(self):
        # 'air_layer': True なら通気層
        layer = layer_input_from_dict(_air_layer_dict())
        assert isinstance(layer, AirLayerInput)
        assert layer == AirLayerInput.from_dict(_air_layer_dict())

    def test_air_layer_false(self):
        # 'air_layer': False は、材料の層としては読まず、エラーにする
        d = _material_layer_dict()
        d['air_layer'] = False
        with pytest.raises(ValueError, match='air_layer'):
            layer_input_from_dict(d)


def _interface_dict() -> dict:
    """境目の正しい入力の dict（層と層の接触）。"""
    return {'g_h': 50.0, 'g_m': 7.75e-7}


class TestInterfaceInput:

    def test_read(self):
        interface = InterfaceInput.from_dict(_interface_dict())
        assert (interface.g_h, interface.g_m) == (50.0, 7.75e-7)

    def test_frozen(self):
        interface = InterfaceInput.from_dict(_interface_dict())
        with pytest.raises(dataclasses.FrozenInstanceError):
            interface.g_h = 9.2

    @pytest.mark.parametrize('key', ['g_h', 'g_m'])
    def test_missing_key(self, key):
        d = _interface_dict()
        del d[key]
        with pytest.raises(ValueError, match=key):
            InterfaceInput.from_dict(d)

    def test_unknown_key(self):
        d = _interface_dict()
        d['cond_hh'] = 9.2
        with pytest.raises(ValueError, match='cond_hh'):
            InterfaceInput.from_dict(d)

    @pytest.mark.parametrize('key, value', [
        ('g_h', 0.0),
        ('g_h', -1.0),
        ('g_m', 0.0),
        ('g_m', -1.0e-11),
    ])
    def test_invalid_value(self, key, value):
        d = _interface_dict()
        d[key] = value
        with pytest.raises(ValueError, match=key):
            InterfaceInput.from_dict(d)


def _surface_dict() -> dict:
    """外表面の正しい入力の dict（東向きの垂直な壁、Fortran の入力の値）。"""
    return {'direction': 'e', 'angle': 90.0, 'albedo': 0.3, 'absorptance': 0.9, 'emissivity': 0.9, 'eva_height': 4.2, 'exposure_factor': 1.0, 'deposition_factor': 1.0}


class TestSurfaceInput:

    def test_read(self):
        surface = SurfaceInput.from_dict(_surface_dict())
        assert surface.direction == Direction.E
        assert (surface.angle, surface.albedo, surface.absorptance, surface.emissivity, surface.eva_height) == (90.0, 0.3, 0.9, 0.9, 4.2)

    def test_frozen(self):
        surface = SurfaceInput.from_dict(_surface_dict())
        with pytest.raises(dataclasses.FrozenInstanceError):
            surface.angle = 45.0

    @pytest.mark.parametrize('direction, angle', [('top', 0.0), ('bottom', 180.0), ('s', 30.0)])
    def test_direction_angle(self, direction, angle):
        # 真上・真下を向いた面と、傾いた屋根（南向き 30 度）は作れる
        d = _surface_dict()
        d.update(direction=direction, angle=angle)
        assert SurfaceInput.from_dict(d).angle == angle

    @pytest.mark.parametrize('key', ['direction', 'angle', 'albedo', 'absorptance', 'emissivity', 'eva_height'])
    def test_missing_key(self, key):
        d = _surface_dict()
        del d[key]
        with pytest.raises(ValueError, match=key):
            SurfaceInput.from_dict(d)

    def test_unknown_key(self):
        # 今の InputWall の打ち間違いのキー 'absoption' はエラー
        d = _surface_dict()
        d['absoption'] = 0.9
        with pytest.raises(ValueError, match='absoption'):
            SurfaceInput.from_dict(d)

    @pytest.mark.parametrize('changes, key', [
        ({'direction': 'east'}, 'direction'),
        ({'direction': 'top', 'angle': 90.0}, 'angle'),
        ({'direction': 'bottom', 'angle': 0.0}, 'angle'),
        ({'angle': 0.0}, 'angle'),
        ({'angle': 180.0}, 'angle'),
        ({'albedo': -0.1}, 'albedo'),
        ({'absorptance': 1.1}, 'absorptance'),
        ({'emissivity': -0.1}, 'emissivity'),
        ({'eva_height': 0.0}, 'eva_height'),
    ])
    def test_invalid_value(self, changes, key):
        d = _surface_dict()
        d.update(changes)
        with pytest.raises(ValueError, match=key):
            SurfaceInput.from_dict(d)


def _wall_kwargs() -> dict:
    """壁の正しい入力の値（今の main.py・Fortran の入力の壁 7 層、境目 8 個）。"""

    def m(name, thickness, n_div, temperature, rh):
        return MaterialLayerInput(source='whamcavity', name=name, thickness=thickness, n_div=n_div,
                                  initial_temperature=temperature, initial_rh=rh)

    layers = (
        m('サイディング', 0.016, 5, 15.0, 50.0),
        AirLayerInput(thickness=0.025, alpha=0.085, initial_temperature=15.0, initial_rh=50.0),
        m('石膏ボード', 0.042, 5, 10.0, 50.0),
        m('構造用合板1', 0.009, 6, 10.0, 80.0),
        m('グラスウール1', 0.1, 5, 10.0, 50.0),
        m('構造用合板1', 0.012, 6, 10.0, 80.0),
        m('石膏ボード', 0.042, 5, 10.0, 50.0),
    )
    # 層と層の境目は水膜の左右 2 組（デフォルトの値と、防湿層だけ右を指定した値）
    def f(left, right):
        return FilmInterfaceInput(left=InterfaceInput(g_h=left[0], g_m=left[1]), right=InterfaceInput(g_h=right[0], g_m=right[1]))

    contact, material_surface, air_surface = (100.0, 1.55e-6), (1000.0, 1.0e-5), (9.2, 3.43e-8)
    interfaces = (
        InterfaceInput(g_h=22.4, g_m=2.0e-11),
        f(material_surface, air_surface),
        f(air_surface, material_surface),
        f(contact, contact),
        f(contact, contact),
        f(contact, (100.0, 2.0e-11)),
        f(contact, contact),
        InterfaceInput(g_h=9.2, g_m=2.6e-10),
    )
    return dict(surface=SurfaceInput.from_dict(_surface_dict()), layers=layers, interfaces=interfaces,
                len_vertical=7.0, len_horizontal=0.42)


class TestWallInput:

    def test_values(self):
        wall = WallInput(**_wall_kwargs())
        assert (len(wall.layers), len(wall.interfaces)) == (7, 8)
        assert isinstance(wall.layers[1], AirLayerInput)
        assert (wall.len_vertical, wall.len_horizontal) == (7.0, 0.42)

    def test_frozen(self):
        wall = WallInput(**_wall_kwargs())
        with pytest.raises(dataclasses.FrozenInstanceError):
            wall.len_vertical = 3.0

    def test_list_not_allowed(self):
        # 層・境目をリストで渡すとエラー（タプルで渡す）
        kwargs = _wall_kwargs()
        kwargs['layers'] = list(kwargs['layers'])
        with pytest.raises(ValueError, match='タプル'):
            WallInput(**kwargs)

    def test_no_layer(self):
        kwargs = _wall_kwargs()
        kwargs.update(layers=(), interfaces=(kwargs['interfaces'][0],))
        with pytest.raises(ValueError, match='1 枚も'):
            WallInput(**kwargs)

    def test_number_of_interfaces(self):
        # 境目が 1 つ足りない
        kwargs = _wall_kwargs()
        kwargs['interfaces'] = kwargs['interfaces'][:-1]
        with pytest.raises(ValueError, match='層の数 \\+ 1'):
            WallInput(**kwargs)

    @pytest.mark.parametrize('index', [0, -1])
    def test_air_layer_at_surface(self, index):
        # いちばん外側・いちばん内側に通気層を置くとエラー
        kwargs = _wall_kwargs()
        layers = list(kwargs['layers'])
        layers[index] = kwargs['layers'][1]
        kwargs['layers'] = tuple(layers)
        with pytest.raises(ValueError, match='いちばん外側'):
            WallInput(**kwargs)

    def test_air_layers_in_a_row(self):
        # 通気層が 2 枚続くとエラー（2 番目の層を通気層に替える）
        kwargs = _wall_kwargs()
        layers = list(kwargs['layers'])
        layers[2] = kwargs['layers'][1]
        kwargs['layers'] = tuple(layers)
        with pytest.raises(ValueError, match='どちらも通気層'):
            WallInput(**kwargs)

    @pytest.mark.parametrize('key, value', [
        ('len_vertical', 0.0),
        ('len_horizontal', -1.0),
    ])
    def test_invalid_value(self, key, value):
        kwargs = _wall_kwargs()
        kwargs[key] = value
        with pytest.raises(ValueError, match=key):
            WallInput(**kwargs)

    @pytest.mark.parametrize('key, value, match', [
        ('surface', {'direction': 'e'}, 'surface'),
        ('layers', ('サイディング',), 'layers'),
        ('interfaces', ({'g_h': 9.2, 'g_m': 3.43e-8},) * 8, 'interfaces'),
    ])
    def test_invalid_type(self, key, value, match):
        # 部品の入力の代わりに dict や文字列を渡すとエラー
        kwargs = _wall_kwargs()
        kwargs[key] = value
        with pytest.raises(ValueError, match=match):
            WallInput(**kwargs)


def _wall_dict() -> dict:
    """壁の正しい入力の dict（_wall_kwargs と同じ壁）。"""

    def m(name, thickness, n_div, temperature, rh):
        return {'source': 'whamcavity', 'name': name, 'thickness': thickness, 'n_div': n_div,
                'initial_temperature': temperature, 'initial_rh': rh}

    return {
        'surface': _surface_dict(),
        'layers': [
            m('サイディング', 0.016, 5, 15.0, 50.0),
            {'air_layer': True, 'thickness': 0.025, 'alpha': 0.085, 'initial_temperature': 15.0, 'initial_rh': 50.0},
            m('石膏ボード', 0.042, 5, 10.0, 50.0),
            m('構造用合板1', 0.009, 6, 10.0, 80.0),
            m('グラスウール1', 0.1, 5, 10.0, 50.0),
            m('構造用合板1', 0.012, 6, 10.0, 80.0),
            m('石膏ボード', 0.042, 5, 10.0, 50.0),
        ],
        'interfaces': [
            {'g_h': 22.4, 'g_m': 2.0e-11},
            {}, {}, {}, {},
            {'right': {'g_h': 100.0, 'g_m': 2.0e-11}},
            {},
            {'g_h': 9.2, 'g_m': 2.6e-10},
        ],
        'len_vertical': 7.0,
        'len_horizontal': 0.42,
    }


class TestWallInputFromDict:

    def test_read(self):
        # dict から作ったものと、直接作ったものが同じ
        assert WallInput.from_dict(_wall_dict()) == WallInput(**_wall_kwargs())

    def test_missing_key(self):
        d = _wall_dict()
        del d['len_horizontal']
        with pytest.raises(ValueError, match='len_horizontal'):
            WallInput.from_dict(d)

    def test_unknown_key(self):
        d = _wall_dict()
        d['blank'] = 0.0
        with pytest.raises(ValueError, match='blank'):
            WallInput.from_dict(d)

    def test_layers_not_list(self):
        # 層を dict 1 つだけで書くとエラー
        d = _wall_dict()
        d['layers'] = d['layers'][0]
        with pytest.raises(ValueError, match='layers'):
            WallInput.from_dict(d)

    def test_error_location(self):
        # 入れ子の中のエラーには、どこの入力かが付く
        d = _wall_dict()
        d['layers'][3]['thickness'] = 0.0
        with pytest.raises(ValueError, match='layers の 3 番目.*thickness'):
            WallInput.from_dict(d)

    def test_nested_not_dict(self):
        d = _wall_dict()
        d['layers'][2] = '石膏ボード'
        with pytest.raises(ValueError, match='layers の 2 番目'):
            WallInput.from_dict(d)

    def test_check_works(self):
        # 境目が 1 つ足りない dict でもエラー（WallInput の確認が dict からでも働く）
        d = _wall_dict()
        d['interfaces'] = d['interfaces'][:-1]
        with pytest.raises(ValueError, match='interfaces'):
            WallInput.from_dict(d)


class TestFilmInterface:

    def test_defaults(self):
        # 省いた側は、境目の両側の層の種類で決まるデフォルトの値（材料と通気層、通気層と材料、材料と材料）
        wall = WallInput.from_dict(_wall_dict())
        g = lambda x: (x.g_h, x.g_m)
        assert (g(wall.interfaces[1].left), g(wall.interfaces[1].right)) == ((1000.0, 1.0e-5), (9.2, 3.43e-8))
        assert (g(wall.interfaces[2].left), g(wall.interfaces[2].right)) == ((9.2, 3.43e-8), (1000.0, 1.0e-5))
        assert (g(wall.interfaces[3].left), g(wall.interfaces[3].right)) == ((100.0, 1.55e-6), (100.0, 1.55e-6))

    def test_one_side(self):
        # 片側だけ書くと、書いた側はその値、省いた側はデフォルトの値
        wall = WallInput.from_dict(_wall_dict())
        assert (wall.interfaces[5].left.g_m, wall.interfaces[5].right.g_m) == (1.55e-6, 2.0e-11)

    def test_contact_in_series(self):
        # 材料と材料のデフォルトを直列にすると、熱 50 W/(m2 K)・湿気 7.75e-7 kg/(m2 s Pa)
        wall = WallInput.from_dict(_wall_dict())
        x = wall.interfaces[3]
        assert 1.0 / (1.0 / x.left.g_h + 1.0 / x.right.g_h) == pytest.approx(50.0)
        assert 1.0 / (1.0 / x.left.g_m + 1.0 / x.right.g_m) == pytest.approx(7.75e-7)

    def test_one_set_not_allowed(self):
        # 層と層の境目を 1 組（g_h・g_m）で書くとエラー（水膜を置くので左右 2 組）
        d = _wall_dict()
        d['interfaces'][3] = {'g_h': 50.0, 'g_m': 7.75e-7}
        with pytest.raises(ValueError, match='interfaces の 3 番目.*左右 2 組'):
            WallInput.from_dict(d)

    def test_error_location(self):
        # 左右の中のエラーには、どちらの側かが付く
        d = _wall_dict()
        d['interfaces'][3] = {'left': {'g_h': 0.0, 'g_m': 1.0e-6}}
        with pytest.raises(ValueError, match='interfaces の 3 番目.*left.*g_h'):
            WallInput.from_dict(d)

    @pytest.mark.parametrize('index', [0, 3])
    def test_wrong_type(self, index):
        # 外表面に FilmInterfaceInput、層と層の境目に InterfaceInput を渡すとエラー
        kwargs = _wall_kwargs()
        interfaces = list(kwargs['interfaces'])
        interfaces[index] = kwargs['interfaces'][3] if index == 0 else kwargs['interfaces'][0]
        kwargs['interfaces'] = tuple(interfaces)
        with pytest.raises(ValueError, match=f'{index} 番目'):
            WallInput(**kwargs)


class TestRainFraction:

    def test_default_zero(self):
        # 省くと 0（浸入しない）
        wall = WallInput.from_dict(_wall_dict())
        assert all(x.rain_fraction == 0.0 for x in wall.interfaces[1:-1])

    def test_read(self):
        d = _wall_dict()
        d['interfaces'][2]['rain_fraction'] = 0.01
        assert WallInput.from_dict(d).interfaces[2].rain_fraction == 0.01

    @pytest.mark.parametrize('value', [-0.1, 1.5])
    def test_out_of_range(self, value):
        d = _wall_dict()
        d['interfaces'][2]['rain_fraction'] = value
        with pytest.raises(ValueError, match='rain_fraction'):
            WallInput.from_dict(d)

    def test_total_over_1(self):
        # 合計が 1 を超えるとエラー
        d = _wall_dict()
        d['interfaces'][1]['rain_fraction'] = 0.6
        d['interfaces'][2]['rain_fraction'] = 0.6
        with pytest.raises(ValueError, match='合計'):
            WallInput.from_dict(d)

    def test_total_over_ashrae(self):
        # 合計が 1%（ASHRAE Standard 160 のデフォルト）を超えると warning、1% ちょうどなら出さない
        d = _wall_dict()
        d['interfaces'][2]['rain_fraction'] = 0.02
        with pytest.warns(UserWarning, match='ASHRAE'):
            WallInput.from_dict(d)
        d['interfaces'][2]['rain_fraction'] = 0.01
        with warnings.catch_warnings():
            warnings.simplefilter('error')
            WallInput.from_dict(d)

    @pytest.mark.parametrize('angle', [80.0, 100.0])
    def test_not_vertical(self, angle):
        # 雨水が浸入するのに、鉛直（90 ± 5 度）でなければエラー
        d = _wall_dict()
        d['surface']['angle'] = angle
        d['interfaces'][2]['rain_fraction'] = 0.01
        with pytest.raises(ValueError, match='鉛直'):
            WallInput.from_dict(d)

    @pytest.mark.parametrize('angle', [85.0, 93.0])
    def test_nearly_vertical(self, angle):
        # 90 ± 5 度の範囲で 90 度でなければ warning
        d = _wall_dict()
        d['surface']['angle'] = angle
        d['interfaces'][2]['rain_fraction'] = 0.01
        with pytest.warns(UserWarning, match='鉛直な壁として'):
            WallInput.from_dict(d)

    def test_roof_without_rain(self):
        # 雨水が浸入しなければ、傾いた面でもよい
        d = _wall_dict()
        d['surface']['angle'] = 30.0
        WallInput.from_dict(d)

    @pytest.mark.parametrize('key', ['exposure_factor', 'deposition_factor'])
    def test_surface_factor_required(self, key):
        d = _wall_dict()
        del d['surface'][key]
        with pytest.raises(ValueError, match=key):
            WallInput.from_dict(d)

