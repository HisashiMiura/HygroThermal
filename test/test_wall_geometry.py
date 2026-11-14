import dataclasses
import math

import pytest

from modules.direction import Direction
from modules.material_loaders import load_materials
from modules.materials import Materials
from modules.state import State
from modules.air_water_properties import C_A, C_W, R_V, T_0, get_rho_a
from modules import wall_geometry
from modules.wall_geometry import MaterialNode, AirNode, FilmNode, _make_nodes, _make_links, material_finder, Link, Interface, Surface, WallGeometry, Outdoor, Indoor, OUTDOOR, INDOOR
from modules.wall_input import MaterialLayerInput, AirLayerInput, InterfaceInput, FilmInterfaceInput, WallInput


def _plywood():
    """構造用合板1（旧 Excel）。"""
    return Materials(load_materials('whamcavity')).get_material(name='構造用合板1')


class TestMaterialNode:

    @pytest.mark.parametrize('d_left, d_right, width', [
        (0.001, 0.001, 0.002),   # 層の中
        (0.0, 0.001, 0.001),     # 層の左の端
        (0.001, 0.0, 0.001),     # 層の右の端
    ])
    def test_width(self, d_left, d_right, width):
        node = MaterialNode(d_left=d_left, d_right=d_right, material=_plywood(), layer_index=0)
        assert node.width == pytest.approx(width)

    def test_heat_capacity(self):
        # 幅 2 mm × 容積比熱
        m = _plywood()
        node = MaterialNode(d_left=0.001, d_right=0.001, material=m, layer_index=0)
        assert node.heat_capacity_dry == pytest.approx(0.002 * m.c_s * m.rho_s)

    @pytest.mark.parametrize('d_left, d_right', [
        (0.001, 0.001),   # 層の中
        (0.0, 0.001),     # 層の左の端（左は 0）
        (0.001, 0.0),     # 層の右の端（右は 0）
    ])
    def test_r_h(self, d_left, d_right):
        # 端までの熱抵抗 = 距離 / 熱伝導率
        m = _plywood()
        node = MaterialNode(d_left=d_left, d_right=d_right, material=m, layer_index=0)
        assert node.r_h_left == pytest.approx(d_left / m.lambda_h)
        assert node.r_h_right == pytest.approx(d_right / m.lambda_h)

    def test_frozen(self):
        node = MaterialNode(d_left=0.001, d_right=0.001, material=_plywood(), layer_index=0)
        with pytest.raises(dataclasses.FrozenInstanceError):
            node.d_left = 0.0

    @pytest.mark.parametrize('d_left, d_right', [(-0.001, 0.001), (0.001, -0.001), (0.0, 0.0)])
    def test_invalid_distance(self, d_left, d_right):
        with pytest.raises(ValueError, match='d_left'):
            MaterialNode(d_left=d_left, d_right=d_right, material=_plywood(), layer_index=0)


class TestAirNode:

    def test_values(self):
        # 厚さ 25 mm、α = 0.085
        node = AirNode(thickness=0.025, alpha=0.085, layer_index=1)
        assert (node.width, node.alpha) == (0.025, 0.085)

    def test_r_h(self):
        # よく混ざった空気なので、端までの熱抵抗は 0
        node = AirNode(thickness=0.025, alpha=0.085, layer_index=1)
        assert (node.r_h_left, node.r_h_right) == (0.0, 0.0)

    @pytest.mark.parametrize('thickness, alpha, key', [(0.0, 0.085, 'thickness'), (0.025, -0.1, 'alpha')])
    def test_invalid_value(self, thickness, alpha, key):
        with pytest.raises(ValueError, match=key):
            AirNode(thickness=thickness, alpha=alpha, layer_index=1)


class TestFilmNode:

    def _film(self) -> FilmNode:
        # 通気層の表面のデフォルト（材料の側 1.0e-5、空気の側 3.43e-8）
        return FilmNode(layer_index=0, g_m_left=1.0e-5, g_m_right=3.43e-8)

    def test_no_thickness(self):
        # 厚さがないので、幅・熱容量・端までの抵抗は 0
        f = self._film()
        s = State(t=293.15, mu=-1000.0)
        assert (f.width, f.heat_capacity(s), f.r_h_left, f.r_h_right, f.r_m_left(s), f.r_m_right(s)) == (0.0,) * 6

    def test_k(self):
        # k = τ (G_m,L + G_m,R) p_vs(T_ref) / (R_v T_ref)、τ = 3600 s、T_ref = 293.15 K
        from modules.air_water_properties import get_p_vs
        f = self._film()
        assert f.k == pytest.approx(3600.0 * (1.0e-5 + 3.43e-8) * get_p_vs(293.15) / (R_V * 293.15))

    @pytest.mark.parametrize('s, mu, w_per_k', [(-500.0, -500.0, 0.0), (0.0, 0.0, 0.0), (800.0, 0.0, 800.0)])
    def test_state(self, s, mu, w_per_k):
        # 水膜の状態量 s < 0：水膜なし（μ = s、水分量 0）。s ≥ 0：水膜あり（μ = 0、水分量 k s）
        f = self._film()
        assert FilmNode.mu(s) == mu
        assert f.water_amount(State(t=293.15, mu=s)) == pytest.approx(f.k * w_per_k)

    def test_link_needs_interface(self):
        # 水膜の質点に面するつながりには、境目が必要
        a = _node('構造用合板1', 0.0009, 0.0, 0)
        with pytest.raises(ValueError, match='水膜'):
            Link(left=a, right=self._film(), interface=None)

    def test_link_same_layer_index(self):
        # 水膜の質点の層の番号は左の層と同じだが、境目のあるつながりにできる
        a = _node('構造用合板1', 0.0009, 0.0, 0)
        link = Link(left=a, right=self._film(), interface=Interface(r_h=1 / 1000.0, r_m=1 / 1.0e-5))
        assert link.r_h == pytest.approx(1 / 1000.0)


def _layer(name: str, thickness: float, n_div: int) -> MaterialLayerInput:
    """旧 Excel の材料の層の入力（初期の温湿度は 10 deg.C・50%）。"""
    return MaterialLayerInput(source='whamcavity', name=name, thickness=thickness, n_div=n_div,
                              initial_temperature=10.0, initial_rh=50.0)


def _air_layer() -> AirLayerInput:
    """通気層 25 mm、α = 0.085。"""
    return AirLayerInput(thickness=0.025, alpha=0.085, initial_temperature=15.0, initial_rh=50.0)


def _main_wall_layers() -> tuple:
    """Fortran の入力の壁 7 層。"""
    return (
        _layer('サイディング', 0.016, 5), _air_layer(), _layer('石膏ボード', 0.042, 5), _layer('構造用合板1', 0.009, 6),
        _layer('グラスウール1', 0.1, 5), _layer('構造用合板1', 0.012, 6), _layer('石膏ボード', 0.042, 5),
    )


def _surf(g_h: float) -> InterfaceInput:
    """外表面・内表面の入力の境目（熱コンダクタンス g_h、湿気は 7.75e-7）。"""
    return InterfaceInput(g_h=g_h, g_m=7.75e-7)


def _film(g_h: float, g_m: float = 7.75e-7) -> FilmInterfaceInput:
    """層と層の入力の境目。水膜の左右に 2 倍ずつ割り振り、直列にすると g_h・g_m になる。"""
    half = InterfaceInput(g_h=2.0 * g_h, g_m=2.0 * g_m)
    return FilmInterfaceInput(left=half, right=half)


def _one_layer_interfaces() -> tuple:
    return (_surf(22.4), _surf(9.2))


class TestMakeNodes:

    def test_one_layer(self):
        # 合板 9 mm・分割数 6：質点 6 個、両端は表面の側までの距離が 0
        nodes = _make_nodes((_layer('構造用合板1', 0.009, 6),), _one_layer_interfaces(), find_material=material_finder())
        assert len(nodes) == 6
        assert nodes[0].d_left == 0.0
        assert nodes[-1].d_right == 0.0
        assert sum(n.width for n in nodes) == pytest.approx(0.009)

    def test_interior_width(self):
        # 層の中の質点の幅は dx = 9 mm / 5 = 1.8 mm、両端は dx / 2
        nodes = _make_nodes((_layer('構造用合板1', 0.009, 6),), _one_layer_interfaces(), find_material=material_finder())
        assert [n.width for n in nodes] == pytest.approx([0.0009, 0.0018, 0.0018, 0.0018, 0.0018, 0.0009])

    def test_two_layers_boundary(self):
        # 合板とグラスウールの境目に、同じ位置の質点が 3 つ並ぶ（合板の端、水膜、グラスウールの端。どれも境目の側までの距離が 0）
        layers = (_layer('構造用合板1', 0.009, 6), _layer('グラスウール1', 0.1, 5))
        nodes = _make_nodes(layers, (_surf(22.4), _film(50.0, 7.75e-7), _surf(9.2)), find_material=material_finder())
        left, film, right = nodes[5], nodes[6], nodes[7]
        assert (left.d_right, film.width, right.d_left) == (0.0, 0.0, 0.0)
        assert isinstance(film, FilmNode)
        assert (left.layer_index, film.layer_index, right.layer_index) == (0, 0, 1)
        assert (film.g_m_left, film.g_m_right) == (1.55e-6, 1.55e-6)
        assert (left.material.name, right.material.name) == ('構造用合板1', 'グラスウール1')

    def test_air_layer(self):
        # 通気層は空気の質点 1 つで、流量係数は入力の α のまま
        layers = (_layer('サイディング', 0.016, 5), _air_layer(), _layer('石膏ボード', 0.042, 5))
        nodes = _make_nodes(layers, (_surf(22.4), _film(9.2), _film(9.2), _surf(9.2)), find_material=material_finder())
        air = nodes[6]
        assert isinstance(air, AirNode)
        assert (air.thickness, air.alpha) == (0.025, 0.085)
        assert air.layer_index == 1

    def test_main_wall(self):
        # 今の壁 7 層：質点 39 個（材料と通気層の 5 + 1 + 5 + 6 + 5 + 6 + 5 = 33 個と、層と層の境目の水膜 6 個）
        nodes = _make_nodes(_main_wall_layers(), _main_wall_interfaces(), find_material=material_finder())
        assert len(nodes) == 39
        assert [i for i, n in enumerate(nodes) if isinstance(n, FilmNode)] == [5, 7, 13, 20, 26, 33]

    def test_width_of_each_layer(self):
        # 各層の質点の幅の合計が、その層の厚さ
        layers = _main_wall_layers()
        nodes = _make_nodes(layers, _main_wall_interfaces(), find_material=material_finder())
        for i, layer in enumerate(layers):
            assert sum(n.width for n in nodes if n.layer_index == i) == pytest.approx(layer.thickness)

    def test_material_not_found(self):
        layers = (_layer('構造用合板1', 0.009, 6), _layer('構造用合板3', 0.009, 6))
        with pytest.raises(ValueError, match="層 1 番目.*構造用合板3"):
            _make_nodes(layers, (_surf(22.4), _film(50.0), _surf(9.2)), find_material=material_finder())


class TestMaterialFinder:

    def test_load_once(self, monkeypatch):
        # 同じ出どころを何回探しても、Excel を読み込むのは 1 回だけ
        calls = []

        def counting_load_materials(source):
            calls.append(source)
            return load_materials(source)

        monkeypatch.setattr(wall_geometry, 'load_materials', counting_load_materials)
        find = material_finder()
        find('whamcavity', '構造用合板1')
        find('whamcavity', 'グラスウール1')
        assert calls == ['whamcavity']


def _node(name: str, d_left: float, d_right: float, layer_index: int) -> MaterialNode:
    """旧 Excel の材料の質点。"""
    return MaterialNode(d_left=d_left, d_right=d_right, material=material_finder()('whamcavity', name), layer_index=layer_index)


def _interface(g_h: float) -> Interface:
    """熱コンダクタンス g_h の境目（湿気は 7.75e-7）。計算の側の境目なので抵抗で持つ。"""
    return Interface(r_h=1.0 / g_h, r_m=1.0 / 7.75e-7)


class TestLinkHeat:

    def test_interior(self):
        # 合板の中の 2 つの質点（左右とも 0.9 mm）：熱抵抗 = 1.8 mm / λ
        a, b = _node('構造用合板1', 0.0009, 0.0009, 0), _node('構造用合板1', 0.0009, 0.0009, 0)
        link = Link(left=a, right=b, interface=None)
        assert link.r_h == pytest.approx(0.0018 / a.material.lambda_h)

    def test_layer_to_layer(self):
        # 合板の右の端と GW の左の端：境目の 1/50 だけ（1 回だけ入る）
        a, b = _node('構造用合板1', 0.0009, 0.0, 0), _node('グラスウール1', 0.0, 0.0125, 1)
        link = Link(left=a, right=b, interface=_interface(50.0))
        assert link.r_h == pytest.approx(1 / 50.0)

    def test_outside_surface(self):
        # 外気とサイディングの左の端：外表面の伝達 1/22.4 だけ
        link = Link(left=OUTDOOR, right=_node('サイディング', 0.0, 0.002, 0), interface=_interface(22.4))
        assert link.r_h == pytest.approx(1 / 22.4)

    def test_inside_surface(self):
        # 石膏ボードの右の端と室内：内表面の伝達 1/9.2 だけ
        link = Link(left=_node('石膏ボード', 0.00525, 0.0, 6), right=INDOOR, interface=_interface(9.2))
        assert link.r_h == pytest.approx(1 / 9.2)

    def test_air_layer(self):
        # サイディングの右の端と空気の質点：通気層の中の伝達 1/9.2 だけ
        air = AirNode(thickness=0.025, alpha=0.085, layer_index=1)
        link = Link(left=_node('サイディング', 0.002, 0.0, 0), right=air, interface=_interface(9.2))
        assert link.r_h == pytest.approx(1 / 9.2)

    def test_outdoor_indoor_resistance(self):
        # 外気・室内は質点ではないので、端までの抵抗は 0
        assert (OUTDOOR.r_h_right, INDOOR.r_h_left) == (0.0, 0.0)
        assert OUTDOOR == Outdoor() and INDOOR == Indoor()

    def test_indoor_on_left(self):
        with pytest.raises(ValueError, match='left'):
            Link(left=INDOOR, right=_node('サイディング', 0.0, 0.002, 0), interface=_interface(22.4))

    def test_outdoor_on_right(self):
        with pytest.raises(ValueError, match='right'):
            Link(left=_node('石膏ボード', 0.00525, 0.0, 6), right=OUTDOOR, interface=_interface(9.2))

    def test_surface_without_interface(self):
        with pytest.raises(ValueError, match='表面の伝達'):
            Link(left=OUTDOOR, right=_node('サイディング', 0.0, 0.002, 0), interface=None)

    def test_same_layer_with_interface(self):
        a, b = _node('構造用合板1', 0.0009, 0.0009, 0), _node('構造用合板1', 0.0009, 0.0009, 0)
        with pytest.raises(ValueError, match='同じ層'):
            Link(left=a, right=b, interface=_interface(50.0))

    def test_different_layers_without_interface(self):
        a, b = _node('構造用合板1', 0.0009, 0.0, 0), _node('グラスウール1', 0.0, 0.0125, 1)
        with pytest.raises(ValueError, match='間のつながり'):
            Link(left=a, right=b, interface=None)

    def test_frozen(self):
        link = Link(left=OUTDOOR, right=_node('サイディング', 0.0, 0.002, 0), interface=_interface(22.4))
        with pytest.raises(dataclasses.FrozenInstanceError):
            link.interface = None


def _lambda(name: str) -> float:
    """旧 Excel の材料の熱伝導率, W/(m K)。"""
    return material_finder()('whamcavity', name).lambda_h


def _main_wall_interfaces() -> tuple:
    """Fortran の入力の壁の境目 8 個（外表面、通気層の両側、接触、防湿層、内表面。層と層の境目は水膜の左右 2 組）。"""
    return WallInput.from_dict(_main_wall_dict()).interfaces


class TestMakeLinks:

    def test_one_layer_total(self):
        # 合板 9 mm・分割数 6：つながり 7 個、熱抵抗の合計 = 1/22.4 + 0.009/λ + 1/9.2
        nodes = _make_nodes((_layer('構造用合板1', 0.009, 6),), _one_layer_interfaces(), find_material=material_finder())
        links = _make_links(nodes, _one_layer_interfaces())
        assert len(links) == 7
        assert sum(l.r_h for l in links) == pytest.approx(1 / 22.4 + 0.009 / _lambda('構造用合板1') + 1 / 9.2)

    def test_two_layers_total(self):
        # 合板 9 mm + GW 100 mm、接触 50（水膜の左右に 1/100 ずつ）：接触抵抗の合計 1/50 が 1 回だけ入る
        layers = (_layer('構造用合板1', 0.009, 6), _layer('グラスウール1', 0.1, 5))
        interfaces = (_surf(22.4), _film(50.0), _surf(9.2))
        nodes = _make_nodes(layers, interfaces, find_material=material_finder())
        links = _make_links(nodes, interfaces)
        expected = 1 / 22.4 + 0.009 / _lambda('構造用合板1') + 1 / 50.0 + 0.1 / _lambda('グラスウール1') + 1 / 9.2
        assert sum(l.r_h for l in links) == pytest.approx(expected)

    def test_air_layer_total(self):
        # サイディング + 通気層 + 石膏ボード：通気層の両側は 1/9.2（水膜の左右の合計）、空気の中の抵抗は 0
        layers = (_layer('サイディング', 0.016, 5), _air_layer(), _layer('石膏ボード', 0.042, 5))
        interfaces = (_surf(22.4), _film(9.2), _film(9.2), _surf(9.2))
        nodes = _make_nodes(layers, interfaces, find_material=material_finder())
        links = _make_links(nodes, interfaces)
        expected = 1 / 22.4 + 0.016 / _lambda('サイディング') + 1 / 9.2 + 1 / 9.2 + 0.042 / _lambda('石膏ボード') + 1 / 9.2
        assert sum(l.r_h for l in links) == pytest.approx(expected)

    def test_main_wall_number(self):
        # 今の壁 7 層・質点 39 個：つながり 40 個
        nodes = _make_nodes(_main_wall_layers(), _main_wall_interfaces(), find_material=material_finder())
        assert len(_make_links(nodes, _main_wall_interfaces())) == 40

    def test_each_interface_used_once(self):
        # 境目を持つつながりは 14 個（表面 2 個と、層と層の境目 6 個の水膜の左右）。
        # 入力の境目が外気側から順に（層と層の境目は左、右の順に）1 回ずつ、抵抗（1 / コンダクタンス）に変えて使われる
        interfaces = _main_wall_interfaces()
        nodes = _make_nodes(_main_wall_layers(), interfaces, find_material=material_finder())
        used = [l.interface for l in _make_links(nodes, interfaces) if l.interface is not None]
        expected = [interfaces[0]] + [x for f in interfaces[1:-1] for x in (f.left, f.right)] + [interfaces[-1]]
        assert len(used) == 14
        assert [(u.r_h, u.r_m) for u in used] == pytest.approx([(1 / x.g_h, 1 / x.g_m) for x in expected])

    def test_order(self):
        # k 番目のつながりは、質点 k - 1 と質点 k の間。両端は外気と室内
        nodes = _make_nodes(_main_wall_layers(), _main_wall_interfaces(), find_material=material_finder())
        links = _make_links(nodes, _main_wall_interfaces())
        assert links[0].left is OUTDOOR and links[0].right is nodes[0]
        assert links[-1].left is nodes[-1] and links[-1].right is INDOOR
        for k in range(1, len(nodes)):
            assert links[k].left is nodes[k - 1] and links[k].right is nodes[k]


def _state(rh: float) -> State:
    """20 deg.C・相対湿度 rh % の状態。"""
    return State.init(t=20.0 + T_0, rh=rh)


def _interface_m(g_m: float) -> Interface:
    """湿気コンダクタンス g_m の境目（熱は 50）。計算の側の境目なので抵抗で持つ。"""
    return Interface(r_h=1.0 / 50.0, r_m=1.0 / g_m)


class TestLinkMoisture:

    def test_interior(self):
        # 合板の中、左右とも RH 50%：0.0018 / 湿気伝導率(50%)
        a, b = _node('構造用合板1', 0.0009, 0.0009, 0), _node('構造用合板1', 0.0009, 0.0009, 0)
        link = Link(left=a, right=b, interface=None)
        s = _state(50.0)
        assert link.r_m(s, s) == pytest.approx(0.0018 / a.material.get_lambda_dsh_m(rh=s.rh))

    def test_interior_different_states(self):
        # 左 RH 50%・右 RH 90%：左右の半分は、それぞれの側の状態で計算する
        a, b = _node('構造用合板1', 0.0009, 0.0009, 0), _node('構造用合板1', 0.0009, 0.0009, 0)
        link = Link(left=a, right=b, interface=None)
        s50, s90 = _state(50.0), _state(90.0)
        m = a.material
        assert link.r_m(s50, s90) == pytest.approx(0.0009 / m.get_lambda_dsh_m(rh=s50.rh) + 0.0009 / m.get_lambda_dsh_m(rh=s90.rh))

    @pytest.mark.parametrize('g_m', [7.75e-7, 2.0e-11])
    def test_layer_to_layer(self, g_m):
        # 合板の右の端と GW の左の端：境目の 1 / g_m だけ（接触 7.75e-7、防湿層 2.0e-11）
        a, b = _node('構造用合板1', 0.0009, 0.0, 0), _node('グラスウール1', 0.0, 0.0125, 1)
        link = Link(left=a, right=b, interface=_interface_m(g_m))
        s = _state(50.0)
        assert link.r_m(s, s) == pytest.approx(1 / g_m)

    def test_outside_surface(self):
        link = Link(left=OUTDOOR, right=_node('サイディング', 0.0, 0.002, 0), interface=_interface_m(2.0e-11))
        s = _state(50.0)
        assert link.r_m(s, s) == pytest.approx(1 / 2.0e-11)

    def test_inside_surface(self):
        link = Link(left=_node('石膏ボード', 0.00525, 0.0, 6), right=INDOOR, interface=_interface_m(2.6e-10))
        s = _state(50.0)
        assert link.r_m(s, s) == pytest.approx(1 / 2.6e-10)

    def test_air_layer(self):
        air = AirNode(thickness=0.025, alpha=0.085, layer_index=1)
        link = Link(left=_node('サイディング', 0.002, 0.0, 0), right=air, interface=_interface_m(3.43e-8))
        s = _state(50.0)
        assert link.r_m(s, s) == pytest.approx(1 / 3.43e-8)


class TestLinkLiquid:

    def test_interface(self):
        # 境目のあるつながりは、液水は移動しない（無限大）
        a, b = _node('構造用合板1', 0.0009, 0.0, 0), _node('グラスウール1', 0.0, 0.0125, 1)
        link = Link(left=a, right=b, interface=_interface_m(7.75e-7))
        s = _state(99.0)
        assert link.r_l(s, s) == float('inf')

    def test_interior_no_liquid(self):
        # 構造用合板1 は RH 50% では液水伝導率 0：無限大
        a, b = _node('構造用合板1', 0.0009, 0.0009, 0), _node('構造用合板1', 0.0009, 0.0009, 0)
        link = Link(left=a, right=b, interface=None)
        s = _state(50.0)
        assert link.r_l(s, s) == float('inf')

    def test_interior_liquid(self):
        # 構造用合板1 は RH 99% では液水が動く：0.0018 / 液水伝導率
        a, b = _node('構造用合板1', 0.0009, 0.0009, 0), _node('構造用合板1', 0.0009, 0.0009, 0)
        link = Link(left=a, right=b, interface=None)
        s = _state(99.0)
        lambda_l = a.material.get_lambda_dsh_mu_l(rh=s.rh, t=s.t)
        assert lambda_l > 0.0
        assert link.r_l(s, s) == pytest.approx(0.0018 / lambda_l)


class TestInterface:

    def test_values(self):
        interface = Interface(r_h=1 / 50.0, r_m=1 / 7.75e-7)
        assert (interface.r_h, interface.r_m) == (1 / 50.0, 1 / 7.75e-7)

    @pytest.mark.parametrize('r_h, r_m, key', [(0.0, 1.0e6, 'r_h'), (0.02, -1.0, 'r_m')])
    def test_invalid_value(self, r_h, r_m, key):
        with pytest.raises(ValueError, match=key):
            Interface(r_h=r_h, r_m=r_m)

    def test_frozen(self):
        interface = Interface(r_h=1 / 50.0, r_m=1 / 7.75e-7)
        with pytest.raises(dataclasses.FrozenInstanceError):
            interface.r_h = 0.1


class TestSurface:

    def test_values(self):
        surface = Surface(direction=Direction.E, angle=90.0, albedo=0.3, absorptance=0.9, emissivity=0.9, eva_height=4.2, exposure_factor=1.0, deposition_factor=1.0)
        assert surface.direction == Direction.E
        assert (surface.angle, surface.albedo, surface.absorptance, surface.emissivity, surface.eva_height) == (90.0, 0.3, 0.9, 0.9, 4.2)

    def test_frozen(self):
        surface = Surface(direction=Direction.E, angle=90.0, albedo=0.3, absorptance=0.9, emissivity=0.9, eva_height=4.2, exposure_factor=1.0, deposition_factor=1.0)
        with pytest.raises(dataclasses.FrozenInstanceError):
            surface.angle = 45.0


def _main_wall_dict() -> dict:
    """Fortran の入力の壁の dict（test_wall_input.py の _wall_dict と同じ壁）。"""

    def m(name, thickness, n_div, temperature, rh):
        return {'source': 'whamcavity', 'name': name, 'thickness': thickness, 'n_div': n_div,
                'initial_temperature': temperature, 'initial_rh': rh}

    return {
        'surface': {'direction': 'e', 'angle': 90.0, 'albedo': 0.3, 'absorptance': 0.9, 'emissivity': 0.9, 'eva_height': 4.2, 'exposure_factor': 1.0, 'deposition_factor': 1.0},
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
            {'left': {'g_h': 100.0, 'g_m': 1.55e-6}, 'right': {'g_h': 100.0, 'g_m': 2.0e-11}},
            {},
            {'g_h': 9.2, 'g_m': 2.6e-10},
        ],
        'len_vertical': 7.0,
        'len_horizontal': 0.42,
    }


def _main_wall_geometry() -> WallGeometry:
    return WallGeometry.from_input(WallInput.from_dict(_main_wall_dict()))


class TestWallGeometry:

    def test_main_wall(self):
        # 今の壁 7 層：質点 39 個（水膜 6 個を含む）、つながり 40 個
        g = _main_wall_geometry()
        assert (len(g.nodes), len(g.links)) == (39, 40)

    def test_total_r_h(self):
        # 熱抵抗の合計 = 外表面 + 各層の d/λ + 通気層の両側（表面の伝達と材料の側の水膜のデフォルト）+ 接触 4 か所 + 内表面
        g = _main_wall_geometry()
        layers = [('サイディング', 0.016), ('石膏ボード', 0.042), ('構造用合板1', 0.009),
                  ('グラスウール1', 0.1), ('構造用合板1', 0.012), ('石膏ボード', 0.042)]
        expected = 1 / 22.4 + sum(d / _lambda(name) for name, d in layers) + 2 * (1 / 9.2 + 1 / 1000.0) + 4 / 50.0 + 1 / 9.2
        assert sum(l.r_h for l in g.links) == pytest.approx(expected)

    def test_area(self):
        assert _main_wall_geometry().area == pytest.approx(7.0 * 0.42)

    def test_copied_values(self):
        # 外表面と壁の大きさが、入力と同じ値で写される
        g = _main_wall_geometry()
        assert g.surface == Surface(direction=Direction.E, angle=90.0, albedo=0.3, absorptance=0.9, emissivity=0.9, eva_height=4.2, exposure_factor=1.0, deposition_factor=1.0)
        assert (g.len_vertical, g.len_horizontal) == (7.0, 0.42)

    @pytest.mark.parametrize('direction, angle, expected', [
        ('e', 90.0, 7.0),      # 垂直な壁：縦の長さと同じ
        ('s', 30.0, 3.5),      # 傾斜角 30 deg の屋根：7 × sin 30 deg
        ('top', 0.0, 0.0),     # 上向きの水平：0
        ('bottom', 180.0, 0.0),  # 下向きの水平：0
    ])
    def test_height(self, direction, angle, expected):
        # 壁の下端と上端の高さの差 H = 縦の長さ × sin(傾斜角)
        d = _main_wall_dict()
        d['surface'] = dict(d['surface'], direction=direction, angle=angle)
        g = WallGeometry.from_input(WallInput.from_dict(d))
        assert g.height == pytest.approx(expected, abs=1e-12)

    def test_find_material(self):
        # 材料を探す関数を渡すと、その関数が使われる
        called = []
        excel = material_finder()

        def recording_find(source, name):
            called.append(name)
            return excel(source, name)

        WallGeometry.from_input(WallInput.from_dict(_main_wall_dict()), find_material=recording_find)
        assert called == ['サイディング', '石膏ボード', '構造用合板1', 'グラスウール1', '構造用合板1', '石膏ボード']

    def test_material_not_found(self):
        d = _main_wall_dict()
        d['layers'][3]['name'] = '構造用合板3'
        with pytest.raises(ValueError, match='層 3 番目'):
            WallGeometry.from_input(WallInput.from_dict(d))

    def test_links_too_few(self):
        # 直接作るときに、つながりが 1 つ足りないとエラー
        g = _main_wall_geometry()
        with pytest.raises(ValueError, match='つながりの数'):
            WallGeometry(nodes=g.nodes, links=g.links[:-1], surface=g.surface,
                         len_vertical=g.len_vertical, len_horizontal=g.len_horizontal)

    def test_links_wrong_order(self):
        # 直接作るときに、つながりの順番を入れ替えるとエラー
        g = _main_wall_geometry()
        links = list(g.links)
        links[1], links[2] = links[2], links[1]
        with pytest.raises(ValueError, match='質点 0 番目'):
            WallGeometry(nodes=g.nodes, links=tuple(links), surface=g.surface,
                         len_vertical=g.len_vertical, len_horizontal=g.len_horizontal)


class TestWaterAmountAndHeatCapacity:

    def test_material_water_amount(self):
        # 構造用合板1、20 deg.C・RH 50%：乾燥密度 × u(50%) × 幅
        m = _plywood()
        node = MaterialNode(d_left=0.0009, d_right=0.0009, material=m, layer_index=0)
        state = _state(50.0)
        assert node.water_amount(state) == pytest.approx(m.rho_s * m.get_u(50.0) * 0.0018)

    def test_air_water_amount(self):
        # 通気層 25 mm、20 deg.C・50%：水蒸気の密度 p_v / (R_v T) × 厚さ（約 0.0002 kg/m2）
        node = AirNode(thickness=0.025, alpha=0.085, layer_index=1)
        state = _state(50.0)
        assert node.water_amount(state) == pytest.approx(state.p_v / (R_V * state.t) * 0.025)
        assert node.water_amount(state) == pytest.approx(0.0002, rel=0.1)

    def test_material_heat_capacity(self):
        # 乾燥した値 + 水の比熱 × 水分量
        node = MaterialNode(d_left=0.0009, d_right=0.0009, material=_plywood(), layer_index=0)
        state = _state(80.0)
        assert node.heat_capacity(state) == pytest.approx(node.heat_capacity_dry + C_W * node.water_amount(state))
        assert node.heat_capacity(state) > node.heat_capacity_dry

    def test_air_heat_capacity(self):
        # 空気の質点：空気の比熱 × その温度の空気の密度 × 厚さ（水蒸気の分は足さない）
        node = AirNode(thickness=0.025, alpha=0.085, layer_index=1)
        state = _state(80.0)
        assert node.heat_capacity(state) == pytest.approx(C_A * get_rho_a(state.t) * 0.025)

        # 20 deg.C では、およそ 1005 × 1.204 × 0.025 ≈ 30.3 J/(m2 K)
        assert node.heat_capacity(state) == pytest.approx(30.3, rel=1e-2)

        # 温度が高いほど、空気の密度が小さいので熱容量も小さい
        assert node.heat_capacity(State.init(t=40.0 + T_0, rh=80.0)) < node.heat_capacity(state)


class TestVentilationRate:

    def test_hand_calculation(self):
        # 厚さ 25 mm、α 0.085、高さ 7 m、縦 7 m、外気 0 deg.C、通気層 20 deg.C：
        # q_v = α d / L_v × √(2 g H |ρ_o - ρ_i| / ρ_o)、ρ = 353 / T
        node = AirNode(thickness=0.025, alpha=0.085, layer_index=1)
        rho_o, rho_i = 353.0 / 273.15, 353.0 / 293.15
        expected = 0.085 * 0.025 / 7.0 * (2 * 9.8 * 7.0 * (rho_o - rho_i) / rho_o) ** 0.5
        assert node.ventilation_rate(293.15, 273.15, height=7.0, len_vertical=7.0) == pytest.approx(expected)
        # 目安：√ の中は 2 × 9.8 × 7 × 0.068 ≈ 9.4 で、q_v ≈ 3.0e-4 × 3.06 ≈ 9.3e-4 m3/(m2 s)（通気層の中の風速は α × 3.06 ≈ 0.26 m/s）
        assert expected == pytest.approx(9.3e-4, rel=0.01)

    def test_same_temperature(self):
        # 外気と同じ温度なら浮力がないので 0
        node = AirNode(thickness=0.025, alpha=0.085, layer_index=1)
        assert node.ventilation_rate(283.15, 283.15, height=7.0, len_vertical=7.0) == 0.0

    def test_direction_does_not_matter(self):
        # 通気層が外気より冷たくても（夜）、温度差の大きさが同じなら同じ量（向きは考えない）
        node = AirNode(thickness=0.025, alpha=0.085, layer_index=1)
        warm = node.ventilation_rate(293.15, 283.15, height=7.0, len_vertical=7.0)
        cold = node.ventilation_rate(273.15, 283.15, height=7.0, len_vertical=7.0)
        assert warm == pytest.approx(cold, rel=0.1)

    def test_no_opening(self):
        # 流量係数 0（換気のない空気層）なら 0
        node = AirNode(thickness=0.025, alpha=0.0, layer_index=1)
        assert node.ventilation_rate(293.15, 273.15, height=7.0, len_vertical=7.0) == 0.0
