import numpy as np

from modules.wall_builder import build_wall
from modules.wall_geometry import AirNode, material_finder
from modules.wall_input import WallInput


def _main_wall_dict() -> dict:
    """Fortran の入力の壁の dict（test_wall_geometry.py の _main_wall_dict と同じ壁）。"""

    def m(name, thickness, n_div, temperature, rh):
        return {'source': 'whamcavity', 'name': name, 'thickness': thickness, 'n_div': n_div,
                'initial_temperature': temperature, 'initial_rh': rh}

    return {
        'surface': {'direction': 'e', 'angle': 90.0, 'albedo': 0.3, 'absorptance': 0.9, 'emissivity': 0.9, 'exposure_factor': 1.0, 'deposition_factor': 1.0},
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


def _build():
    return build_wall(WallInput.from_dict(_main_wall_dict()))


class TestBuildWall:

    def test_size(self):
        # 今の壁：質点 39 個（水膜 6 個を含む）、状態の長さも 39
        geometry, state = _build()
        assert (len(geometry.nodes), len(state)) == (39, 39)

    def test_initial_temperature(self):
        # サイディング（質点 0〜4）・その右の水膜（5）・通気層（6）・その右の水膜（7）は 15 deg.C、それ以外は 10 deg.C
        # （水膜の質点の初期値は、左の層の値）
        _, state = _build()
        expected = [15.0] * 8 + [10.0] * 31
        np.testing.assert_allclose(state.theta_is, expected)

    def test_initial_rh(self):
        # 合板（質点 14〜19、27〜32）とその右の水膜（20、33）は 80%、それ以外は 50%
        _, state = _build()
        expected = np.full(39, 50.0)
        expected[14:21] = 80.0
        expected[27:34] = 80.0
        np.testing.assert_allclose(state.rh_is, expected)

    def test_air_node(self):
        # 質点 6 は空気の質点で、初期値は通気層の値（15 deg.C・50%）
        geometry, state = _build()
        assert isinstance(geometry.nodes[6], AirNode)
        state_i = state.state(6)
        assert (round(state_i.theta, 10), round(state_i.rh, 10)) == (15.0, 50.0)

    def test_find_material(self):
        # 渡した材料を探す関数が、壁の形を作るときに使われる
        called = []
        excel = material_finder()

        def recording_find(source, name):
            called.append(name)
            return excel(source, name)

        build_wall(WallInput.from_dict(_main_wall_dict()), find_material=recording_find)
        assert len(called) == 6
