import numpy as np
import pytest

from modules import decay
from modules.balance import liquid_flow_ks, vapour_flow_ks, water_amount_is
from modules.boundary import Boundaries, IndoorBoundary, OutdoorBoundary
from modules.simulation import simulate
from modules.solver import MOISTURE_RESIDUAL_TOL
from modules.wall_builder import build_wall
from modules.wall_input import WallInput

_SURFACE = {'direction': 'e', 'angle': 90.0, 'albedo': 0.3, 'absorptance': 0.9, 'emissivity': 0.9, 'eva_height': 4.2,
            'exposure_factor': 1.0, 'deposition_factor': 1.0}


def _wall(decay_layers=(0, 2)):
    """3 層の壁：構造用合板1（分割数 4、層 0）＋ 通気層（層 1）＋ 木材（分割数 5、層 2）。

    質点は 0〜3 が合板、4 が水膜、5 が空気、6 が水膜、7〜11 が木材。decay_layers の層で木材腐朽を求める。
    """
    m = lambda name, d, n, l: {'source': 'whamcavity', 'name': name, 'thickness': d, 'n_div': n,
                               'initial_temperature': 20.0, 'initial_rh': 50.0, 'wood_decay': l in decay_layers}
    g, _ = build_wall(WallInput.from_dict({
        'surface': _SURFACE,
        'layers': [m('構造用合板1', 0.012, 4, 0),
                   {'air_layer': True, 'thickness': 0.02, 'alpha': 0.0, 'initial_temperature': 20.0, 'initial_rh': 50.0},
                   m('木材', 0.02, 5, 2)],
        'interfaces': [{'g_h': 22.4, 'g_m': 2.0e-11}, {}, {}, {'g_h': 9.2, 'g_m': 2.6e-10}],
        'len_vertical': 1.0, 'len_horizontal': 1.0,
    }))
    return g


class TestFunctions:

    @pytest.mark.parametrize('theta, expected', [(-10.0, 100.0), (0.0, 100.0), (10.0, 95.0), (15.0, 92.5), (25.0, 92.5)])
    def test_rh_critical(self, theta, expected):
        # 15 deg.C 以下は 100 - 0.5 θ（上限 100%）、それを超えると 92.5%
        assert decay.rh_critical(theta) == pytest.approx(expected)

    def test_reaction_rate(self):
        # 20 deg.C で約 1.33e-8 1/s（齋藤ら 2008 の式 (17)）、1.3〜2.7 deg.C で式は負になるが 0 にする
        assert decay.reaction_rate(20.0) == pytest.approx(1.3297e-8, rel=1e-4)
        assert decay.reaction_rate(2.0) == 0.0
        assert decay.reaction_rate(3.0) > 0.0

    @pytest.mark.parametrize('theta, rh, months', [
        (5.0, 98.0, 4.4640),   # -g / f が正：その値
        (20.0, 98.5, 0.5),     # -g / f が正だが小さい（約 0.006 か月）：下限 0.5 か月
        (25.0, 98.0, 0.5),     # -g / f が -0.59 より大きい負（約 -0.30）：0 とみなし、下限 0.5 か月
        (30.0, 100.0, 100.0),  # -g / f が -0.59 より小さい：発芽しないとみなす（100 か月）
    ])
    def test_initial_response_time(self, theta, rh, months):
        assert decay.initial_response_time(theta, rh) == pytest.approx(months * 30 * 86400, rel=1e-4)


class TestDecayNodes:

    def test_layers(self):
        # 合板（層 0）と木材（層 2）で求める：表面は各層の両端、隣は同じ層の中だけ
        nodes = decay.DecayNodes.from_geometry(_wall())
        assert list(np.flatnonzero(nodes.is_target_is)) == [0, 1, 2, 3, 7, 8, 9, 10, 11]
        assert list(np.flatnonzero(nodes.is_surface_is)) == [0, 3, 7, 11]
        assert (nodes.i_left_is[0], nodes.i_right_is[0]) == (-1, 1)
        assert (nodes.i_left_is[3], nodes.i_right_is[3]) == (2, -1)
        assert (nodes.i_left_is[7], nodes.i_right_is[7]) == (-1, 8)
        assert nodes.any

    def test_none(self):
        # 木材腐朽を求める層がなければ、どの質点も求めない
        nodes = decay.DecayNodes.from_geometry(_wall(decay_layers=()))
        assert not nodes.any


class TestAdvance:

    def _run(self, theta, rh, days):
        """木材（層 2）で求め、一様な温度・相対湿度で、1 日ずつ days 日進める。"""
        g = _wall(decay_layers=(2,))
        nodes = decay.DecayNodes.from_geometry(g)
        n = len(g.nodes)
        state = decay.DecayState.initial(n)
        for _ in range(days):
            state = decay.advance(state, nodes, np.full(n, theta), np.full(n, rh), 86400.0)
        return state

    def test_germination(self):
        # 20 deg.C・98.5% では、発芽までの時間は 15 日：15 日目までは発芽せず、16 日目に表面の質点だけ発芽する
        assert not self._run(20.0, 98.5, 15).germinated_is.any()
        state = self._run(20.0, 98.5, 16)
        assert list(np.flatnonzero(state.germinated_is)) == [7, 11]

    def test_spread(self):
        # 発芽した日から表面で腐朽が進み、次の日から 1 日に 1 つずつ内部へ広がる（隣の質量減少率が始めに 0 より大きいとき）
        k = decay.reaction_rate(20.0) * 86400.0
        state = self._run(20.0, 98.5, 17)
        np.testing.assert_allclose(state.mass_loss_is[7:12], [2 * k, k, 0.0, k, 2 * k])
        assert state.mass_loss_is[:7].sum() == 0.0

    def test_dry_no_decay(self):
        # 相対湿度が φ_g（98%）より低ければ、発芽しても腐朽は進まない（20 deg.C・95% は φ_c 92.5% 以上なので発芽はする）
        state = self._run(20.0, 95.0, 60)
        assert state.germinated_is[7]
        assert state.mass_loss_is.sum() == 0.0

    def test_reset(self):
        # 相対湿度が φ_c より低い日があると、続いた時間は 0 に戻る
        g = _wall(decay_layers=(2,))
        nodes = decay.DecayNodes.from_geometry(g)
        n = len(g.nodes)
        state = decay.DecayState.initial(n)
        for rh in [98.5] * 10 + [80.0] + [98.5] * 10:
            state = decay.advance(state, nodes, np.full(n, 20.0), np.full(n, rh), 86400.0)
        assert state.time_wet_is[7] == pytest.approx(10 * 86400.0)
        assert not state.germinated_is.any()

    def test_mass_loss_max(self):
        # 質量減少率は 0.6 で止まる
        state = self._run(28.0, 99.0, 800)
        assert state.mass_loss_is.max() == pytest.approx(decay.MASS_LOSS_MAX)


class TestWaterGeneration:

    def test_value(self):
        # 腐朽で生まれる水 = h ρ_w ΔL / Δt × 幅（齋藤ら 2008 の式 (16)）
        g = _wall(decay_layers=(2,))
        nodes = decay.DecayNodes.from_geometry(g)
        n = len(g.nodes)
        state_nmns = decay.DecayState.initial(n)
        mass_loss_is = np.zeros(n)
        mass_loss_is[7] = 1.0e-4
        state_n = decay.DecayState(time_wet_is=np.zeros(n), germinated_is=np.zeros(n, dtype=bool), mass_loss_is=mass_loss_is)
        j_is = decay.water_generation_is(nodes, state_n, state_nmns, 3600.0)
        assert j_is[7] == pytest.approx(0.319 * 998.0 * 1.0e-4 / 3600.0 * g.nodes[7].width)
        assert np.count_nonzero(j_is) == 1


def _panel(wood_decay=True):
    """厚さ 20 mm の木材の単板（分割数 11）。外気側はほぼ断湿、室内側は 20 deg.C・98.5%（齋藤ら 2008 の case1 に近い）。"""
    g, state_0 = build_wall(WallInput.from_dict({
        'surface': _SURFACE,
        'layers': [{'source': 'whamcavity', 'name': '木材', 'thickness': 0.02, 'n_div': 11,
                    'initial_temperature': 20.0, 'initial_rh': 97.0, 'wood_decay': wood_decay}],
        'interfaces': [{'g_h': 9.3, 'g_m': 1.0e-16}, {'g_h': 9.3, 'g_m': 3.43e-8}],
        'len_vertical': 1.0, 'len_horizontal': 1.0,
    }))
    b = Boundaries(outdoor=OutdoorBoundary.from_theta_rh(theta=20.0, rh=50.0), indoor=IndoorBoundary.from_theta_rh(theta=20.0, rh=98.5))
    return g, state_0, (lambda time: b)


class TestSimulation:

    def test_decay_and_conservation(self):
        # 120 日（時間刻み 1 日）進めると、室内側の表面から腐朽が進む。
        # 水分は、腐朽で生まれた水を含めて、許容値の分（ステップ数 × 質点の数 × 許容値 × Δt）まで保存される
        g, state_0, boundaries_at = _panel()
        dt, n_steps = 86400.0, 120
        r = simulate(g, state_0, dt=dt, n_steps=n_steps, boundaries_at=boundaries_at)
        assert r.mass_loss_is_ns[-1, -1] > 0.0
        assert r.mass_loss_is_ns[0, -1] == 0.0
        assert r.decay_last.germinated_is[-1]

        # 外気・室内から入った水分（ステップごとの、終わりの状態の流れ）
        absorbed = 0.0
        for n in range(1, n_steps + 1):
            j_ks = (vapour_flow_ks(g, r.t_is_ns[:, n], r.mu_is_ns[:, n], boundaries_at(0).outdoor, boundaries_at(0).indoor)
                    + liquid_flow_ks(g, r.t_is_ns[:, n], r.mu_is_ns[:, n], boundaries_at(0).outdoor, boundaries_at(0).indoor))
            absorbed += (j_ks[0] - j_ks[-1]) * dt

        # 腐朽で生まれた水：ステップ n + 1 に、ステップ n - 1 から n までの増えの分が入るので、最後のステップの 1 つ前までの質量減少率の分
        width_is = np.array([node.width for node in g.nodes])
        generated = (0.319 * 998.0 * r.mass_loss_is_ns[:, -2] * width_is).sum()
        assert generated > 0.0

        stored = (water_amount_is(g, r.t_is_ns[:, -1], r.mu_is_ns[:, -1]) - water_amount_is(g, state_0.t_is, state_0.mu_is)).sum()
        assert abs(stored - absorbed - generated) <= n_steps * len(g.nodes) * MOISTURE_RESIDUAL_TOL * dt

    def test_generation_wets(self):
        # 腐朽で生まれた水の分、木材腐朽を求めない場合より水分量が多い
        g, state_0, boundaries_at = _panel()
        g_off, state_0_off, _ = _panel(wood_decay=False)
        r = simulate(g, state_0, dt=86400.0, n_steps=120, boundaries_at=boundaries_at)
        r_off = simulate(g_off, state_0_off, dt=86400.0, n_steps=120, boundaries_at=boundaries_at)
        w = water_amount_is(g, r.t_is_ns[:, -1], r.mu_is_ns[:, -1]).sum()
        w_off = water_amount_is(g_off, r_off.t_is_ns[:, -1], r_off.mu_is_ns[:, -1]).sum()
        assert w > w_off
        assert r_off.decay_last is None
        assert not r_off.mass_loss_is_ns.any()
