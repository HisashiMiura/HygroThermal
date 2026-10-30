import numpy as np
import pytest

from modules.balance import heat_capacity_is, heat_flow_ks, liquid_flow_ks, vapour_flow_ks, ventilation_rate_is, water_amount_is
from modules.boundary import IndoorBoundary, OutdoorBoundary, SolarRadiation
from modules.newton import NewtonNotConverged
from modules.solver import HEAT_RESIDUAL_TOL, MOISTURE_RESIDUAL_TOL, Condensation, step
from modules.state import State
from modules.thermo_dynamics import ATP, GOFF, ROW_CP, RW, get_wp
from modules.wall_builder import build_wall
from modules.wall_geometry import AirNode, WallGeometry
from modules.wall_input import WallInput
from modules.wall_state import WallState

_SURFACE = {'direction': 'e', 'angle': 90.0, 'albedo': 0.3, 'absorptance': 0.9, 'emissivity': 0.9, 'eva_height': 4.2}


def _outdoor(theta: float, rh: float) -> OutdoorBoundary:
    """外気温 theta（deg.C）・相対湿度 rh（%）の外気の境界条件。"""
    return OutdoorBoundary.from_theta_rh(theta=theta, rh=rh)


def _indoor(theta: float, rh: float) -> IndoorBoundary:
    """室温 theta（deg.C）・相対湿度 rh（%）の室内の境界条件。"""
    return IndoorBoundary.from_theta_rh(theta=theta, rh=rh)


def _one_layer(name: str, thickness: float, n_div: int, theta: float, rh: float) -> tuple[WallGeometry, WallState]:
    """1 層の壁（外表面 22.4・2.0e-11、内表面 9.2・2.6e-10）。初期は theta（deg.C）・rh（%）で一様。"""
    return build_wall(WallInput.from_dict({
        'surface': _SURFACE,
        'layers': [{'source': 'whamcavity', 'name': name, 'thickness': thickness, 'n_div': n_div,
                    'initial_temperature': theta, 'initial_rh': rh}],
        'interfaces': [{'cond_h': 22.4, 'cond_m': 2.0e-11}, {'cond_h': 9.2, 'cond_m': 2.6e-10}],
        'len_vertical': 1.0, 'len_horizontal': 1.0, 'height': 1.0,
    }))


def _main_wall() -> tuple[WallGeometry, WallState]:
    """今の main.py・Fortran の入力の壁 7 層（通気層を含む。換気はまだ扱わない）。"""

    def m(name, thickness, n_div, temperature, rh):
        return {'source': 'whamcavity', 'name': name, 'thickness': thickness, 'n_div': n_div,
                'initial_temperature': temperature, 'initial_rh': rh}

    return build_wall(WallInput.from_dict({
        'surface': _SURFACE,
        'layers': [
            m('サイディング', 0.016, 5, 15.0, 50.0),
            {'air_layer': True, 'thickness': 0.025, 'alpha': 0.085, 'initial_temperature': 15.0, 'initial_rh': 50.0},
            m('石膏ボード', 0.042, 5, 10.0, 50.0),
            m('構造用合板1', 0.009, 6, 10.0, 80.0),
            m('グラスウール1', 0.1, 5, 10.0, 50.0),
            m('構造用合板1', 0.012, 6, 10.0, 80.0),
            m('石膏ボード', 0.042, 5, 10.0, 50.0),
        ],
        'interfaces': [{'cond_h': h, 'cond_m': w} for h, w in (
            (22.4, 2.0e-11), (9.2, 3.43e-8), (9.2, 3.43e-8), (50.0, 7.75e-7),
            (50.0, 7.75e-7), (50.0, 2.0e-11), (50.0, 7.75e-7), (9.2, 2.6e-10),
        )],
        'len_vertical': 7.0, 'len_horizontal': 0.42, 'height': 7.0,
    }))


def _moisture_flow_ks(g, state, outdoor, indoor):
    """つながりごとの水分の流れ（水蒸気 + 液水）, kg/(m2 s), [k]。"""
    return (vapour_flow_ks(g, state.t_is, state.mu_is, outdoor, indoor)
            + liquid_flow_ks(g, state.t_is, state.mu_is, outdoor, indoor))


def _total_heat_flow_ks(g, state, outdoor, indoor):
    """つながりごとの熱の流れ（熱伝導 + 水蒸気が運ぶ潜熱）, W/m2, [k]。"""
    return (heat_flow_ks(g, state.t_is, outdoor, indoor)
            + RW * vapour_flow_ks(g, state.t_is, state.mu_is, outdoor, indoor))


class TestStep:

    def test_uniform(self):
        # 外気・室内・壁がすべて 20 deg.C・50% なら、1 ステップ進めても変わらない
        g, state_n = _one_layer('構造用合板1', 0.009, 6, 20.0, 50.0)
        state, _ = step(g, state_n, dt=180.0, outdoor=_outdoor(20.0, 50.0), indoor=_indoor(20.0, 50.0))
        np.testing.assert_allclose(state.theta_is, 20.0)
        np.testing.assert_allclose(state.rh_is, 50.0)

    def test_iterations(self):
        # 非線形だが、数回（5 回以下）の修正で収束する
        g, state_n = _one_layer('構造用合板1', 0.009, 6, 10.0, 50.0)
        _, m = step(g, state_n, dt=180.0, outdoor=_outdoor(0.0, 80.0), indoor=_indoor(20.0, 50.0))
        assert 1 <= m <= 5

    def test_conservation(self):
        # 1 ステップで壁にたまった水分 = (外気から入った水分 - 室内へ出た水分) × Δt。
        # 熱も同じ（潜熱を含む。たまった熱は、次の時刻の熱容量 × 温度の変化）。
        # ずれてよいのは、ニュートン法の許容値の分（質点の数 × 許容値 × Δt）まで
        g, state_n = _one_layer('構造用合板1', 0.009, 6, 10.0, 50.0)
        dt = 180.0
        outdoor, indoor = _outdoor(0.0, 80.0), _indoor(20.0, 50.0)
        state, _ = step(g, state_n, dt=dt, outdoor=outdoor, indoor=indoor)

        w_stored = (water_amount_is(g, state.t_is, state.mu_is) - water_amount_is(g, state_n.t_is, state_n.mu_is)).sum()
        j_ks = _moisture_flow_ks(g, state, outdoor, indoor)
        assert w_stored == pytest.approx((j_ks[0] - j_ks[-1]) * dt, rel=0.0, abs=len(g.nodes) * MOISTURE_RESIDUAL_TOL * dt)

        h_stored = (heat_capacity_is(g, state.t_is, state.mu_is) * (state.t_is - state_n.t_is)).sum()
        q_ks = _total_heat_flow_ks(g, state, outdoor, indoor)
        assert h_stored == pytest.approx((q_ks[0] - q_ks[-1]) * dt, rel=0.0, abs=len(g.nodes) * HEAT_RESIDUAL_TOL * dt)

    def test_absorption(self):
        # 構造用合板1 の 1 層、20 deg.C 一定、初期 50%、外気・室内とも 70%。
        # 1 年（Δt = 1 日）進めると、どの質点も 70% になり、水分量の増えた分が、吸った水分の合計と一致する
        g, state_0 = _one_layer('構造用合板1', 0.009, 6, 20.0, 50.0)
        outdoor, indoor = _outdoor(20.0, 70.0), _indoor(20.0, 70.0)
        dt, n_steps = 86400.0, 365
        state, absorbed = state_0, 0.0
        for _ in range(n_steps):
            state, _ = step(g, state, dt=dt, outdoor=outdoor, indoor=indoor)
            j_ks = _moisture_flow_ks(g, state, outdoor, indoor)
            absorbed += (j_ks[0] - j_ks[-1]) * dt
        np.testing.assert_allclose(state.rh_is, 70.0, atol=1e-4)
        np.testing.assert_allclose(state.theta_is, 20.0, atol=1e-6)
        w_increase = (water_amount_is(g, state.t_is, state.mu_is) - water_amount_is(g, state_0.t_is, state_0.mu_is)).sum()
        assert w_increase > 0.0
        assert w_increase == pytest.approx(absorbed, rel=0.0, abs=n_steps * len(g.nodes) * MOISTURE_RESIDUAL_TOL * dt)

    def test_steady_vapour(self):
        # 湿気伝導率が一定の構造用合板2、20 deg.C 一定、外気 30%・室内 70%。
        # 1 年（Δt = 1 日）進めると定常になり、どのつながりの水蒸気の流れも (p_v,外 - p_v,内) / 透湿抵抗の合計
        g, state = _one_layer('構造用合板2', 0.012, 5, 20.0, 50.0)
        outdoor, indoor = _outdoor(20.0, 30.0), _indoor(20.0, 70.0)
        for _ in range(365):
            state, _ = step(g, state, dt=86400.0, outdoor=outdoor, indoor=indoor)
        state_any = State.init(t=20.0 + ATP, rh=50.0)
        r_m_total = sum(link.r_m(state_any, state_any) for link in g.links)
        j_v_ks = vapour_flow_ks(g, state.t_is, state.mu_is, outdoor, indoor)
        np.testing.assert_allclose(j_v_ks, (outdoor.state.p_v - indoor.state.p_v) / r_m_total, rtol=1e-4)

    def test_steady_main_wall(self):
        # 今の壁 7 層、外気 0 deg.C・80%、室内 20 deg.C・50%（日射なし。通気層は温度差で換気される）。
        # 1 年（Δt = 1 日）進めると、熱は定常になる：通気層より外側・内側それぞれで、どのつながりの熱の流れ（潜熱を含む）も同じ値になり、
        # 通気層の空気の質点では、熱伝導で入る熱 - 出ていく熱 + 換気で入る熱 = 0（空気の質点は相変化しないので潜熱は数えない）。
        # 水分は両側の湿気の抵抗が大きく、1 年では定常にならないので、結露していない（100% 未満）ことだけ確かめる
        g, state = _main_wall()
        outdoor, indoor = _outdoor(0.0, 80.0), _indoor(20.0, 50.0)
        for _ in range(365):
            state, _ = step(g, state, dt=86400.0, outdoor=outdoor, indoor=indoor)
        assert state.rh_is.max() < 100.0

        # 通気層の空気の質点の番号（つながり i_air は外側から入る所、i_air + 1 は内側へ出る所）
        i_air = next(i for i, node in enumerate(g.nodes) if isinstance(node, AirNode))
        q_ks = _total_heat_flow_ks(g, state, outdoor, indoor)
        np.testing.assert_allclose(q_ks[:i_air + 1], q_ks[:i_air + 1].mean(), rtol=1e-5)
        np.testing.assert_allclose(q_ks[i_air + 1:], q_ks[i_air + 1:].mean(), rtol=1e-5)

        # 通気層の熱の収支（熱伝導だけ + 換気）。換気で外気に熱が逃げるので、内側から入る熱のほうが外側へ出る熱より大きい
        q_h_ks = heat_flow_ks(g, state.t_is, outdoor, indoor)
        q_v_is = ventilation_rate_is(g, state.t_is, outdoor)
        q_vent = ROW_CP * q_v_is[i_air] * (outdoor.state.t - state.t_is[i_air])
        assert q_v_is[i_air] > 0.0
        assert q_h_ks[i_air] - q_h_ks[i_air + 1] + q_vent == pytest.approx(0.0, abs=1e-4)

    def test_condensation(self):
        # 構造用合板1 の 1 層、外気 -5 deg.C・90%、室内 25 deg.C・95%。
        # 室内側から入る水蒸気を外表面から逃がせず、外側の質点（0 番目）が飽和すると、結露として止まる
        g, state = _one_layer('構造用合板1', 0.009, 6, 0.0, 95.0)
        outdoor, indoor = _outdoor(-5.0, 90.0), _indoor(25.0, 95.0)
        with pytest.raises(Condensation, match='質点 0（層 0'):
            for _ in range(24 * 30):
                state, _ = step(g, state, dt=3600.0, outdoor=outdoor, indoor=indoor)

    def test_condensation_is_not_converged(self):
        # 結露は、収束しなかったことの一種（NewtonNotConverged で受けても捕まえられる）
        assert issubclass(Condensation, NewtonNotConverged)

    def test_sol_air_temperature(self):
        # 日射ありの計算と、日射なしで外気温を q_sol / α_o（相当外気温度）にした計算は、毎ステップ一致する。
        # 外気の水蒸気圧はそろえる。外表面の質点の熱の収支が、式の上で同じになるため
        g, state_a = _one_layer('構造用合板1', 0.009, 6, 10.0, 50.0)
        state_b = state_a
        solar = SolarRadiation(sin_h=0.5, cos_h=3 ** 0.5 / 2, sin_a=0.0, cos_a=1.0, i_dn=600.0, i_sky=150.0, r_eff=40.0)
        q_sol = g.surface.solar_heat_gain(solar)
        assert q_sol > 0.0

        # 外表面の熱伝達率 α_o, W/(m2 K)（外気と質点 0 のつながりの熱抵抗の逆数）
        alpha_o = 1.0 / g.links[0].r_h

        outdoor_a = OutdoorBoundary(state=_outdoor(0.0, 80.0).state, solar=solar)
        t_sol_air = outdoor_a.state.t + q_sol / alpha_o
        outdoor_b = OutdoorBoundary(state=State.init(t=t_sol_air, rh=100.0 * outdoor_a.state.p_v / GOFF(t_sol_air)[1]))
        assert outdoor_b.state.p_v == pytest.approx(outdoor_a.state.p_v)
        indoor = _indoor(20.0, 50.0)

        # 比べるための、日射なし・外気温もそのままの計算
        state_c = state_a
        outdoor_c = OutdoorBoundary(state=outdoor_a.state)

        for _ in range(20):
            state_a, _ = step(g, state_a, dt=180.0, outdoor=outdoor_a, indoor=indoor)
            state_b, _ = step(g, state_b, dt=180.0, outdoor=outdoor_b, indoor=indoor)
            state_c, _ = step(g, state_c, dt=180.0, outdoor=outdoor_c, indoor=indoor)
        np.testing.assert_allclose(state_a.t_is, state_b.t_is, rtol=0.0, atol=1e-6)
        np.testing.assert_allclose(state_a.rh_is, state_b.rh_is, rtol=0.0, atol=1e-6)

        # 日射の分だけ、外表面の質点は日射なしより温かい
        assert state_a.t_is[0] > state_c.t_is[0] + 1.0


def _air_layer_wall(alpha: float) -> tuple[WallGeometry, WallState]:
    """3 層の壁：構造用合板2 12 mm ＋ 通気層 25 mm（流量係数 alpha）＋ 構造用合板2 12 mm。初期は 20 deg.C・50%。"""
    m = {'source': 'whamcavity', 'name': '構造用合板2', 'thickness': 0.012, 'n_div': 2, 'initial_temperature': 20.0, 'initial_rh': 50.0}
    return build_wall(WallInput.from_dict({
        'surface': _SURFACE,
        'layers': [m, {'air_layer': True, 'thickness': 0.025, 'alpha': alpha, 'initial_temperature': 20.0, 'initial_rh': 50.0}, m],
        'interfaces': [{'cond_h': 22.4, 'cond_m': 2.0e-11}, {'cond_h': 9.2, 'cond_m': 3.43e-8},
                       {'cond_h': 9.2, 'cond_m': 3.43e-8}, {'cond_h': 9.2, 'cond_m': 2.6e-10}],
        'len_vertical': 7.0, 'len_horizontal': 0.42, 'height': 7.0,
    }))


class TestVentilation:

    def test_large_ventilation(self):
        # 換気量をとても大きくすると（流量係数 100）、通気層の温度は外気温に近づく。換気なし（0）では近づかない
        outdoor, indoor = _outdoor(0.0, 50.0), _indoor(20.0, 50.0)
        theta_air = {}
        for alpha in (0.0, 100.0):
            g, state = _air_layer_wall(alpha)
            for _ in range(48):
                state, _ = step(g, state, dt=1800.0, outdoor=outdoor, indoor=indoor)
            theta_air[alpha] = state.theta_is[2]
        assert abs(theta_air[100.0] - 0.0) < 0.5
        assert theta_air[0.0] > 5.0

    def test_drying(self):
        # 乾いた外気（20 deg.C・20%）で換気すると、通気層の湿度が、換気なしより下がる。
        # 温度差がないと浮力がないので、外気を少し冷たく（15 deg.C）して換気させる
        outdoor, indoor = _outdoor(15.0, 20.0), _indoor(20.0, 50.0)
        rh_air = {}
        for alpha in (0.0, 0.085):
            g, state = _air_layer_wall(alpha)
            for _ in range(24):
                state, _ = step(g, state, dt=3600.0, outdoor=outdoor, indoor=indoor)
            rh_air[alpha] = state.rh_is[2]
        assert rh_air[0.085] < rh_air[0.0] - 1.0

    def test_conservation(self):
        # 1 ステップで壁にたまった水分・熱 = 外気から入った分 - 室内へ出た分 + 換気で入った分（換気量は前の時刻の温度で決める）
        g, state_n = _air_layer_wall(0.085)
        dt = 600.0
        outdoor, indoor = _outdoor(0.0, 80.0), _indoor(20.0, 50.0)
        state, _ = step(g, state_n, dt=dt, outdoor=outdoor, indoor=indoor)
        q_v_is = ventilation_rate_is(g, state_n.t_is, outdoor)
        rho_v_o = outdoor.state.p_v / (461.5 * outdoor.state.t)
        rho_v_is = np.array([state.state(i).p_v for i in range(len(state))]) / (461.5 * state.t_is)

        w_stored = (water_amount_is(g, state.t_is, state.mu_is) - water_amount_is(g, state_n.t_is, state_n.mu_is)).sum()
        j_ks = _moisture_flow_ks(g, state, outdoor, indoor)
        j_vent = (q_v_is * (rho_v_o - rho_v_is)).sum()
        assert w_stored == pytest.approx((j_ks[0] - j_ks[-1] + j_vent) * dt, rel=0.0, abs=len(g.nodes) * MOISTURE_RESIDUAL_TOL * dt)

        # 熱：凝縮熱は材料の質点だけで数えるので、各質点の凝縮熱の合計を足す
        h_stored = (heat_capacity_is(g, state.t_is, state.mu_is) * (state.t_is - state_n.t_is)).sum()
        q_h_ks = heat_flow_ks(g, state.t_is, outdoor, indoor)
        j_v_ks = vapour_flow_ks(g, state.t_is, state.mu_is, outdoor, indoor)
        is_material_is = np.array([not isinstance(node, AirNode) for node in g.nodes])
        q_cond = (RW * (j_v_ks[:-1] - j_v_ks[1:]))[is_material_is].sum()
        q_vent = (ROW_CP * q_v_is * (outdoor.state.t - state.t_is)).sum()
        assert h_stored == pytest.approx((q_h_ks[0] - q_h_ks[-1] + q_cond + q_vent) * dt, rel=0.0, abs=len(g.nodes) * HEAT_RESIDUAL_TOL * dt)
