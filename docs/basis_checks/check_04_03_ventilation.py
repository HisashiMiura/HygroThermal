"""04_03「通気層の換気のモデル」：表 1（壁 1 m2 あたりの換気量の例）と、8.2 節（換気量をニュートン法の中で計算し直す場合）。"""

import numpy as np

import _common as c
import modules.simulation as simulation
import modules.solver as solver
from modules.balance import ventilation_rate_is
from modules.newton import NewtonNotConverged


def section_8_2():
    """8.2 節：換気量を、ステップの終わりの温度からニュートン法の中で計算し直すと、どこで止まるか。"""
    h0, m0, step0 = solver.heat_residual_is, solver.moisture_residual_is, simulation.step
    count = [0]

    def h_res(geometry, t_npls, mu_npls, t_n, dt, outdoor_npls, indoor_npls, q_sol, q_v, *a):
        return h0(geometry, t_npls, mu_npls, t_n, dt, outdoor_npls, indoor_npls, q_sol,
                  ventilation_rate_is(geometry, t_npls, outdoor_npls), *a)

    def m_res(geometry, t_npls, mu_npls, t_n, mu_n, dt, outdoor_npls, indoor_npls, q_v, *a):
        return m0(geometry, t_npls, mu_npls, t_n, mu_n, dt, outdoor_npls, indoor_npls,
                  ventilation_rate_is(geometry, t_npls, outdoor_npls), *a)

    def step(geometry, state_n, dt, outdoor_n, outdoor_npls, indoor_npls, *a):
        count[0] += 1
        try:
            return step0(geometry, state_n, dt, outdoor_n, outdoor_npls, indoor_npls, *a)
        except NewtonNotConverged:
            i_air = next(i for i, n in enumerate(geometry.nodes) if type(n).__name__ == 'AirNode')
            print(f'8.2 節　{count[0]} ステップ目で止まった：ステップの初めの通気層の空気 {state_n.t_is[i_air] - 273.15:.2f} deg.C、'
                  f'外気 {outdoor_n.state.t - 273.15:.2f} deg.C、ステップの終わりの外気 {outdoor_npls.state.t - 273.15:.2f} deg.C')
            raise

    try:
        with c.patched((solver, 'heat_residual_is', h_res), (solver, 'moisture_residual_is', m_res), (simulation, 'step', step)):
            c.ex.run(3600.0, 2)
        print('8.2 節　最後まで進んだ')
    except NewtonNotConverged:
        pass
    it = c.run_example('base')['iterations_ns']
    print(f'  今の方法（換気量をステップの中で一定）：修正の回数 平均 {it.mean():.2f} 回、最大 {it.max()} 回')


def table_1():
    """表 1：壁 1 m2 あたりの換気量の例（流量係数 0.085、通気層の厚さ 25 mm、縦の長さと高さの差 7 m、外気 0 deg.C）。"""
    from modules.wall_geometry import AirNode
    node = AirNode(thickness=0.025, alpha=0.085, layer_index=1)
    print('表 1　壁 1 m2 あたりの換気量の例')
    for d in (1.0, 5.0, 10.0):
        q = node.ventilation_rate(273.15 + d, 273.15, 7.0, 7.0)
        print(f'  温度の差 {d:.0f} K：{q:.5f} m3/(m2 s)、流速 {q * 7.0 / 0.025:.2f} m/s、換気回数 {q / 0.025 * 3600.0:.0f} 回/h')


def main():
    table_1()
    section_8_2()


if __name__ == '__main__':
    main()
