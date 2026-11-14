"""04_05「数値解法の安定性」：表 1（許容値）、表 2（修正の行き過ぎを抑える仕組み）、表 3（修正の回数の下限）、6.2 節（修正の回数）。"""

import numpy as np

import _common as c
import modules.newton as newton
import modules.solver as solver
from modules.air_water_properties import L_V
from modules.boundary import IndoorBoundary, OutdoorBoundary
from modules.wall_builder import build_wall
from modules.wall_input import WallInput

_SURFACE = {'direction': 'e', 'angle': 90.0, 'albedo': 0.3, 'absorptance': 0.9, 'emissivity': 0.9, 'eva_height': 4.2,
            'exposure_factor': 1.0, 'deposition_factor': 1.0}


def table_1():
    """表 1：許容値を変えた場合の、今の許容値との差と修正の回数。"""
    base = c.run_example('base')
    y2 = c.year2(base)
    rh_b = c.rh_is_ns(base['t_is_ns'][:, y2], base['mu_is_ns'][:, y2])
    print('表 1　許容値を変えた場合（年間の計算、2 年目）')
    for eps in (1.0e-4, 1.0e-6, 1.0e-8):
        d = base if eps == 1.0e-6 else c.run_example(
            f'tolerance_{eps:.0e}', patches=((solver, 'HEAT_RESIDUAL_TOL', eps), (solver, 'MOISTURE_RESIDUAL_TOL', eps / L_V)))
        it = d['iterations_ns']
        if eps == 1.0e-6:
            print(f'  {eps:.0e} W/m2（今の値）：修正の回数 平均 {it.mean():.2f} 回・最大 {it.max()} 回')
            continue
        d_t = np.abs(d['t_is_ns'][:, y2] - base['t_is_ns'][:, y2])
        d_rh = np.abs(c.rh_is_ns(d['t_is_ns'][:, y2], d['mu_is_ns'][:, y2]) - rh_b)
        print(f'  {eps:.0e} W/m2：温度の差 {d_t.max():.1e} K・{d_t.mean():.1e} K、相対湿度の差 {d_rh.max():.1e}%・{d_rh.mean():.1e}%、'
              f'修正の回数 平均 {it.mean():.2f} 回・最大 {it.max()} 回')


def section_6_2():
    """6.2 節：年間の計算の修正の回数の分布。"""
    it = c.run_example('base')['iterations_ns']
    counts = np.bincount(it)
    print('6.2 節　年間の計算の修正の回数（ステップ数）：' + '、'.join(f'{m} 回が {n}' for m, n in enumerate(counts) if n))


def _brick(name: str):
    """れんが 1 層（厚さ 20 mm、分割数 20、20 deg.C・50%）。"""
    return build_wall(WallInput.from_dict({
        'surface': _SURFACE,
        'layers': [{'source': 'mews_ashrae', 'name': name, 'thickness': 0.02, 'n_div': 20,
                    'initial_temperature': 20.0, 'initial_rh': 50.0}],
        'interfaces': [{'g_h': 22.4, 'g_m': 3.43e-8}, {'g_h': 9.2, 'g_m': 2.6e-10}],
        'len_vertical': 1.0, 'len_horizontal': 1.0}))


def table_2():
    """表 2：乾いたれんがの壁に湿った外気を当てたときの、仕組みの有無と計算の結果。"""
    solve = newton.newton_solve

    def variant(stop: bool, line_search: bool):
        def f(res, x0, band, tol, max_iter=20, min_iter=0, x_max_ps=None, x_stop_ps=None):
            with c.patched((newton, 'MAX_HALVING', newton.MAX_HALVING if line_search else 0)):
                return solve(res, x0, band, tol, max_iter, min_iter, x_max_ps, x_stop_ps if stop else None)
        return f

    variants = (('ふつうのニュートン法', False, False), ('直線探索だけ', False, True),
                ('またがせない値だけ', True, False), ('またがせない値と直線探索（今の方法）', True, True))
    names = ('ASHRAE Clay Brick', 'MEWS Brick II - Red Matt Clay Brick', 'MEWS Brick VI - Calcium Silicate Brick')
    indoor = IndoorBoundary.from_theta_rh(theta=20.0, rh=50.0)
    print('表 2　乾いたれんがの壁に湿った外気を当てたとき（12 の条件のうち、最後まで進んだ数）')
    for label, stop, line_search in variants:
        ok, max_it = 0, 0
        with c.patched((solver, 'newton_solve', variant(stop, line_search))):
            for dt in (3600.0, 86400.0):
                for name in names:
                    for rh in (97.0, 99.0):
                        outdoor = OutdoorBoundary.from_theta_rh(theta=20.0, rh=rh)
                        g, st = _brick(name)
                        try:
                            for _ in range(int(30 * 86400 / dt)):
                                st, m = solver.step(g, st, dt, outdoor, outdoor, indoor)
                                max_it = max(max_it, m) if label.endswith('（今の方法）') else max_it
                            ok += 1
                        except newton.NewtonNotConverged:
                            pass
        print(f'  {label}：{ok} / 12' + (f'（修正の回数の最大 {max_it} 回）' if label.endswith('（今の方法）') else ''))


def table_3():
    """表 3：修正の回数の下限と、1 年後の相対湿度（構造用合板 1 層、相対湿度 50% から 70% への吸湿）。"""
    solve = newton.newton_solve
    outdoor = OutdoorBoundary.from_theta_rh(theta=20.0, rh=70.0)
    indoor = IndoorBoundary.from_theta_rh(theta=20.0, rh=70.0)
    print('表 3　修正の回数の下限と、1 年後の相対湿度')
    for min_iter in (0, 1):
        def f(res, x0, band, tol, max_iter=20, min_iter_given=None, x_max_ps=None, x_stop_ps=None, **kw):
            return solve(res, x0, band, tol, max_iter, min_iter, x_max_ps, x_stop_ps)
        g, st = build_wall(WallInput.from_dict({
            'surface': _SURFACE,
            'layers': [{'source': 'whamcavity', 'name': '構造用合板1', 'thickness': 0.009, 'n_div': 6,
                        'initial_temperature': 20.0, 'initial_rh': 50.0}],
            'interfaces': [{'g_h': 22.4, 'g_m': 2.0e-11}, {'g_h': 9.2, 'g_m': 2.6e-10}],
            'len_vertical': 1.0, 'len_horizontal': 1.0}))
        zero, first = 0, None
        with c.patched((solver, 'newton_solve', f)):
            for day in range(365):
                st, m = solver.step(g, st, 86400.0, outdoor, outdoor, indoor)
                if m == 0:
                    zero += 1
                    first = first or day + 1
        rh = st.rh_is
        print(f'  下限 {min_iter} 回：1 年後の相対湿度 {rh.min():.5f}%〜{rh.max():.5f}%、'
              f'修正せずに進んだステップ {zero}' + (f'（{first} 日目から）' if first else ''))


def main():
    table_1()
    section_6_2()
    table_2()
    table_3()


if __name__ == '__main__':
    main()
