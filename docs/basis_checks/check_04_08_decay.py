"""04_08「木材腐朽のモデル」：表 1〜表 5。

表 1：途中で条件が変わる場合の、発芽までの経過（案 A と案 B）
表 2：発芽までの時間の式の値 t'_i
表 3：反応速度定数の式の値
表 4：腐朽で生まれる水の影響（厚さ 50 mm の木材の単板、20 deg.C・98.5% 一定、1 年）
表 5：判定に使う値（1 日の平均か、時間刻みごとか）と、1 年後の表面の質量減少率（相対湿度が 1 日の周期で上下する単板）
"""

import numpy as np
from scipy.optimize import brentq

import _common as c
from modules import decay
from modules.boundary import Boundaries, IndoorBoundary, OutdoorBoundary
from modules.material_loaders import load_materials
from modules.simulation import simulate
from modules.wall_builder import build_wall
from modules.wall_input import WallInput


def _f_g(theta: float, rh: float) -> tuple[float, float]:
    """発芽までの時間の式の係数 f と g（「木材腐朽」の式 (2)・(3)）。"""
    f = 0.1384 * theta + 0.4370 * rh - 42.9450 + decay.WS * (0.0340 * theta - 0.0210 * rh + 1.7210)
    g = -2.2270 * theta - 0.0347 * rh + 0.0244 * theta * rh + decay.WS * (-0.504 * theta + 0.0096 * rh + 0.0047 * theta * rh)
    return f, g


def tables_1_to_3():
    days = lambda theta, rh: decay.initial_response_time(theta, rh) / 86400.0
    t1, t2 = days(20.0, 95.0), days(15.0, 98.0)
    print(f'表 1　20 deg.C・95% の発芽までの時間 {t1:.1f} 日、15 deg.C・98% は {t2:.1f} 日')
    print(f'  案 A：{t2:.1f} 日で発芽。案 B：10 日で {10 / t1:.2f}、発芽は {10 + (1 - 10 / t1) * t2:.1f} 日')
    t93, t98 = days(20.0, 93.0), days(20.0, 98.0)
    print(f'  4.3.2 節の例：20 deg.C・93% は {t93 / 30:.0f} か月、98% は {t98:.1f} 日。案 B で 300 日の進み具合 {300 / t93:.2f}、'
          f'98% になってから {(1 - 300 / t93) * t98:.1f} 日で発芽')
    print('表 2　発芽までの時間の式の値 t\'_i（月）')
    for theta in (5.0, 10.0, 15.0, 20.0, 25.0, 30.0):
        cells = []
        for rh in (93.0, 95.0, 98.0, 100.0):
            f, g = _f_g(theta, rh)
            cells.append(f'{-g / f:.2f}' + ('（f < 0）' if f < 0 else '（g > 0）' if g > 0 else ''))
        print(f'  {theta:.0f} deg.C：' + ' | '.join(cells))
    print('表 3　反応速度定数の式の値（1e-10 1/s。負の値も示す）')
    k = lambda theta: (2.77 - 3.23 * theta + 0.865 * theta ** 2 - 0.0189 * theta ** 3)
    print('  ' + '、'.join(f'{theta:.0f} deg.C：{k(theta):.2f}' for theta in (0, 1, 2, 3, 5, 10, 20, 30, 40)))
    roots = sorted(r.real for r in np.roots([-0.0189, 0.865, -3.23, 2.77]) if abs(r.imag) < 1e-12)
    print(f'  式が 0 になる温度：' + '、'.join(f'{r:.2f}' for r in roots) + ' deg.C')


def _panel(rh_at, wood_decay: bool):
    """厚さ 50 mm の木材の単板を、時間刻み 1 時間で 1 年計算する（片側はほぼ断湿、もう片側は 20 deg.C・相対湿度 rh_at(時刻)）。"""
    m = next(x for x in load_materials('whamcavity') if x.name == '木材')
    rh_0 = brentq(lambda r: m.get_u(r) - 0.2, 50.0, 99.99)
    g, st = build_wall(WallInput.from_dict({
        'surface': {'direction': 'e', 'angle': 90.0, 'albedo': 0.3, 'absorptance': 0.9, 'emissivity': 0.9, 'eva_height': 4.2,
                    'exposure_factor': 1.0, 'deposition_factor': 1.0},
        'layers': [{'source': 'whamcavity', 'name': '木材', 'thickness': 0.05, 'n_div': 51, 'initial_temperature': 20.0,
                    'initial_rh': rh_0, 'wood_decay': wood_decay}],
        'interfaces': [{'g_h': 9.3, 'g_m': 1.0e-16}, {'g_h': 9.3, 'g_m': 3.43e-8}],
        'len_vertical': 1.0, 'len_horizontal': 1.0}))
    outdoor = OutdoorBoundary.from_theta_rh(theta=20.0, rh=50.0)
    r = simulate(g, st, 3600.0, 8760,
                 lambda t: Boundaries(outdoor=outdoor, indoor=IndoorBoundary.from_theta_rh(theta=20.0, rh=rh_at(t))))
    return g, r


def _post(g, r, daily: bool) -> np.ndarray:
    """腐朽を求めない計算の結果から、腐朽の状態だけを進めた質量減少率（腐朽で生まれる水なし）, [i]。"""
    nodes = decay.DecayNodes(
        is_target_is=np.ones(len(g.nodes), dtype=bool),
        is_surface_is=np.array([i in (0, len(g.nodes) - 1) for i in range(len(g.nodes))]),
        i_left_is=np.arange(len(g.nodes)) - 1,
        i_right_is=np.array([i + 1 if i + 1 < len(g.nodes) else -1 for i in range(len(g.nodes))]),
        width_is=np.array([n.width for n in g.nodes]))
    theta = r.t_is_ns[:, 1:] - 273.15
    rh = c.rh_is_ns(r.t_is_ns[:, 1:], r.mu_is_ns[:, 1:])
    state = decay.DecayState.initial(len(g.nodes))
    if daily:
        for d in range(365):
            sl = slice(24 * d, 24 * (d + 1))
            state = decay.advance(state, nodes, theta[:, sl].mean(axis=1), rh[:, sl].mean(axis=1), 86400.0)
    else:
        for n in range(theta.shape[1]):
            state = decay.advance(state, nodes, theta[:, n], rh[:, n], 3600.0)
    return state.mass_loss_is


def _moisture_content(g, r) -> np.ndarray:
    """1 年後の含水率, kg/kg, [i]。"""
    return np.array([n.water_amount(r.state_last.state(i)) / n.width / n.material.rho_s for i, n in enumerate(g.nodes)])


def table_4():
    print('表 4　腐朽で生まれる水の影響（1 年後）')
    g, off = _panel(lambda t: 98.5, wood_decay=False)
    _, on = _panel(lambda t: 98.5, wood_decay=True)
    for label, mass_loss, u in (('入れない', _post(g, off, daily=False), _moisture_content(g, off)),
                                ('入れる', on.mass_loss_is_ns[:, -1], _moisture_content(g, on))):
        print(f'  {label}：表面の質量減少率 {mass_loss[-1]:.3f}、腐朽が進んだ質点 {(mass_loss > 0).sum()}、'
              f'含水率の最大 {u.max():.3f} kg/kg、表面から 10 mm {u[-11]:.3f} kg/kg')


def table_5():
    print('表 5　判定に使う値と、1 年後の表面の質点の質量減少率')
    for mean, amp in ((95.0, 4.5), (96.5, 3.5)):
        rh_at = lambda t, mean=mean, amp=amp: mean + amp * np.sin(2 * np.pi * t / 86400.0)
        g, off = _panel(rh_at, wood_decay=False)
        _, on = _panel(rh_at, wood_decay=True)
        rh_s = c.rh_is_ns(off.t_is_ns[-1, -24 * 5:], off.mu_is_ns[-1, -24 * 5:])
        print(f'  {mean} ± {amp}%（表面の質点 {rh_s.min():.1f}〜{rh_s.max():.1f}%、平均 {rh_s.mean():.1f}%）：'
              f'1 日の平均で判定 {_post(g, off, daily=True)[-1]:.3f}、時間刻みごとに判定 {on.mass_loss_is_ns[-1, -1]:.3f}')


def main():
    tables_1_to_3()
    table_4()
    table_5()


if __name__ == '__main__':
    main()
