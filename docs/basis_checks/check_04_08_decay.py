"""04_08「木材腐朽のモデル」：表 1〜表 5。

表 1：途中で条件が変わる場合の、発芽までの経過（案 A と案 B）
表 2：発芽までの時間の式の値 t'_i
表 3：反応速度定数の式の値
表 4：腐朽で生まれる水の影響（厚さ 50 mm の木材の単板、20 deg.C・98.5% 一定、1 年）
表 5：判定に使う値（1 日の平均か、時間刻みごとか）と、1 年後の表面の質量減少率（相対湿度が 1 日の周期で上下する単板）
表 6：係数 WS の値と、発芽までの時間の式の値 t'_i
図 1：反応速度定数の式と、齋藤ら（2008）の図 8 から読み取った実験の値（docs/source/contents/figures/decay_reaction_rate.svg）
"""

import os

import numpy as np
from scipy.optimize import brentq

import _common as c
from modules import decay
from modules.boundary import Boundaries, IndoorBoundary, OutdoorBoundary
from modules.material_loaders import load_materials
from modules.simulation import simulate
from modules.wall_builder import build_wall
from modules.wall_input import WallInput


def _f_g(theta: float, rh: float, ws: float = decay.WS) -> tuple[float, float]:
    """発芽までの時間の式の係数 f と g（「木材腐朽」の式 (2)・(3)）。"""
    f = 0.1384 * theta + 0.4370 * rh - 42.9450 + ws * (0.0340 * theta - 0.0210 * rh + 1.7210)
    g = -2.2270 * theta - 0.0347 * rh + 0.0244 * theta * rh + ws * (-0.504 * theta + 0.0096 * rh + 0.0047 * theta * rh)
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
        'surface': {'direction': 'e', 'angle': 90.0, 'albedo': 0.3, 'absorptance': 0.9, 'emissivity': 0.9,
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


def table_6():
    print('表 6　係数 WS の値と、発芽までの時間の式の値 t\'_i（月）')
    for theta, rh in ((5.0, 98.0), (10.0, 98.0), (15.0, 95.0), (20.0, 95.0), (20.0, 98.0)):
        cells = []
        for ws in (-1.0, 0.0, 1.0):
            f, g = _f_g(theta, rh, ws)
            cells.append(f'WS = {ws:.0f}：{-g / f:.2f}')
        print(f'  {theta:.0f} deg.C・{rh:.0f}%：' + '、'.join(cells))


# 齋藤ら（2008）の図 8 から目で読み取った、反応速度定数の実験の値（温度 deg.C、1e-10 1/s）
_KM_READ = ((0.0, 0.0), (5.0, 7.0), (10.0, 44.0), (20.0, 130.0), (30.0, 153.0), (40.0, 30.0))

_FIGURE = os.path.join(c.ROOT, 'docs', 'source', 'contents', 'figures', 'decay_reaction_rate.svg')


def figure_1():
    """図 1：反応速度定数の式（3 次式と、負を 0 で切ったもの）と、実験の値の読み取り。0〜6 deg.C の拡大を右上に入れる。"""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    text, text2, grid, surf = '#0b0b0b', '#52514e', '#e4e3df', '#fcfcfb'
    blue, orange, gray = '#2a78d6', '#eb6834', '#9a9893'
    plt.rcParams.update({'font.family': 'Hiragino Sans', 'font.size': 10, 'axes.edgecolor': text2, 'axes.labelcolor': text,
                         'xtick.color': text2, 'ytick.color': text2, 'svg.fonttype': 'path', 'svg.hashsalt': 'decay_km'})
    cubic = lambda t: 2.77 - 3.23 * t + 0.865 * t ** 2 - 0.0189 * t ** 3
    theta = np.linspace(0.0, 45.0, 901)
    fig, ax = plt.subplots(figsize=(9, 4.4), facecolor=surf)
    ax.set_facecolor(surf)
    ax.axvspan(40.0, 45.0, color=grid, alpha=0.6, linewidth=0)
    ax.text(42.5, 185, '腐朽が進まない\n（40 deg.C を超える）', ha='center', va='top', fontsize=8, color=text2)
    ax.axhline(0.0, color=text2, linewidth=0.8)
    ax.plot(theta, cubic(theta), color=blue, linewidth=2, label='式 (18) の 3 次式', zorder=3)
    ax.plot(*zip(*_KM_READ), linestyle='none', marker='o', markersize=7, color=orange, markeredgecolor=surf,
            markeredgewidth=1.2, label='実験の値（齋藤らの図 8 からの読み取り）', zorder=4)
    ax.set_xlim(0.0, 45.0)
    ax.set_ylim(-20.0, 190.0)
    ax.set_xlabel('温度, deg.C')
    ax.set_ylabel('反応速度定数, 10$^{-10}$ 1/s')
    ax.grid(True, color=grid, linewidth=0.8)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    ax.legend(loc='lower center', bbox_to_anchor=(0.62, 0.16), frameon=False, fontsize=9, labelcolor=text)

    # 0〜6 deg.C の拡大：3 次式が 1.29〜2.72 deg.C で負になり、そこを 0 で切る
    inset = ax.inset_axes([0.07, 0.6, 0.24, 0.34])
    inset.set_facecolor(surf)
    small = np.linspace(0.0, 6.0, 601)
    roots = sorted(r.real for r in np.roots([-0.0189, 0.865, -3.23, 2.77]) if abs(r.imag) < 1e-12 and r.real < 10)
    inset.axvspan(roots[0], roots[1], color=grid, alpha=0.9, linewidth=0)
    inset.axhline(0.0, color=text2, linewidth=0.8)
    inset.plot(small, cubic(small), color=blue, linewidth=1.2, linestyle='--')
    inset.plot(small, np.maximum(cubic(small), 0.0), color=blue, linewidth=2)
    inset.plot(*zip(*[p for p in _KM_READ if p[0] <= 6.0]), linestyle='none', marker='o', markersize=5, color=orange,
               markeredgecolor=surf, markeredgewidth=1)
    inset.set_xlim(0.0, 6.0)
    inset.set_ylim(-1.0, 8.0)
    inset.tick_params(labelsize=7)
    inset.set_title('0〜6 deg.C の拡大', fontsize=8, color=text2, loc='left')
    print(f'  式が負になる温度：{roots[0]:.2f}〜{roots[1]:.2f} deg.C')
    for side in ('top', 'right'):
        inset.spines[side].set_visible(False)
    fig.tight_layout()
    fig.savefig(_FIGURE, metadata={'Date': None})
    print(f'図 1 を書き出しました：{os.path.relpath(_FIGURE, c.ROOT)}')


def main():
    figure_1()
    tables_1_to_3()
    table_6()
    table_4()
    table_5()


if __name__ == '__main__':
    main()
