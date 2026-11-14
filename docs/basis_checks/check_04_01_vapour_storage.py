"""04_01「支配方程式と離散化」5.3.3 節：空隙にたまる水蒸気の量の変化を省く近似の影響（年間の計算）。

近似をしない計算として、水分の収支と熱の収支に、ステップの初めと終わりの空隙の水蒸気の量の差の項を足す。
  gw ：グラスウールの質点だけ、空隙率 0.99 で足す
  all：すべての材料の質点に、空隙率 1 で足す（上限の見積もり）
図 1（docs/source/contents/figures/vapour_storage_rh_diff.svg）も描き直す。
"""

import os

import numpy as np

import _common as c
import modules.simulation as simulation
import modules.solver as solver
from modules.air_water_properties import L_V, R_V
from modules.state import State
from modules.wall_geometry import FilmNode, MaterialNode

_FIGURE = os.path.join(c.ROOT, 'docs', 'source', 'contents', 'figures', 'vapour_storage_rh_diff.svg')


def _patches(case: str) -> tuple:
    """近似をしない計算の、差し替える部品。"""

    def psi_width_is(geometry):
        """空隙率 × 幅, m, [i]。"""
        return np.array([
            0.0 if not isinstance(n, MaterialNode)
            else (1.0 if case == 'all' else 0.99 if 'グラスウール' in n.material.name else 0.0) * n.width
            for n in geometry.nodes])

    def vapour_is(geometry, t_is, mu_is):
        """空隙の水蒸気の量, kg/m2, [i]。"""
        rho_v = np.array([State(t=t, mu=min(m, 0.0)).p_v for t, m in zip(t_is, mu_is)]) / (R_V * np.asarray(t_is))
        return psi_width_is(geometry) * rho_v

    saved = {}
    step0, m0, h0 = solver.step, solver.moisture_residual_is, solver.heat_residual_is

    def step(geometry, state_n, dt, *a, **k):
        saved['v_n'] = vapour_is(geometry, state_n.t_is, state_n.mu_is)
        return step0(geometry, state_n, dt, *a, **k)

    def m_res(geometry, t_npls, mu_npls, t_n, mu_n, dt, *a, **k):
        return m0(geometry, t_npls, mu_npls, t_n, mu_n, dt, *a, **k) + (vapour_is(geometry, t_npls, mu_npls) - saved['v_n']) / dt

    def h_res(geometry, t_npls, mu_npls, t_n, dt, *a, **k):
        # 空隙にたまる水蒸気の分は相変化しないので、凝縮熱から引く
        return h0(geometry, t_npls, mu_npls, t_n, dt, *a, **k) + L_V * (vapour_is(geometry, t_npls, mu_npls) - saved['v_n']) / dt

    return ((simulation, 'step', step), (solver, 'moisture_residual_is', m_res), (solver, 'heat_residual_is', h_res))


def _figure(geometry, rh_diff_max: dict):
    """図 1：質点ごとの相対湿度の差の最大。"""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    keep = [i for i, n in enumerate(geometry.nodes) if not isinstance(n, FilmNode)]
    x_is = c.ex._node_positions(geometry)[keep] * 1000
    text, text2, grid, surf = '#0b0b0b', '#52514e', '#e4e3df', '#fcfcfb'
    plt.rcParams.update({'font.family': 'Hiragino Sans', 'font.size': 10, 'axes.edgecolor': text2, 'axes.labelcolor': text,
                         'xtick.color': text2, 'ytick.color': text2, 'svg.fonttype': 'path', 'svg.hashsalt': 'vapour_storage'})
    fig, ax = plt.subplots(figsize=(9, 4.2), facecolor=surf)
    ax.set_facecolor(surf)
    for case, label, color, marker in (('gw', 'グラスウールの層だけ近似をしない（空隙率 0.99）', '#2a78d6', 'o'),
                                       ('all', 'すべての材料の層で近似をしない（空隙率 1）', '#eb6834', 's')):
        ax.plot(x_is, rh_diff_max[case][keep], color=color, linewidth=2, marker=marker, markersize=5, label=label,
                zorder=3, markeredgecolor=surf, markeredgewidth=1)
    edges = np.cumsum([0.0] + [n.width for n in geometry.nodes]) * 1000
    bounds = [edges[0]] + [edges[i] for i in range(1, len(geometry.nodes))
                           if geometry.nodes[i].layer_index != geometry.nodes[i - 1].layer_index] + [edges[-1]]
    for b in bounds[1:-1]:
        ax.axvline(b, color=text2, linewidth=0.6, linestyle=':')
    names = ['サイディング', '通気層', 'せっこう\nボード', '合板', 'グラスウール', '合板', 'せっこう\nボード']
    ax.set_ylim(0, 1.0)
    for a, b, nm in zip(bounds[:-1], bounds[1:], names):
        ax.text((a + b) / 2, 0.97, nm, ha='center', va='top', fontsize=8, color=text2)
    ax.set_xlim(edges[0], edges[-1])
    ax.set_xlabel('外表面からの距離, mm（点線は層の境目）')
    ax.set_ylabel('相対湿度の差の最大, %')
    ax.grid(True, axis='y', color=grid, linewidth=0.8)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    ax.legend(loc='upper left', bbox_to_anchor=(0.0, 0.86), frameon=False, fontsize=9, labelcolor=text)
    fig.tight_layout()
    # 書き出すたびにファイルが変わらないように、日時を書かない（図の部品の名前は svg.hashsalt で固定する）
    fig.savefig(_FIGURE, metadata={'Date': None})
    print(f'図 1 を書き出しました：{os.path.relpath(_FIGURE, c.ROOT)}')


def main():
    geometry = c.example_geometry()
    base = c.run_example('base')
    y2 = c.year2(base)
    rh_base = c.rh_is_ns(base['t_is_ns'][:, y2], base['mu_is_ns'][:, y2])
    material = np.array([isinstance(n, MaterialNode) for n in geometry.nodes])
    i_gw_in = max(i for i, n in enumerate(geometry.nodes) if isinstance(n, MaterialNode) and 'グラスウール' in n.material.name)

    rh_diff_max = {}
    for case in ('gw', 'all'):
        d = c.run_example(f'vapour_storage_{case}', patches=_patches(case))
        rh = c.rh_is_ns(d['t_is_ns'][:, y2], d['mu_is_ns'][:, y2])
        d_rh = np.abs(rh - rh_base)
        d_t = np.abs(d['t_is_ns'][:, y2] - base['t_is_ns'][:, y2])
        rh_diff_max[case] = d_rh.max(axis=1)
        print(f'[{case}]')
        print(f'  温度の差（すべての質点）：最大 {d_t.max():.3f} K、平均 {d_t.mean():.4f} K')
        print(f'  相対湿度の差の最大：{d_rh.max():.3f}%（質点 {d_rh.max(axis=1).argmax()}）')
        print(f'  年平均の相対湿度の差の最大（材料の質点）：{np.abs(rh.mean(axis=1) - rh_base.mean(axis=1))[material].max():.3f}%')
        print(f'  グラスウールの室内側の端（質点 {i_gw_in}）の 2 年目の最高：近似あり {rh_base[i_gw_in].max():.1f}%、'
              f'近似なし {rh[i_gw_in].max():.1f}%')
        print(f'  修正の回数：平均 {d["iterations_ns"].mean():.2f} 回（近似あり {base["iterations_ns"].mean():.2f} 回）')
    _figure(geometry, rh_diff_max)


if __name__ == '__main__':
    main()
