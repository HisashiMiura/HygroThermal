"""例の入力のファイル examples/main_wall.json（Fortran の入力の壁 7 層）を、時間刻みと年数だけ変えて計算し、結果の図を描く。

使い方（リポジトリの一番上のフォルダで）：
    python examples/one_year_main_wall.py            # Δt = 3600 秒
    python examples/one_year_main_wall.py --dt 900   # Δt = 900 秒（Δt で答えが変わらないかを比べる）

結果の数値（npz）と図（png）は examples/output/ に書き出す（git には入れない）。
壁・外気・室内は入力のファイルのとおりで、時間刻みは --dt、年数は --years で変える（記録は 1 時間ごと）。
通気層の換気、境目の水膜、材料の中の自由水、木材腐朽（構造用合板の 2 層）を含む。雨水を浸入させる境目はない。
"""

import argparse
import os
import sys
import time

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from modules.boundary import Boundaries
from modules.indoor_climate import CosineIndoor
from modules.outdoor_climate import WeatherOutdoor
from modules.output import node_positions, rh_is_ns
from modules.simulation import SimulationResult, simulate
from modules.simulation_input import SimulationInput
from modules.air_water_properties import T_0
from modules.time_axis import SECONDS_PER_DAY, SECONDS_PER_HOUR, SECONDS_PER_YEAR
from modules.wall_builder import build_wall
from modules.wall_geometry import WallGeometry
from modules.wall_input import WallInput

_HERE = os.path.dirname(__file__)
_OUTPUT = os.path.join(_HERE, 'output')

# 例の入力のファイル
_INPUT = os.path.join(_HERE, 'main_wall.json')

# 図の色（データの図の決まりに沿った、色覚の違いでも見分けやすい順の 5 色）と、文字・線の色
_SERIES = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4']
_TEXT = '#0b0b0b'
_TEXT_2 = '#52514e'
_GRID = '#e4e3df'
_SURFACE = '#fcfcfb'


def _main_wall_input() -> WallInput:
    """例の入力のファイルの壁（Fortran の入力の壁 7 層。構造用合板の 2 層で木材腐朽を求める）。"""
    return SimulationInput.from_json(_INPUT).wall


def _daily_mean(v_ns: np.ndarray, per_day: int) -> np.ndarray:
    """1 時間ごとの値を、1 日の平均にする（最初の記録（時刻 0）は除く）。"""
    v = v_ns[1:]
    return v[: len(v) // per_day * per_day].reshape(-1, per_day).mean(axis=1)


def _spread(y_ls: list[float], gap: float) -> list[float]:
    """線の右端に付ける名前の高さを、重ならないように gap 以上ずつ離す（順番は変えない）。"""
    order = np.argsort(y_ls)
    y_new = np.array(y_ls, dtype=float)
    for a, b in zip(order[:-1], order[1:]):
        y_new[b] = max(y_new[b], y_new[a] + gap)
    return list(y_new)


def run(dt: float, years: int) -> tuple[WallGeometry, SimulationResult, WeatherOutdoor, CosineIndoor]:
    """壁を作って、years 年間計算する。"""

    geometry, state_0 = build_wall(_main_wall_input())
    sim_input = SimulationInput.from_json(_INPUT)
    weather = sim_input.outdoor.load()
    indoor = sim_input.indoor

    def boundaries_at(time_s: float) -> Boundaries:
        return Boundaries(outdoor=weather.get_outdoor_boundary(time_s), indoor=indoor.get_indoor_boundary(time_s))

    n_steps = int(round(years * SECONDS_PER_YEAR / dt))
    record_every = int(round(SECONDS_PER_HOUR / dt))
    result = simulate(geometry, state_0, dt=dt, n_steps=n_steps, boundaries_at=boundaries_at, record_every=record_every)
    return geometry, result, weather, indoor


def plot(geometry: WallGeometry, result: SimulationResult, weather: WeatherOutdoor, indoor: CosineIndoor, tag: str):
    """図を 3 枚描く：2 年目の温度・相対湿度の 1 日の平均の変化、質点ごとの 2 年目の相対湿度の範囲。"""

    plt.rcParams.update({'font.family': 'Hiragino Sans', 'font.size': 10, 'axes.edgecolor': _TEXT_2,
                         'axes.labelcolor': _TEXT, 'xtick.color': _TEXT_2, 'ytick.color': _TEXT_2})

    # 2 年目の記録（1 時間ごと）
    time_ns = result.time_ns
    year2 = (time_ns > SECONDS_PER_YEAR) & (time_ns <= 2 * SECONDS_PER_YEAR)
    time_2 = np.concatenate(([SECONDS_PER_YEAR], time_ns[year2]))
    sel = np.concatenate(([np.argmax(year2) - 1], np.flatnonzero(year2)))
    t_is = result.t_is_ns[:, sel]
    rh_is = rh_is_ns(t_is, result.mu_is_ns[:, sel])

    # 外気と室内（記録と同じ時刻）
    outdoor_ns = [weather.get_outdoor_boundary(t) for t in time_2]
    indoor_ns = [indoor.get_indoor_boundary(t) for t in time_2]

    # 図にする場所（外気・通気層・外側の構造用合板のグラスウール側・内側の構造用合板の外側（防湿層の室内側）・室内）
    places = [
        ('外気', np.array([b.state.theta for b in outdoor_ns]), np.array([b.state.rh for b in outdoor_ns])),
        ('通気層（質点 6）', t_is[6] - T_0, rh_is[6]),
        ('外側の合板・GW 側（質点 19）', t_is[19] - T_0, rh_is[19]),
        ('内側の合板・外側（質点 27）', t_is[27] - T_0, rh_is[27]),
        ('室内', np.array([b.state.theta for b in indoor_ns]), np.array([b.state.rh for b in indoor_ns])),
    ]
    per_day = int(round(SECONDS_PER_DAY / SECONDS_PER_HOUR))
    day_ds = np.arange(1, 366)

    for k, (title, unit, file) in enumerate((('温度（2 年目、1 日の平均）', 'deg.C', 'temperature'),
                                              ('相対湿度（2 年目、1 日の平均）', '%', 'rh'))):
        fig, ax = plt.subplots(figsize=(10, 5), facecolor=_SURFACE)
        ax.set_facecolor(_SURFACE)
        v_ds_list = [_daily_mean(theta_ns if k == 0 else rh_ns, per_day) for _, theta_ns, rh_ns in places]
        for (name, _, _), v_ds, color in zip(places, v_ds_list, _SERIES):
            ax.plot(day_ds[: len(v_ds)], v_ds, color=color, linewidth=1.5, label=name)

        # 線の右端に名前を付ける（重ならないように、縦の幅の 4% ずつ離す）
        y_min, y_max = ax.get_ylim()
        y_labels = _spread([v_ds[-1] for v_ds in v_ds_list], gap=0.04 * (y_max - y_min))
        for (name, _, _), v_ds, y in zip(places, v_ds_list, y_labels):
            ax.annotate(name, (day_ds[len(v_ds) - 1], v_ds[-1]), xytext=(day_ds[len(v_ds) - 1] + 4, y),
                        textcoords='data', va='center', fontsize=8, color=_TEXT_2)
        ax.set_title(f'{title}　Δt = {tag} 秒', loc='left', color=_TEXT)
        ax.set_xlabel('通算日（1 月 1 日が 1）')
        ax.set_ylabel(unit)
        ax.set_xlim(1, 420)
        ax.grid(True, color=_GRID, linewidth=0.8)
        for s in ('top', 'right'):
            ax.spines[s].set_visible(False)
        ax.legend(loc='upper left', frameon=False, fontsize=8, labelcolor=_TEXT)
        fig.tight_layout()
        fig.savefig(os.path.join(_OUTPUT, f'{file}_dt{tag}.png'), dpi=150)
        plt.close(fig)

    # 質点ごとの 2 年目の相対湿度（最低・平均・最高）を、壁の厚さ方向に並べる
    x_is = node_positions(geometry) * 1000.0
    fig, ax = plt.subplots(figsize=(10, 5), facecolor=_SURFACE)
    ax.set_facecolor(_SURFACE)
    ax.fill_between(x_is, rh_is.min(axis=1), rh_is.max(axis=1), color=_SERIES[0], alpha=0.25, linewidth=0, label='最低〜最高')
    ax.plot(x_is, rh_is.mean(axis=1), color=_SERIES[0], linewidth=2, marker='o', markersize=4, label='平均')
    edges = np.cumsum([0.0] + [n.width for n in geometry.nodes]) * 1000.0
    for i in range(1, len(geometry.nodes)):
        if geometry.nodes[i].layer_index != geometry.nodes[i - 1].layer_index:
            ax.axvline(edges[i], color=_TEXT_2, linewidth=0.6, linestyle=':')
    ax.set_title(f'質点ごとの相対湿度（2 年目）　Δt = {tag} 秒', loc='left', color=_TEXT)
    ax.set_xlabel('外表面からの距離, mm（点線は層の境目）')
    ax.set_ylabel('%')
    ax.grid(True, color=_GRID, linewidth=0.8)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    ax.legend(loc='lower right', frameon=False, labelcolor=_TEXT)
    fig.tight_layout()
    fig.savefig(os.path.join(_OUTPUT, f'rh_profile_dt{tag}.png'), dpi=150)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--dt', type=float, default=3600.0, help='時間刻み, s')
    parser.add_argument('--years', type=int, default=2, help='計算する年数')
    args = parser.parse_args()

    os.makedirs(_OUTPUT, exist_ok=True)
    tag = f'{args.dt:g}'

    start = time.time()
    geometry, result, weather, indoor = run(args.dt, args.years)
    print(f'計算の時間 {time.time() - start:.0f} 秒、ステップ数 {len(result.iterations_ns)}、'
          f'修正の回数 平均 {result.iterations_ns.mean():.2f}・最大 {result.iterations_ns.max()}')

    np.savez(os.path.join(_OUTPUT, f'result_dt{tag}.npz'), time_ns=result.time_ns, t_is_ns=result.t_is_ns,
             mu_is_ns=result.mu_is_ns, iterations_ns=result.iterations_ns)
    plot(geometry, result, weather, indoor, tag)
    print(f'図と数値を {_OUTPUT} に書き出しました。')


if __name__ == '__main__':
    main()
