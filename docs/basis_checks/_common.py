"""根拠のページの試算に共通の部品。

例の壁（examples/one_year_main_wall.py の 7 層の壁）の 2 年間の計算を、計算の部品を一時的に差し替えて回し、
結果を docs/basis_checks/_out/ に保存して使い回す。保存した結果があれば、計算し直さない（--force で計算し直す）。
"""

import contextlib
import os
import sys
import time

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(_HERE, '..', '..'))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, 'examples'))

import one_year_main_wall as ex  # noqa: E402
from modules.state import State  # noqa: E402
from modules.time_axis import SECONDS_PER_YEAR  # noqa: E402
from modules.wall_builder import build_wall  # noqa: E402

# 結果を保存するフォルダ（git に入れない）
OUT = os.path.join(_HERE, '_out')

# 計算し直すか（コマンドラインに --force があれば計算し直す）
FORCE = '--force' in sys.argv


@contextlib.contextmanager
def patched(*items):
    """(対象, 属性の名前, 値) の組を、with の中だけ差し替える。"""
    saved = [(obj, name, getattr(obj, name)) for obj, name, _ in items]
    try:
        for obj, name, value in items:
            setattr(obj, name, value)
        yield
    finally:
        for obj, name, value in saved:
            setattr(obj, name, value)


def run_example(name: str, dt: float = 3600.0, years: int = 2, patches: tuple = ()) -> dict:
    """例の壁を years 年間計算した結果（保存した結果があれば、それを読む）。

    Args:
        name: 結果の名前（_out/{name}.npz に保存する）
        dt: 時間刻み, s
        years: 計算する年数
        patches: patched に渡す (対象, 属性の名前, 値) の組の並び

    Returns:
        time_ns・t_is_ns・mu_is_ns・iterations_ns・mass_loss_is_ns の dict（記録は 1 時間ごと）
    """
    path = os.path.join(OUT, f'{name}.npz')
    if os.path.exists(path) and not FORCE:
        return dict(np.load(path))
    os.makedirs(OUT, exist_ok=True)
    start = time.time()
    with patched(*patches):
        _, r, _, _ = ex.run(dt, years)
    d = {'time_ns': r.time_ns, 't_is_ns': r.t_is_ns, 'mu_is_ns': r.mu_is_ns, 'iterations_ns': r.iterations_ns,
         'mass_loss_is_ns': r.mass_loss_is_ns, 'seconds': np.array(time.time() - start)}
    np.savez(path, **d)
    print(f'  （{name} を計算しました：{d["seconds"]:.0f} 秒）')
    return d


def example_geometry():
    """例の壁の形。"""
    geometry, _ = build_wall(ex._main_wall_input())
    return geometry


def year2(d: dict) -> np.ndarray:
    """2 年目の記録か, [n]。"""
    return (d['time_ns'] > SECONDS_PER_YEAR) & (d['time_ns'] <= 2 * SECONDS_PER_YEAR)


def rh_is_ns(t_is_ns: np.ndarray, mu_is_ns: np.ndarray) -> np.ndarray:
    """相対湿度, %（状態量が正の質点は、水分化学ポテンシャル 0 として 100%）。"""
    return np.vectorize(lambda t, mu: State(t=t, mu=min(mu, 0.0)).rh)(t_is_ns, mu_is_ns)


def daily_max(x_n: np.ndarray) -> np.ndarray:
    """1 時間ごとの 1 年分（8760 個）から、1 日の最大, [日]。"""
    return x_n[:8760].reshape(365, 24).max(axis=1)
