"""計算の結果の書き出し：CSV（表計算のソフトで見る用）と HDF5（値を 2 進数のまま残す用）。

CSV は BOM 付きの UTF-8 で書く（Excel で開いても日本語が化けないように）。小数は、読み直したときに同じ値に戻る、いちばん短い書き方（Python の repr）で書く。
HDF5 には、CSV と同じ値に加えて、水分化学ポテンシャルの欄の状態量、修正の回数、計算の入力、止まったかどうかを入れ、
それぞれに単位と意味を属性として付ける。
"""

import csv
import os
from dataclasses import dataclass

import h5py
import numpy as np

from modules.air_water_properties import T_0
from modules.state import State
from modules.wall_geometry import AirNode, FilmNode, MaterialNode, WallGeometry
from modules.simulation import SimulationResult

# 結果のファイルの名前
NODES_CSV = 'nodes.csv'
TEMPERATURE_CSV = 'temperature.csv'
RH_CSV = 'relative_humidity.csv'
MASS_LOSS_CSV = 'mass_loss.csv'
RESULT_H5 = 'result.h5'

# CSV の文字コード（BOM 付きの UTF-8。Excel で開いても日本語が化けないように）
_CSV_ENCODING = 'utf-8-sig'


def _f(v: float) -> str:
    """CSV の小数の書き方（読み直したときに同じ値に戻る、いちばん短い書き方）。"""
    return repr(float(v))


@dataclass(frozen=True, kw_only=True)
class RunStatus:
    """計算が最後まで進んだか、止まったか。"""

    # 最後まで進んだら True
    completed: bool

    # 止まった理由（最後まで進んだら空）
    message: str = ''

    # 解けなかったステップの終わりの時刻, s（最後まで進んだら None）
    time_failed: float | None = None


def node_kind(node) -> str:
    """質点の種類：'material'（材料）、'air'（通気層の空気）、'film'（境目の水膜）。"""
    if isinstance(node, MaterialNode):
        return 'material'
    if isinstance(node, AirNode):
        return 'air'
    if isinstance(node, FilmNode):
        return 'film'
    raise TypeError(f'知らない種類の質点です（{type(node).__name__}）。')


def node_positions(geometry: WallGeometry) -> np.ndarray:
    """各質点の位置（外表面からの距離）, m, [i]。材料の質点はコントロールボリュームの左の端からの距離、空気の質点は層の真ん中、水膜の質点は境目。"""
    x_is, x_edge = [], 0.0
    for node in geometry.nodes:
        x_is.append(x_edge + (node.d_left if isinstance(node, MaterialNode) else node.width / 2))
        x_edge += node.width
    return np.array(x_is)


def rh_is_ns(t_is_ns: np.ndarray, mu_is_ns: np.ndarray) -> np.ndarray:
    """温度と、水分化学ポテンシャルの欄の値から、相対湿度, %, [i, n]。

    材料・水膜の質点の欄は状態量 s で、水分化学ポテンシャルは μ = min(s, 0)。相対湿度は 100% を超えない。
    """
    return np.vectorize(lambda t, s: State(t=t, mu=min(s, 0.0)).rh)(t_is_ns, mu_is_ns)


def _node_rows(geometry: WallGeometry) -> list[dict]:
    """質点の表の行。"""
    x_is = node_positions(geometry)
    return [{
        'index': i,
        'layer': node.layer_index,
        'kind': node_kind(node),
        'material': node.material.name if isinstance(node, MaterialNode) else '',
        'x': x_is[i],
        'width': node.width,
        'wood_decay': isinstance(node, MaterialNode) and node.wood_decay,
    } for i, node in enumerate(geometry.nodes)]


def _write_table(path: str, time_ns: np.ndarray, v_is_ns: np.ndarray):
    """行が記録した時刻、列が質点の表を書く（1 列目は経過時間, s）。"""
    header = ['time_s'] + [f'node_{i}' for i in range(v_is_ns.shape[0])]
    with open(path, 'w', encoding=_CSV_ENCODING, newline='') as f:
        w = csv.writer(f)
        w.writerow(header)
        for n, t in enumerate(time_ns):
            w.writerow([_f(t)] + [_f(v) for v in v_is_ns[:, n]])


def write_csv(out_dir: str, geometry: WallGeometry, result: SimulationResult):
    """CSV を書く：質点の表、温度（deg.C）、相対湿度（%）、木材腐朽の質量減少率（求める層があるときだけ）。

    Args:
        out_dir: 書き出すフォルダ（なければ作る）
        geometry: 壁の形
        result: 計算の結果
    """

    os.makedirs(out_dir, exist_ok=True)

    # 質点の表
    rows = _node_rows(geometry)
    with open(os.path.join(out_dir, NODES_CSV), 'w', encoding=_CSV_ENCODING, newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        for row in rows:
            w.writerow({**row, 'x': _f(row['x']), 'width': _f(row['width'])})

    # 温度, deg.C、相対湿度, %
    _write_table(os.path.join(out_dir, TEMPERATURE_CSV), result.time_ns, result.t_is_ns - T_0)
    _write_table(os.path.join(out_dir, RH_CSV), result.time_ns, rh_is_ns(result.t_is_ns, result.mu_is_ns))

    # 木材腐朽の質量減少率, -（求める層があるときだけ）
    if result.decay_last is not None:
        _write_table(os.path.join(out_dir, MASS_LOSS_CSV), result.time_ns, result.mass_loss_is_ns)


def _dataset(group, name: str, data, unit: str, description: str):
    """データを書き、単位と意味を属性に付ける。"""
    d = group.create_dataset(name, data=data)
    d.attrs['unit'] = unit
    d.attrs['description'] = description


def write_hdf5(path: str, geometry: WallGeometry, result: SimulationResult, status: RunStatus, input_text: str = ''):
    """HDF5 を書く。

    Args:
        path: 書き出すファイルの場所
        geometry: 壁の形
        result: 計算の結果
        status: 計算が最後まで進んだか、止まったか
        input_text: 計算の入力のファイルの中身（そのまま残す）
    """

    with h5py.File(path, 'w') as f:

        # 計算が最後まで進んだか、止まったか、計算の入力
        f.attrs['completed'] = status.completed
        f.attrs['message'] = status.message
        f.attrs['time_failed'] = np.nan if status.time_failed is None else status.time_failed
        f.attrs['input'] = input_text

        # 質点の表
        rows = _node_rows(geometry)
        g = f.create_group('nodes')
        _dataset(g, 'layer', np.array([r['layer'] for r in rows]), '-', '層の番号（0 から。水膜の質点は左の層）')
        _dataset(g, 'kind', np.array([r['kind'] for r in rows], dtype=h5py.string_dtype()), '-',
                 '質点の種類（material：材料、air：通気層の空気、film：境目の水膜）')
        _dataset(g, 'material', np.array([r['material'] for r in rows], dtype=h5py.string_dtype()), '-', '材料の名前（材料の質点だけ）')
        _dataset(g, 'x', np.array([r['x'] for r in rows]), 'm', '外表面からの位置')
        _dataset(g, 'width', np.array([r['width'] for r in rows]), 'm', 'コントロールボリュームの幅')
        _dataset(g, 'wood_decay', np.array([r['wood_decay'] for r in rows]), '-', '木材腐朽を求める質点なら True')

        # 記録した時刻ごとの値。質点ごとの配列の形は [i, n]（i：質点、n：記録の番号）
        _dataset(f, 'time', result.time_ns, 's', '記録した時刻（計算を始めてからの経過時間）')
        _dataset(f, 'temperature', result.t_is_ns, 'K', '温度 [i, n]')
        _dataset(f, 'mu', result.mu_is_ns, 'J/kg',
                 '水分化学ポテンシャルの欄の値 [i, n]（材料・水膜の質点は状態量 s で、水分化学ポテンシャルは min(s, 0)）')
        _dataset(f, 'relative_humidity', rh_is_ns(result.t_is_ns, result.mu_is_ns), '%', '相対湿度 [i, n]')
        _dataset(f, 'mass_loss', result.mass_loss_is_ns, '-', '木材腐朽の質量減少率 [i, n]（求めない質点は 0）')
        _dataset(f, 'iterations', result.iterations_ns, '-', '各ステップのニュートン法の修正の回数（解けたステップだけ）')
