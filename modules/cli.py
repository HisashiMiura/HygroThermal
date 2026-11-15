"""コマンド hygrothermal：入力のファイル（JSON）を読んで計算し、結果を書き出す。

使い方：
    hygrothermal input.json                 # 結果は input.json と同じフォルダの input_output/ に書く
    hygrothermal input.json -o result/      # 結果を result/ に書く
    python -m modules.cli input.json        # パッケージを入れずに使うとき（リポジトリの一番上のフォルダで）

計算が止まったとき（結露など）は、止まるまでの結果を書き出し、理由を表示して、終了コード 1 で終わる。
入力がおかしいときは、何も書き出さずに、理由を表示して、終了コード 2 で終わる。
"""

import argparse
import os
import sys
import time

from modules.boundary import Boundaries
from modules.newton import NewtonNotConverged
from modules.output import RESULT_H5, RunStatus, write_csv, write_hdf5
from modules.simulation import SimulationResult, simulate
from modules.simulation_input import SimulationInput
from modules.solver import PoresFilled
from modules.wall_builder import build_wall
from modules.wall_geometry import WallGeometry


def run(input_path: str) -> tuple[WallGeometry, SimulationResult, RunStatus]:
    """入力のファイルを読んで計算する。

    Args:
        input_path: 入力のファイルの場所

    Returns:
        壁の形、計算の結果（止まったときは止まるまで）、計算が最後まで進んだか

    Raises:
        ValueError: 入力がおかしいとき
    """

    sim_input = SimulationInput.from_json(input_path)
    geometry, state_0 = build_wall(sim_input.wall)
    weather = sim_input.outdoor.load()
    indoor = sim_input.indoor

    def boundaries_at(time_s: float) -> Boundaries:
        return Boundaries(outdoor=weather.get_outdoor_boundary(time_s), indoor=indoor.get_indoor_boundary(time_s))

    settings = sim_input.simulation
    try:
        result = simulate(geometry, state_0, dt=settings.dt, n_steps=settings.n_steps, boundaries_at=boundaries_at,
                          record_every=settings.record_every)
        return geometry, result, RunStatus(completed=True)
    except (NewtonNotConverged, PoresFilled) as e:
        status = RunStatus(completed=False, message=f'{type(e).__name__}: {e}', time_failed=e.time_failed)
        return geometry, e.partial_result, status


def write(out_dir: str, geometry: WallGeometry, result: SimulationResult, status: RunStatus, input_path: str):
    """結果を書き出す（CSV と HDF5）。

    Args:
        out_dir: 書き出すフォルダ
        geometry: 壁の形
        result: 計算の結果
        status: 計算が最後まで進んだか
        input_path: 入力のファイルの場所（中身を HDF5 に残す）
    """

    write_csv(out_dir, geometry, result)
    with open(input_path, encoding='utf-8') as f:
        input_text = f.read()
    write_hdf5(os.path.join(out_dir, RESULT_H5), geometry, result, status, input_text)


def _default_out_dir(input_path: str) -> str:
    """結果を書くフォルダの既定：入力のファイルと同じフォルダの「入力のファイルの名前_output」。"""
    stem = os.path.splitext(os.path.basename(input_path))[0]
    return os.path.join(os.path.dirname(os.path.abspath(input_path)), f'{stem}_output')


def main(argv: list[str] | None = None) -> int:
    """コマンドの入口。

    Args:
        argv: コマンドの引数（None ならコマンドラインから）

    Returns:
        終了コード（0：最後まで進んだ、1：計算が止まった、2：入力がおかしい）
    """

    parser = argparse.ArgumentParser(prog='hygrothermal', description='壁の熱・水分同時移動の計算')
    parser.add_argument('input', help='入力のファイル（JSON）')
    parser.add_argument('-o', '--output', help='結果を書き出すフォルダ（省くと、入力のファイルと同じフォルダの「名前_output」）')
    args = parser.parse_args(argv)
    out_dir = args.output or _default_out_dir(args.input)

    start = time.time()
    try:
        geometry, result, status = run(args.input)
    except (ValueError, OSError) as e:
        print(f'入力がおかしいので、計算できません：{e}', file=sys.stderr)
        return 2

    write(out_dir, geometry, result, status, args.input)
    seconds = time.time() - start
    if status.completed:
        print(f'計算が終わりました（{seconds:.0f} 秒、{len(result.iterations_ns)} ステップ）。結果：{out_dir}')
        return 0
    print(f'経過時間 {status.time_failed:.0f} s（{status.time_failed / 86400:.2f} 日）のステップで計算が止まりました：{status.message}',
          file=sys.stderr)
    print(f'止まるまでの結果を書き出しました：{out_dir}', file=sys.stderr)
    return 1


if __name__ == '__main__':
    sys.exit(main())
