import csv
import json
import os

import h5py
import numpy as np
import pytest

from modules import cli
from modules.air_water_properties import T_0
from modules.output import rh_is_ns

_ROOT = os.path.join(os.path.dirname(__file__), '..')


def _input(porosity: float | None = None, wood_decay: bool = True) -> dict:
    """構造用合板 9 mm 1 層を、例の気象データで 1 年（1 日刻み）計算する入力。"""
    layer = {'source': 'whamcavity', 'name': '構造用合板1', 'thickness': 0.009, 'n_div': 4,
             'initial_temperature': 10.0, 'initial_rh': 80.0, 'wood_decay': wood_decay}
    if porosity is not None:
        layer['porosity'] = porosity
    return {
        'wall': {
            'surface': {'direction': 'e', 'angle': 90.0, 'albedo': 0.3, 'absorptance': 0.9, 'emissivity': 0.9,
                        'exposure_factor': 1.0, 'deposition_factor': 1.0},
            'layers': [layer],
            'interfaces': [{'g_h': 22.4, 'g_m': 2.0e-11}, {'g_h': 9.2, 'g_m': 2.6e-10}],
            'len_vertical': 1.0, 'len_horizontal': 1.0},
        'outdoor': {'weather_file': os.path.abspath(os.path.join(_ROOT, 'modules', 'weather.csv')),
                    'latitude': 36.7, 'longitude': 137.21, 'radiation_kind': 'interval', 'radiation_interval': 'before',
                    'radiation_unit': 'MJ/(m2 h)', 'rain_interval': 'before', 'wind_height': 6.5, 'wind_exponent': 0.25},
        'indoor': {'type': 'cosine', 'theta_in_ave': 22.5, 'theta_in_amp': 4.5, 'rh_in': 60.0},
        'simulation': {'dt': 86400, 'years': 1, 'record_interval': 86400},
    }


def _write_input(tmp_path, d: dict) -> str:
    p = tmp_path / 'case.json'
    p.write_text(json.dumps(d, ensure_ascii=False), encoding='utf-8')
    return str(p)


def _read_table(path: str) -> tuple[list[str], np.ndarray]:
    """結果の表（1 行目が見出し）を読む。値の形は [行, 列]。"""
    with open(path, encoding='utf-8-sig', newline='') as f:
        rows = list(csv.reader(f))
    return rows[0], np.array([[float(v) for v in r] for r in rows[1:]])


@pytest.fixture(scope='module')
def completed(tmp_path_factory):
    """最後まで進む計算の、終了コード・結果を書いたフォルダ・計算の結果。"""
    tmp_path = tmp_path_factory.mktemp('completed')
    path = _write_input(tmp_path, _input())
    _, result, _ = cli.run(path)
    code = cli.main([path])
    return code, str(tmp_path / 'case_output'), result


class TestCompleted:

    def test_exit_code(self, completed):
        assert completed[0] == 0

    def test_files(self, completed):
        # 1 層の壁なので、水膜の質点はない（水膜の水分量と、流れ落ちた水の CSV は書かない）
        assert sorted(os.listdir(completed[1])) == ['free_water.csv', 'mass_loss.csv', 'nodes.csv', 'relative_humidity.csv', 'result.h5',
                                                    'temperature.csv', 'water_content.csv']

    def test_nodes_csv(self, completed):
        with open(os.path.join(completed[1], 'nodes.csv'), encoding='utf-8-sig', newline='') as f:
            rows = list(csv.DictReader(f))
        assert [r['kind'] for r in rows] == ['material'] * 4
        assert sum(float(r['width']) for r in rows) == pytest.approx(0.009)
        assert {r['material'] for r in rows} == {'構造用合板1'}

    def test_csv_bom(self, completed):
        # CSV は BOM 付きの UTF-8（Excel で開いても日本語が化けない）
        for name in ('nodes.csv', 'temperature.csv', 'relative_humidity.csv', 'mass_loss.csv'):
            with open(os.path.join(completed[1], name), 'rb') as f:
                assert f.read(3) == b'\xef\xbb\xbf'

    def test_csv_same_values(self, completed):
        # CSV の値は、読み直すと計算の結果と同じ値に戻る
        _, out_dir, result = completed
        header, v = _read_table(os.path.join(out_dir, 'temperature.csv'))
        assert header == ['time_s', 'node_0', 'node_1', 'node_2', 'node_3']
        np.testing.assert_array_equal(v[:, 0], result.time_ns)
        np.testing.assert_array_equal(v[:, 1:].T, result.t_is_ns - T_0)
        _, v = _read_table(os.path.join(out_dir, 'relative_humidity.csv'))
        np.testing.assert_array_equal(v[:, 1:].T, rh_is_ns(result.t_is_ns, result.mu_is_ns))
        _, v = _read_table(os.path.join(out_dir, 'mass_loss.csv'))
        np.testing.assert_array_equal(v[:, 1:].T, result.mass_loss_is_ns)
        assert len(result.time_ns) == 366

    def test_water_content(self, completed):
        # 材料の含水量は、水分量（kg/m2）をコントロールボリュームの幅で割ったもの。HDF5 の水分量と合う
        _, out_dir, result = completed
        header, v = _read_table(os.path.join(out_dir, 'water_content.csv'))
        assert header == ['time_s', 'node_0', 'node_1', 'node_2', 'node_3']
        with h5py.File(os.path.join(out_dir, 'result.h5'), 'r') as f:
            w = f['water_amount'][()]
            width = f['nodes/width'][()]
            assert f['water_amount'].attrs['unit'] == 'kg/m2'
            np.testing.assert_array_equal(f['drained'][()], np.zeros_like(w))
        np.testing.assert_allclose(v[:, 1:].T, w / width[:, None], rtol=1e-15)
        # この計算では自由水はたまらないので、含水量は、相対湿度 100% の平衡含水量（構造用合板1 は 637.5 kg/m3）を超えない
        _, free = _read_table(os.path.join(out_dir, 'free_water.csv'))
        np.testing.assert_array_equal(free[:, 1:], 0.0)
        assert 0.0 < v[:, 1:].min() and v[:, 1:].max() < 637.5

    def test_hdf5(self, completed):
        _, out_dir, result = completed
        with h5py.File(os.path.join(out_dir, 'result.h5'), 'r') as f:
            assert f.attrs['completed']
            assert json.loads(f.attrs['input'])['simulation']['years'] == 1
            np.testing.assert_array_equal(f['time'][()], result.time_ns)
            np.testing.assert_array_equal(f['temperature'][()], result.t_is_ns)
            np.testing.assert_array_equal(f['mu'][()], result.mu_is_ns)
            np.testing.assert_array_equal(f['iterations'][()], result.iterations_ns)
            np.testing.assert_array_equal(f['mass_loss'][()], result.mass_loss_is_ns)
            assert f['temperature'].attrs['unit'] == 'K'
            assert [k.decode() for k in f['nodes/kind'][()]] == ['material'] * 4

    def test_no_mass_loss_csv(self, tmp_path):
        # 木材腐朽を求める層がなければ、質量減少率の CSV は書かない
        path = _write_input(tmp_path, _input(wood_decay=False))
        assert cli.main([path, '-o', str(tmp_path / 'out')]) == 0
        assert 'mass_loss.csv' not in os.listdir(tmp_path / 'out')


class TestStopped:

    def test_partial(self, tmp_path, capsys):
        # 空隙を水で満たす量が小さいと、最初のステップで止まる。止まるまで（最初の状態）を書き出す
        path = _write_input(tmp_path, _input(porosity=0.05))
        assert cli.main([path, '-o', str(tmp_path / 'out')]) == 1
        assert 'PoresFilled' in capsys.readouterr().err
        _, v = _read_table(os.path.join(tmp_path, 'out', 'temperature.csv'))
        assert v.shape[0] == 1
        with h5py.File(os.path.join(tmp_path, 'out', 'result.h5'), 'r') as f:
            assert not f.attrs['completed']
            assert f.attrs['time_failed'] == 86400.0
            assert 'PoresFilled' in f.attrs['message']


class TestInputError:

    def test_bad_input(self, tmp_path, capsys):
        d = _input()
        d['simulation']['dt'] = -1
        path = _write_input(tmp_path, d)
        assert cli.main([path]) == 2
        assert 'dt' in capsys.readouterr().err
        assert not os.path.exists(tmp_path / 'case_output')

    def test_no_file(self, tmp_path, capsys):
        assert cli.main([str(tmp_path / 'none.json')]) == 2
