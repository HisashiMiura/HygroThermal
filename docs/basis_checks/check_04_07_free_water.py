"""04_07「材料の中の自由水のモデル」：表 1・表 2、5.4 節、7 節（空隙率を与えたときに止まる日）。"""

import numpy as np

import _common as c
from modules.air_water_properties import RHO_W
from modules.boundary import IndoorBoundary, OutdoorBoundary
from modules.material_loaders import load_materials
from modules.solver import PoresFilled, step
from modules.wall_builder import build_wall
from modules.wall_input import WallInput


def table_1():
    """表 1：相対湿度 100% の含水量と、乾燥密度と実質の密度から見積もった、空隙を水で満たす量。"""
    ms = {m.name: m for m in load_materials('whamcavity')}
    print('表 1　相対湿度 100% の含水量と、空隙を水で満たす量の見積もり（kg/m3）')
    for name, rho_solid in (('グラスウール1', 2500.0), ('石膏ボード', 2300.0), ('構造用合板1', 1500.0), ('ALC', 2700.0)):
        m = ms[name]
        print(f'  {name}：乾燥密度 {m.rho_s:.0f}、100% の含水量 {m.rho_s * m.get_u(100.0):.0f}、'
              f'空隙を水で満たす量 {(1.0 - m.rho_s / rho_solid) * RHO_W:.0f}')


def table_2_and_5_4():
    """表 2：100% の含水量が水の密度を超える材料。5.4 節：乾燥密度と傾きの積の最小（吸湿する材料）。"""
    print('表 2　相対湿度 100% の含水量が 998 kg/m3 を超える材料')
    slopes = []
    for src in ('whamcavity', 'mews_ashrae'):
        for m in load_materials(src):
            if m.rho_s == 0.0:
                continue
            w_100 = m.rho_s * m.get_u(100.0)
            if w_100 > RHO_W:
                print(f'  {m.name}：{w_100:.0f}')
            slope = m.rho_s * (m.get_u(100.0) - m.get_u(99.99)) / 1.0e-4
            if slope > 0.0:
                slopes.append((slope, m.name))
    print(f'5.4 節　吸湿する材料の、乾燥密度 × 傾きの最小：{min(slopes)[0]:.0f} kg/m3（{min(slopes)[1]}）')


def section_7():
    """7 節：構造用合板 1 層に、外気 -5 deg.C・90%、室内 25 deg.C・95% を 30 日（1 時間刻み）。空隙率による違い。"""
    outdoor = OutdoorBoundary.from_theta_rh(theta=-5.0, rh=90.0)
    indoor = IndoorBoundary.from_theta_rh(theta=25.0, rh=95.0)
    print('7 節　空隙率を与えた場合')
    for porosity in (0.64, 0.864):
        g, st = build_wall(WallInput.from_dict({
            'surface': {'direction': 'e', 'angle': 90.0, 'albedo': 0.3, 'absorptance': 0.9, 'emissivity': 0.9,
                        'exposure_factor': 1.0, 'deposition_factor': 1.0},
            'layers': [{'source': 'whamcavity', 'name': '構造用合板1', 'thickness': 0.009, 'n_div': 6,
                        'initial_temperature': 0.0, 'initial_rh': 95.0, 'porosity': porosity}],
            'interfaces': [{'g_h': 22.4, 'g_m': 2.0e-11}, {'g_h': 9.2, 'g_m': 2.6e-10}],
            'len_vertical': 1.0, 'len_horizontal': 1.0}))
        try:
            for n in range(24 * 30):
                st, _ = step(g, st, 3600.0, outdoor, outdoor, indoor)
            print(f'  空隙率 {porosity}：30 日で上限を超えなかった')
        except PoresFilled as e:
            print(f'  空隙率 {porosity}：{(n + 1) / 24:.1f} 日目に止まった（{e}）')


def main():
    table_1()
    table_2_and_5_4()
    section_7()


if __name__ == '__main__':
    main()
