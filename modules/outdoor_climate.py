"""外気の気候：1 時間ごとの気象データから、経過時間の外気の境界条件を作る（段階 6）。

経過時間の起点は 1 月 1 日 0 時、1 年は 365 日（time_axis.py）。2 年目以降も同じデータを繰り返す。
日射・夜間放射の値がどの時刻の値かは、データによって違うので設定で受け取る（docs/weather_data_time_convention.md）。
今の weather.Weather・OutdoorCondition は、置き換えが終わったら消す。風・雨は段階 7 で足す。
"""

import math
from dataclasses import dataclass
from enum import Enum
from typing import Self

import numpy as np
import pandas as pd

from modules import solar_position
from modules.boundary import OutdoorBoundary, SolarRadiation
from modules.thermo_dynamics import ATP, GOFF
from modules.time_axis import HOURS_PER_YEAR, SECONDS_PER_HOUR, SECONDS_PER_YEAR


class RadiationKind(Enum):
    """日射・夜間放射の値の種類（docs/weather_data_time_convention.md 2 章）。"""

    # 瞬時値（⑥）：「h 時」の欄の値は h 時ちょうどの値。時刻の間を直線で補う
    INSTANT = 'instant'

    # 1 時間の積算値・平均値（①〜⑤）：その 1 時間の中は一定（1 時間の積算値がもとのデータと合う）
    INTERVAL = 'interval'


class RadiationInterval(Enum):
    """「h 時」の欄の値が、どの 1 時間の値か（RadiationKind.INTERVAL のときだけ使う）。

    区間は始まりを含まず、終わりを含む（h 時ちょうどは、h 時で終わる 1 時間に入る）。
    ステップの終わりの時刻の値で解く（後退差分）ので、ステップ (h - Δt, h] には h 時で終わる 1 時間の値を使うため。
    """

    # h - 1 〜 h 時（①・④。気象庁、weather.csv）
    BEFORE = 'before'

    # h - 0.5 〜 h + 0.5 時（③）
    CENTERED = 'centered'

    # h 〜 h + 1 時（②・⑤）
    AFTER = 'after'


class RadiationUnit(Enum):
    """気象データの日射・夜間放射の単位。値は W/m2 に直すときに掛ける数。"""

    # W/m2（平均値・瞬時値）
    W_PER_M2 = 1.0

    # MJ/(m2 h)（1 時間の積算値。weather.csv）
    MJ_PER_M2_H = 1.0e6 / 3600.0


def _readonly_copy(v_hs: np.ndarray, name: str) -> np.ndarray:
    """1 年分（8760 個）の小数の配列の写しを、読み取り専用にして返す。"""

    v_hs = np.array(v_hs, dtype=float)

    # 1 次元で 8760 個でなければエラー
    if v_hs.shape != (HOURS_PER_YEAR,):
        raise ValueError(f'{name} は 1 時間ごと 1 年分（{HOURS_PER_YEAR} 個）の 1 次元の配列にしてください（{v_hs.shape}）。')

    v_hs.flags.writeable = False
    return v_hs


@dataclass(frozen=True, eq=False)
class WeatherOutdoor:
    """1 時間ごとの気象データ（1 年分、8760 行）から、外気の境界条件を作る。

    行 h は 1 月 1 日 0 時から h 時間後（添え字 h）。気温・水蒸気圧は h 時の瞬時値として、時刻の間を直線で補う。
    12 月 31 日 23 時と次の年の 1 月 1 日 0 時の間は、最後の行と最初の行の間を補う。
    numpy の配列どうしの == は True/False の配列を返すので、自動の == は作らない（eq=False）。
    """

    # 緯度, 度（北緯が正）
    latitude: float

    # 経度, 度（東経が正）
    longitude: float

    # 気温, deg.C, [h]
    theta_hs: np.ndarray

    # 水蒸気圧, Pa, [h]（飽和水蒸気圧を超える行は、読み込むときに飽和水蒸気圧に丸めてある）
    p_v_hs: np.ndarray

    # 法線面直達日射量, W/m2, [h]
    i_dn_hs: np.ndarray

    # 水平面天空日射量, W/m2, [h]
    i_sky_hs: np.ndarray

    # 夜間放射量（実効放射量）, W/m2, [h]（負の値もある）
    r_eff_hs: np.ndarray

    # 日射・夜間放射の値の種類
    radiation_kind: RadiationKind

    # 日射・夜間放射の値が、どの 1 時間の値か（radiation_kind が INTERVAL のときだけ使う）
    radiation_interval: RadiationInterval = RadiationInterval.BEFORE

    # 読み込むときに、相対湿度が 100% を超えていて 100% に丸めた行の数
    n_rh_clamped: int = 0

    def __post_init__(self):
        """作るときに値を確かめ、配列を書き換えられないようにする。"""

        # 緯度が -90〜90 度、経度が -180〜180 度でなければエラー
        if not -90.0 <= self.latitude <= 90.0:
            raise ValueError(f'緯度 latitude は -90〜90 度にしてください（{self.latitude}）。')
        if not -180.0 <= self.longitude <= 180.0:
            raise ValueError(f'経度 longitude は -180〜180 度にしてください（{self.longitude}）。')

        # 値の種類・区間の設定が、決めた選択肢でなければエラー
        if not isinstance(self.radiation_kind, RadiationKind):
            raise ValueError(f'radiation_kind は RadiationKind にしてください（{self.radiation_kind}）。')
        if not isinstance(self.radiation_interval, RadiationInterval):
            raise ValueError(f'radiation_interval は RadiationInterval にしてください（{self.radiation_interval}）。')

        # frozen の dataclass の項目は、ふつうの代入ができないので、この書き方で写しに入れ直す
        for name in ('theta_hs', 'p_v_hs', 'i_dn_hs', 'i_sky_hs', 'r_eff_hs'):
            object.__setattr__(self, name, _readonly_copy(getattr(self, name), name))

        # 日射量が負ならエラー（夜間放射は負の値もあるので確かめない）
        if (self.i_dn_hs < 0.0).any() or (self.i_sky_hs < 0.0).any():
            raise ValueError('日射量 i_dn_hs・i_sky_hs は 0 以上にしてください。')

    @classmethod
    def from_csv(
            cls,
            path: str,
            latitude: float,
            longitude: float,
            radiation_kind: RadiationKind,
            radiation_interval: RadiationInterval,
            radiation_unit: RadiationUnit
        ) -> Self:
        """今の気象データの形の CSV（Shift_JIS、weather.csv と同じ列）から作る。

        使う列：気温（deg.C）、絶対湿度（g/kg(DA)）、気圧（Pa）、法線面直達日射量・水平面天空日射量・夜間放射（単位は radiation_unit）。
        水蒸気圧は p_v = p_atm x / (622 + x)。相対湿度が 100% を超える行は、水蒸気圧を飽和水蒸気圧に丸めて、行の数を数える。

        Args:
            path: CSV のファイルの場所
            latitude: 緯度, 度
            longitude: 経度, 度
            radiation_kind: 日射・夜間放射の値の種類
            radiation_interval: 日射・夜間放射の値が、どの 1 時間の値か
            radiation_unit: 日射・夜間放射の単位
        """

        d = pd.read_csv(path, encoding='shift_jis')

        # 気温, deg.C, [h]
        theta_hs = d['気温'].to_numpy(dtype=float)

        # 絶対湿度, g/kg(DA)、気圧, Pa から、水蒸気圧, Pa, [h]
        x_hs = d['絶対湿度'].to_numpy(dtype=float)
        p_v_hs = d['気圧'].to_numpy(dtype=float) * x_hs / (622.0 + x_hs)

        # 飽和水蒸気圧, Pa, [h]
        p_v_sat_hs = np.array([GOFF(theta_h + ATP)[1] for theta_h in theta_hs])

        # 相対湿度が 100% を超える行は、水蒸気圧を飽和水蒸気圧に丸め、その行の数を数える
        over_hs = p_v_hs > p_v_sat_hs
        p_v_hs = np.where(over_hs, p_v_sat_hs, p_v_hs)

        # 日射・夜間放射を W/m2 に直す
        f = radiation_unit.value

        return cls(
            latitude=latitude,
            longitude=longitude,
            theta_hs=theta_hs,
            p_v_hs=p_v_hs,
            i_dn_hs=d['法線面直達日射量'].to_numpy(dtype=float) * f,
            i_sky_hs=d['水平面天空日射量'].to_numpy(dtype=float) * f,
            r_eff_hs=d['夜間放射'].to_numpy(dtype=float) * f,
            radiation_kind=radiation_kind,
            radiation_interval=radiation_interval,
            n_rh_clamped=int(over_hs.sum()),
        )

    def _instant(self, v_hs: np.ndarray, time: float) -> float:
        """瞬時値の並び v_hs の、経過時間 time（s）の値（前後の行の間を直線で補う。最後の行の次は最初の行）。"""

        # 1 年の中での時刻, h（2 年目以降も同じデータを繰り返す）
        hour = (time % SECONDS_PER_YEAR) / SECONDS_PER_HOUR

        # 前の行と次の行の番号、前の行からの割合
        h_0 = int(math.floor(hour)) % HOURS_PER_YEAR
        h_1 = (h_0 + 1) % HOURS_PER_YEAR
        frac = hour - math.floor(hour)

        return float(v_hs[h_0] + frac * (v_hs[h_1] - v_hs[h_0]))

    def _interval(self, v_hs: np.ndarray, time: float) -> float:
        """1 時間の値の並び v_hs の、経過時間 time（s）を含む 1 時間の値（区間は始まりを含まず、終わりを含む）。"""

        # 1 年の中での時刻, h（丸めの誤差で、ちょうどの時刻が次の 1 時間に入らないように、ごく小さく引いてから切り上げる）
        hour = (time % SECONDS_PER_YEAR) / SECONDS_PER_HOUR - 1.0e-9

        # time を含む 1 時間の値が書かれた欄の番号
        match self.radiation_interval:

            # 「h 時」の欄が h - 1 〜 h 時：time を含む 1 時間の終わりの時刻の欄
            case RadiationInterval.BEFORE:
                h = math.ceil(hour)

            # 「h 時」の欄が h - 0.5 〜 h + 0.5 時：30 分ずらして同じように数える
            case RadiationInterval.CENTERED:
                h = math.ceil(hour - 0.5)

            # 「h 時」の欄が h 〜 h + 1 時：time を含む 1 時間の始まりの時刻の欄
            case RadiationInterval.AFTER:
                h = math.ceil(hour) - 1

        return float(v_hs[h % HOURS_PER_YEAR])

    def _radiation(self, v_hs: np.ndarray, time: float) -> float:
        """日射・夜間放射の並び v_hs の、経過時間 time（s）の値（値の種類に合わせて取り出す）。"""

        # 瞬時値なら直線で補い、1 時間の値なら、その 1 時間の中は一定
        if self.radiation_kind == RadiationKind.INSTANT:
            return self._instant(v_hs, time)
        return self._interval(v_hs, time)

    def _get_solar(self, time: float) -> SolarRadiation:
        """経過時間 time（s）の太陽の位置と、日射・夜間放射。太陽の位置は、その時刻で式から求める。"""

        sin_h, cos_h, sin_a, cos_a = solar_position.get_solar_position(
            t=time, phi_lat=self.latitude, lambda_lon=self.longitude)

        return SolarRadiation(
            sin_h=float(sin_h),
            cos_h=float(cos_h),
            sin_a=float(sin_a),
            cos_a=float(cos_a),
            i_dn=self._radiation(self.i_dn_hs, time),
            i_sky=self._radiation(self.i_sky_hs, time),
            r_eff=self._radiation(self.r_eff_hs, time),
        )

    def get_outdoor_boundary(self, time: float) -> OutdoorBoundary:
        """1 月 1 日 0 時からの経過時間 time（s）の外気の境界条件（太陽の位置と日射を含む）。"""

        # 気温, deg.C と水蒸気圧, Pa（どちらも瞬時値として、行の間を直線で補う）
        theta = self._instant(self.theta_hs, time)
        p_v = self._instant(self.p_v_hs, time)

        # 相対湿度, %（行の間を補った時刻でも 100% を超えないように丸める）
        rh = min(100.0 * p_v / GOFF(theta + ATP)[1], 100.0)

        return OutdoorBoundary.from_theta_rh(theta=theta, rh=rh, solar=self._get_solar(time))
