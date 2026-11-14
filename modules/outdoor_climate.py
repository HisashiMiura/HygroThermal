"""外気の条件：1 時間ごとの気象データから、ステップ n の外気の境界条件を作る（計算方法の文書の評価方法「外気の条件」）。

各値のコメントは、1 行目に「名前 記号, 単位」（無次元は -）、2 行目に「式 (k)」と補足を書く。
コメントの「式 (k)」は、計算方法の文書の「外気の条件」（docs/source/contents/03_07_eval_outdoor.md）の式の番号。
変数の名前は、その節の記号に合わせる。ステップ n の値には、記号と同じく _n を付ける（例：t_y,n は t_y_n）。
経過時間の起点は 1 月 1 日 0 時、1 年は 365 日（time_axis.py）。2 年目以降も同じデータを繰り返す。
日射・夜間放射の値がどの時刻の値かは、データによって違うので設定で受け取る（その節の 6 節）。
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
from modules.air_water_properties import T_0, get_p_vs
from modules.time_axis import HOURS_PER_YEAR, SECONDS_PER_HOUR, SECONDS_PER_YEAR


class RadiationKind(Enum):
    """日射・夜間放射の値の種類（計算方法の文書の「外気の条件」の 5 節・6 節）。"""

    # 瞬時値：「h 時」の欄の値は h 時ちょうどの値。時刻の間を直線で補う
    INSTANT = 'instant'

    # 1 時間の積算値・平均値：その 1 時間の中は一定（1 時間の積算値がもとのデータと合う）
    INTERVAL = 'interval'


class RadiationInterval(Enum):
    """「h 時」の欄の値が、どの 1 時間の値か（RadiationKind.INTERVAL のときだけ使う。計算方法の文書の「外気の条件」の 6 節）。

    区間は始まりを含まず、終わりを含む（h 時ちょうどは、h 時で終わる 1 時間に入る）。
    ステップの終わりの時刻の値で解く（後退差分）ので、ステップ (h - Δt, h] には h 時で終わる 1 時間の値を使うため。
    """

    # 前：h - 1 〜 h 時（気象庁、weather.csv）
    BEFORE = 'before'

    # 中央：h - 0.5 〜 h + 0.5 時
    CENTERED = 'centered'

    # 後：h 〜 h + 1 時
    AFTER = 'after'


class RadiationUnit(Enum):
    """気象データの日射・夜間放射の単位。値は W/m2 に直すときに掛ける数 f_u（計算方法の文書の「外気の条件」の 3 節）。"""

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

    # 緯度 φ_lat, deg（北緯が正）
    latitude: float

    # 経度 λ_lon, deg（東経が正）
    longitude: float

    # 行 h の外気の温度 θ_ex,h, deg.C, [h]
    theta_ex_hs: np.ndarray

    # 行 h の外気の水蒸気圧 p_v,ex,h, Pa, [h]（式 (1)。飽和水蒸気圧を超える行は、飽和水蒸気圧に丸めてある）
    p_v_ex_hs: np.ndarray

    # 行 h の法線面直達日射量 I_dn,h, W/m2, [h]
    i_dn_hs: np.ndarray

    # 行 h の水平面天空日射量 I_sky,h, W/m2, [h]
    i_sky_hs: np.ndarray

    # 行 h の夜間放射量（実効放射量） R_eff,h, W/m2, [h]（負の値もある）
    r_eff_hs: np.ndarray

    # 日射・夜間放射の値の種類
    radiation_kind: RadiationKind

    # 日射・夜間放射の値が、どの 1 時間の値か（radiation_kind が INTERVAL のときだけ使う）
    radiation_interval: RadiationInterval = RadiationInterval.BEFORE

    # 行 h の水平面の降水量 r_h,h, mm/h, [h]（1 時間の値。省くと雨なし）
    rain_hs: np.ndarray | None = None

    # 降水量の値が、どの 1 時間の値か
    rain_interval: RadiationInterval = RadiationInterval.BEFORE

    # 行 h の風速の北向きの成分 U_N,h, m/s, [h]（気象データの風速の高さの値。風が北へ吹いていくとき正。省くと無風）
    wind_north_hs: np.ndarray | None = None

    # 行 h の風速の東向きの成分 U_E,h, m/s, [h]（同じく。風が東へ吹いていくとき正。省くと無風）
    wind_east_hs: np.ndarray | None = None

    # 気象データの風速の高さ H_met, m
    wind_height: float = 10.0

    # 風速を高さ 10 m の値に直すべき乗則のべき指数 p_w, -
    wind_exponent: float = 0.25

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

        if not isinstance(self.rain_interval, RadiationInterval):
            raise ValueError(f'rain_interval は RadiationInterval にしてください（{self.rain_interval}）。')

        # 風速の高さが 0 より大きくなければ、べき指数が 0 以上でなければエラー
        if not self.wind_height > 0.0:
            raise ValueError(f'風速の高さ wind_height は 0 より大きくしてください（{self.wind_height}）。')
        if not self.wind_exponent >= 0.0:
            raise ValueError(f'べき指数 wind_exponent は 0 以上にしてください（{self.wind_exponent}）。')

        # 降水量・風速を省いたときは、雨なし・無風（0 の並び）
        for name in ('rain_hs', 'wind_north_hs', 'wind_east_hs'):
            if getattr(self, name) is None:
                object.__setattr__(self, name, np.zeros(HOURS_PER_YEAR))

        # frozen の dataclass の項目は、ふつうの代入ができないので、この書き方で写しに入れ直す
        for name in ('theta_ex_hs', 'p_v_ex_hs', 'i_dn_hs', 'i_sky_hs', 'r_eff_hs', 'rain_hs', 'wind_north_hs', 'wind_east_hs'):
            object.__setattr__(self, name, _readonly_copy(getattr(self, name), name))

        # 降水量が負ならエラー
        if (self.rain_hs < 0.0).any():
            raise ValueError('降水量 rain_hs は 0 以上にしてください。')

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
            radiation_unit: RadiationUnit,
            rain_interval: RadiationInterval,
            wind_height: float,
            wind_exponent: float
        ) -> Self:
        """今の気象データの形の CSV（Shift_JIS、weather.csv と同じ列）から作る。

        使う列：気温（deg.C）、絶対湿度（g/kg(DA)）、気圧（Pa）、法線面直達日射量・水平面天空日射量・夜間放射（単位は radiation_unit）、
        降水量（mm/h）、南北風・東西風（m/s。風が北・東へ吹いていく向きを正とする成分）。
        水蒸気圧は p_v = p_atm x / (622 + x)。相対湿度が 100% を超える行は、水蒸気圧を飽和水蒸気圧に丸めて、行の数を数える。

        Args:
            path: CSV のファイルの場所
            latitude: 緯度, 度
            longitude: 経度, 度
            radiation_kind: 日射・夜間放射の値の種類
            radiation_interval: 日射・夜間放射の値が、どの 1 時間の値か
            radiation_unit: 日射・夜間放射の単位
            rain_interval: 降水量の値が、どの 1 時間の値か
            wind_height: 気象データの風速の高さ, m
            wind_exponent: 風速を高さ 10 m の値に直すべき乗則のべき指数, -

        Returns:
            1 時間ごとの気象データ（行 h の値は、式 (1)〜(4) で求めたもの）
        """

        d = pd.read_csv(path, encoding='shift_jis')

        # 外気の温度 θ_ex,h, deg.C, [h]
        theta_ex_hs = d['気温'].to_numpy(dtype=float)

        # 外気の絶対湿度 x_ex,h, g/kg(DA), [h]
        x_ex_hs = d['絶対湿度'].to_numpy(dtype=float)

        # 気圧 P_ex,h, Pa, [h]
        p_ex_hs = d['気圧'].to_numpy(dtype=float)

        # 飽和水蒸気圧 p_vs(θ_ex,h + 273.15), Pa, [h]
        p_vs_hs = np.array([get_p_vs(theta_ex_h + T_0) for theta_ex_h in theta_ex_hs])

        # 相対湿度が 100% を超える行, [h]
        # 式 (1) の min で飽和水蒸気圧をとる行。行の数を n_rh_clamped に残す
        over_hs = p_ex_hs * x_ex_hs / (622.0 + x_ex_hs) > p_vs_hs

        # 外気の水蒸気圧 p_v,ex,h, Pa, [h]
        # 式 (1)
        p_v_ex_hs = np.minimum(p_ex_hs * x_ex_hs / (622.0 + x_ex_hs), p_vs_hs)

        # 日射・夜間放射を W/m2 に換算する係数 f_u, -
        # 3 節の表
        f_u = radiation_unit.value

        # 法線面直達日射量 I_dn,h, W/m2, [h]
        # 式 (2)。気象データに書かれた値 Î_dn,h に f_u を掛ける
        i_dn_hs = f_u * d['法線面直達日射量'].to_numpy(dtype=float)

        # 水平面天空日射量 I_sky,h, W/m2, [h]
        # 式 (3)。気象データに書かれた値 Î_sky,h に f_u を掛ける
        i_sky_hs = f_u * d['水平面天空日射量'].to_numpy(dtype=float)

        # 夜間放射量 R_eff,h, W/m2, [h]
        # 式 (4)。気象データに書かれた値 R̂_eff,h に f_u を掛ける
        r_eff_hs = f_u * d['夜間放射'].to_numpy(dtype=float)

        # 水平面の降水量 r_h,h, mm/h, [h]
        rain_hs = d['降水量'].to_numpy(dtype=float)

        # 風速の北向き・東向きの成分 U_N,h・U_E,h, m/s, [h]（気象データの風速の高さの値）
        wind_north_hs = d['南北風'].to_numpy(dtype=float)
        wind_east_hs = d['東西風'].to_numpy(dtype=float)

        return cls(
            latitude=latitude,
            longitude=longitude,
            theta_ex_hs=theta_ex_hs,
            p_v_ex_hs=p_v_ex_hs,
            i_dn_hs=i_dn_hs,
            i_sky_hs=i_sky_hs,
            r_eff_hs=r_eff_hs,
            radiation_kind=radiation_kind,
            radiation_interval=radiation_interval,
            rain_hs=rain_hs,
            rain_interval=rain_interval,
            wind_north_hs=wind_north_hs,
            wind_east_hs=wind_east_hs,
            wind_height=wind_height,
            wind_exponent=wind_exponent,
            n_rh_clamped=int(over_hs.sum()),
        )

    def _get_instant(self, x_hs: np.ndarray, t_y_n: float) -> float:
        """瞬時値の並び x_hs の、ステップ n の値（前後の行の間を直線で補う。最後の行の次は最初の行）。

        Args:
            x_hs: 行 h の値 X_h, [h]
            t_y_n: ステップ n の、1 年の中での時刻 t_y,n（式 (5)）, h

        Returns:
            ステップ n の値 X_n

        Notes:
            計算方法の文書の「外気の条件」の 5 節「瞬時値の補い方」の式 (6)〜(9) に対応する。
        """

        # 前の行の番号 h_0,n, -
        # 式 (6)。丸めの誤差で t_y,n が 8760 になっても行の外に出ないように、8760 で割った余りをとる
        h_0_n = int(math.floor(t_y_n)) % HOURS_PER_YEAR

        # 次の行の番号 h_1,n, -
        # 式 (7)
        h_1_n = (h_0_n + 1) % HOURS_PER_YEAR

        # 前の行の時刻から測った割合 f_n, -
        # 式 (8)
        f_n = t_y_n - math.floor(t_y_n)

        # ステップ n の値 X_n
        # 式 (9)
        x_n = float(x_hs[h_0_n] + f_n * (x_hs[h_1_n] - x_hs[h_0_n]))

        return x_n

    def _get_interval(self, x_hs: np.ndarray, t_y_n: float, interval: RadiationInterval) -> float:
        """1 時間の値の並び x_hs の、ステップ n の時刻を含む 1 時間の値（区間は始まりを含まず、終わりを含む）。

        Args:
            x_hs: 行 h の値 X_h, [h]
            t_y_n: ステップ n の、1 年の中での時刻 t_y,n（式 (5)）, h
            interval: 行 h の値が、どの 1 時間の値か

        Returns:
            ステップ n の値 X_n

        Notes:
            計算方法の文書の「外気の条件」の 6 節「1 時間の値の取り出し方」の式 (10)・(11) に対応する。
        """

        # 切り上げに使う、1 年の中での時刻, h
        # 丸めの誤差で、ちょうどの時刻が次の 1 時間に入らないように、ごく小さく引いておく（文書の式にはない、計算の都合）
        t_y_n_ceil = t_y_n - 1.0e-9

        # ステップ n の時刻を含む 1 時間の値が書かれた行の番号 h_n, -
        # 式 (10)
        match interval:

            # 前：h - 1 〜 h 時。時刻を含む 1 時間の終わりの時刻の行
            case RadiationInterval.BEFORE:
                h_n = math.ceil(t_y_n_ceil) % HOURS_PER_YEAR

            # 中央：h - 0.5 〜 h + 0.5 時。30 分ずらして同じように数える
            case RadiationInterval.CENTERED:
                h_n = math.ceil(t_y_n_ceil - 0.5) % HOURS_PER_YEAR

            # 後：h 〜 h + 1 時。時刻を含む 1 時間の始まりの時刻の行
            case RadiationInterval.AFTER:
                h_n = (math.ceil(t_y_n_ceil) - 1) % HOURS_PER_YEAR

            # 上のどれでもなければエラー（RadiationInterval に値を足して、ここを直し忘れたとき）
            case _:
                raise ValueError(f'interval の {interval} の扱いが決まっていません。')

        # ステップ n の値 X_n
        # 式 (11)
        x_n = float(x_hs[h_n])

        return x_n

    def _get_radiation(self, x_hs: np.ndarray, t_y_n: float) -> float:
        """日射・夜間放射の並び x_hs の、ステップ n の値（値の種類に合わせて、式 (6)〜(9) か式 (10)・(11) で求める）。

        Args:
            x_hs: 行 h の日射・夜間放射, W/m2, [h]
            t_y_n: ステップ n の、1 年の中での時刻 t_y,n（式 (5)）, h

        Returns:
            ステップ n の日射・夜間放射, W/m2
        """

        match self.radiation_kind:

            # 瞬時値：前後の行の間を直線で補う
            case RadiationKind.INSTANT:
                return self._get_instant(x_hs, t_y_n)

            # 1 時間の値：その 1 時間の中は一定
            case RadiationKind.INTERVAL:
                return self._get_interval(x_hs, t_y_n, self.radiation_interval)

            # 上のどれでもなければエラー（RadiationKind に値を足して、ここを直し忘れたとき）
            case _:
                raise ValueError(f'radiation_kind の {self.radiation_kind} の扱いが決まっていません。')

    def _get_solar(self, t_n: float, t_y_n: float) -> SolarRadiation:
        """ステップ n の太陽の位置と、日射・夜間放射。

        Args:
            t_n: ステップ n の時刻 t_n（1 月 1 日 0 時からの経過時間。太陽位置に使う）, s
            t_y_n: ステップ n の、1 年の中での時刻 t_y,n（式 (5)。日射・夜間放射の行を探すのに使う）, h

        Returns:
            ステップ n の太陽の位置（「太陽位置」の式で求める）と、日射・夜間放射

        Notes:
            計算方法の文書の「外気の条件」の 8 節「太陽の位置」と、7 節の日射・夜間放射に対応する。
        """

        # 太陽高度・太陽方位角の正弦と余弦 sin h_sun,n・cos h_sun,n・sin a_sun,n・cos a_sun,n, -
        # 8 節。「太陽位置」により、t_y,n ではなく t_n から求める
        sin_h_sun_n, cos_h_sun_n, sin_a_sun_n, cos_a_sun_n = solar_position.get_solar_position(
            t_n=t_n, phi_lat=self.latitude, lambda_lon=self.longitude)

        return SolarRadiation(
            sin_h=sin_h_sun_n,
            cos_h=cos_h_sun_n,
            sin_a=sin_a_sun_n,
            cos_a=cos_a_sun_n,
            # 法線面直達日射量 I_dn,n, W/m2：式 (6)〜(9) か式 (10)・(11)
            i_dn=self._get_radiation(self.i_dn_hs, t_y_n),
            # 水平面天空日射量 I_sky,n, W/m2：式 (6)〜(9) か式 (10)・(11)
            i_sky=self._get_radiation(self.i_sky_hs, t_y_n),
            # 夜間放射量 R_eff,n, W/m2：式 (6)〜(9) か式 (10)・(11)
            r_eff=self._get_radiation(self.r_eff_hs, t_y_n),
        )

    def get_outdoor_boundary(self, t_n: float) -> OutdoorBoundary:
        """ステップ n の外気の境界条件（太陽の位置と日射・夜間放射を含む）。

        Args:
            t_n: ステップ n の時刻 t_n（1 月 1 日 0 時からの経過時間）, s

        Returns:
            ステップ n の外気の境界条件
        """

        # 1 年の中での時刻 t_y,n, h
        # 式 (5)。2 年目以降も同じデータを繰り返す
        t_y_n = (t_n % SECONDS_PER_YEAR) / SECONDS_PER_HOUR

        # 外気の温度 θ_ex,n, deg.C
        # 式 (6)〜(9)（瞬時値）
        theta_ex_n = self._get_instant(self.theta_ex_hs, t_y_n)

        # 外気の水蒸気圧 p_v,ex,n, Pa
        # 式 (6)〜(9)（瞬時値）
        p_v_ex_n = self._get_instant(self.p_v_ex_hs, t_y_n)

        # 外気の相対湿度 φ_ex,n, %
        # 式 (12)。行の間を補った時刻でも 100% を超えないようにする
        rh_ex_n = min(100.0 * p_v_ex_n / get_p_vs(theta_ex_n + T_0), 100.0)

        # 水平面の降水量 r_h,n, mm/h
        # 式 (10)・(11)（1 時間の値）
        r_h_n = self._get_interval(self.rain_hs, t_y_n, self.rain_interval)

        # 高さ 10 m の風速に直す係数 f_w, -
        # 式 (13)
        f_w = (10.0 / self.wind_height) ** self.wind_exponent

        # 高さ 10 m の風速の北向き・東向きの成分 U_10,N,n・U_10,E,n, m/s
        # 式 (6)〜(9)（瞬時値）と式 (14)
        u_10_n_n = f_w * self._get_instant(self.wind_north_hs, t_y_n)
        u_10_e_n = f_w * self._get_instant(self.wind_east_hs, t_y_n)

        return OutdoorBoundary.from_theta_rh(theta=theta_ex_n, rh=rh_ex_n, solar=self._get_solar(t_n, t_y_n),
                                             rain=r_h_n, wind_north=u_10_n_n, wind_east=u_10_e_n)
