"""計算の入力（確認済みの値だけを持つ部品）：壁、外気、室内、計算の設定。

入力のファイル（JSON）または dict から作り、値がそろっていておかしくないことを、作るときに 1 回だけ確かめる。
作ったあとは値を変えない（frozen）。壁の入力は wall_input.WallInput が確かめる。
外気の値の範囲（緯度・経度・風速の高さなど）と気象データの中身は、気象データを読むとき（OutdoorInput.load）に
outdoor_climate.WeatherOutdoor が確かめる。
"""

import json
import os
from dataclasses import dataclass
from typing import Self

from modules.indoor_climate import CosineIndoor
from modules.input_check import check_keys, check_number
from modules.outdoor_climate import RadiationInterval, RadiationKind, RadiationUnit, WeatherOutdoor
from modules.time_axis import SECONDS_PER_HOUR, SECONDS_PER_YEAR
from modules.wall_input import WallInput

# 計算の入力の dict のキー（すべて必ず指定する）
_SIMULATION_INPUT_KEYS = ('wall', 'outdoor', 'indoor', 'simulation')

# 外気の入力の dict のキー（すべて必ず指定する）
_OUTDOOR_KEYS = ('weather_file', 'latitude', 'longitude', 'radiation_kind', 'radiation_interval', 'radiation_unit',
                 'rain_interval', 'wind_height', 'wind_exponent')

# 室内の入力の種類と、種類ごとの dict のキー（'type' のほか）。省いてよいキーは 2 つ目
_INDOOR_KEYS = {
    'cosine': (('theta_in_ave', 'theta_in_amp', 'rh_in'), ('d_peak',)),
}

# 計算の設定の dict のキー（必ず指定するキーと、省いてよいキー）
_SETTINGS_KEYS = ('dt', 'years')
_SETTINGS_OPTIONAL_KEYS = ('record_interval',)

# 日射・夜間放射の単位の、入力での書き方
_RADIATION_UNITS = {'W/m2': RadiationUnit.W_PER_M2, 'MJ/(m2 h)': RadiationUnit.MJ_PER_M2_H}


def _enum_at(enum, value, name: str, where: str):
    """入力の文字列を Enum の値にする。知らない文字列ならエラー。

    Args:
        enum: Enum のクラス
        value: 入力の値
        name: キーの名前
        where: どの入力か
    """

    names = [e.value for e in enum]
    if value not in names:
        raise ValueError(f'{where} の {name} は {names} のどれかにしてください（{value!r}）。')
    return enum(value)


@dataclass(frozen=True, kw_only=True)
class OutdoorInput:
    """外気の入力：気象データのファイルと、その読み方。"""

    # 気象データのファイルの場所（入力のファイルがあるフォルダからの場所を、ここでは絶対パスにして持つ）
    weather_file: str

    # 緯度, 度
    latitude: float

    # 経度, 度
    longitude: float

    # 日射・夜間放射の値の種類
    radiation_kind: RadiationKind

    # 日射・夜間放射の値が、どの 1 時間の値か
    radiation_interval: RadiationInterval

    # 日射・夜間放射の単位
    radiation_unit: RadiationUnit

    # 降水量の値が、どの 1 時間の値か
    rain_interval: RadiationInterval

    # 気象データの風速の高さ, m
    wind_height: float

    # 風速を高さ 10 m の値に直すべき乗則のべき指数, -
    wind_exponent: float

    @classmethod
    def from_dict(cls, d: dict, base_dir: str) -> Self:
        """外気の入力の dict から作る。

        Args:
            d: 外気の入力の dict
            base_dir: 気象データのファイルの場所の起点（入力のファイルがあるフォルダ）
        """

        where = '外気の入力（outdoor）'
        check_keys(d, _OUTDOOR_KEYS, where=where)

        # 気象データのファイルの場所。文字列でなければ、またはファイルがなければエラー
        if not isinstance(d['weather_file'], str):
            raise ValueError(f'{where} の weather_file は文字列にしてください（{d["weather_file"]!r}）。')
        path = os.path.abspath(os.path.join(base_dir, d['weather_file']))
        if not os.path.isfile(path):
            raise ValueError(f'{where} の気象データのファイルがありません（{path}）。')

        # 日射・夜間放射の単位。知らない書き方ならエラー
        if d['radiation_unit'] not in _RADIATION_UNITS:
            raise ValueError(f'{where} の radiation_unit は {list(_RADIATION_UNITS)} のどれかにしてください（{d["radiation_unit"]!r}）。')

        return cls(
            weather_file=path,
            latitude=check_number(d['latitude'], 'latitude', where),
            longitude=check_number(d['longitude'], 'longitude', where),
            radiation_kind=_enum_at(RadiationKind, d['radiation_kind'], 'radiation_kind', where),
            radiation_interval=_enum_at(RadiationInterval, d['radiation_interval'], 'radiation_interval', where),
            radiation_unit=_RADIATION_UNITS[d['radiation_unit']],
            rain_interval=_enum_at(RadiationInterval, d['rain_interval'], 'rain_interval', where),
            wind_height=check_number(d['wind_height'], 'wind_height', where),
            wind_exponent=check_number(d['wind_exponent'], 'wind_exponent', where),
        )

    def load(self) -> WeatherOutdoor:
        """気象データを読む。"""

        return WeatherOutdoor.from_csv(
            self.weather_file,
            latitude=self.latitude,
            longitude=self.longitude,
            radiation_kind=self.radiation_kind,
            radiation_interval=self.radiation_interval,
            radiation_unit=self.radiation_unit,
            rain_interval=self.rain_interval,
            wind_height=self.wind_height,
            wind_exponent=self.wind_exponent,
        )


def _indoor_from_dict(d: dict) -> CosineIndoor:
    """室内の入力の dict から、室内の条件を作る（今は余弦の 1 種類）。

    Args:
        d: 室内の入力の dict（'type' で種類を指定する）
    """

    where = '室内の入力（indoor）'

    # 種類。知らない種類ならエラー
    if d.get('type') not in _INDOOR_KEYS:
        raise ValueError(f'{where} の type は {list(_INDOOR_KEYS)} のどれかにしてください（{d.get("type")!r}）。')
    required, optional = _INDOOR_KEYS[d['type']]
    check_keys(d, ('type',) + required, where=where, optional=optional)

    # いちばん暑い日は整数（CosineIndoor が確かめる）なので、数にはしない
    values = {k: check_number(d[k], k, where) for k in required}
    try:
        return CosineIndoor(**values, **{k: d[k] for k in optional if k in d})
    except ValueError as e:
        raise ValueError(f'{where}: {e}') from None


@dataclass(frozen=True, kw_only=True)
class SimulationSettings:
    """計算の設定：時間刻み、計算する年数、記録の間隔。"""

    # 時間刻み Δt, s
    dt: float

    # 計算する年数（1 年は 365 日）
    years: int

    # 記録の間隔, s（Δt の整数倍）
    record_interval: float = SECONDS_PER_HOUR

    def __post_init__(self):
        """作るときに値を確かめる。"""

        # 時間刻みが正でなければエラー
        if not self.dt > 0.0:
            raise ValueError(f'時間刻み dt は正にしてください（{self.dt}）。')

        # 年数が 1 以上の整数でなければエラー
        if isinstance(self.years, bool) or not isinstance(self.years, int) or self.years < 1:
            raise ValueError(f'計算する年数 years は 1 以上の整数にしてください（{self.years}）。')

        # 計算する期間が時間刻みの整数倍でなければエラー（最後のステップが期間の終わりに合わない）
        if (self.years * SECONDS_PER_YEAR / self.dt) % 1.0 != 0.0:
            raise ValueError(f'計算する期間（{self.years} 年）が、時間刻み dt（{self.dt} s）の整数倍になりません。')

        # 記録の間隔が時間刻みの整数倍でなければエラー
        if not self.record_interval > 0.0 or (self.record_interval / self.dt) % 1.0 != 0.0:
            raise ValueError(f'記録の間隔 record_interval（{self.record_interval} s）は、時間刻み dt（{self.dt} s）の整数倍にしてください。')

    @property
    def n_steps(self) -> int:
        """進めるステップの数。"""
        return int(round(self.years * SECONDS_PER_YEAR / self.dt))

    @property
    def record_every(self) -> int:
        """何ステップごとに記録するか。"""
        return int(round(self.record_interval / self.dt))

    @classmethod
    def from_dict(cls, d: dict) -> Self:
        """計算の設定の dict から作る。"""

        where = '計算の設定（simulation）'
        check_keys(d, _SETTINGS_KEYS, where=where, optional=_SETTINGS_OPTIONAL_KEYS)
        values = {'dt': check_number(d['dt'], 'dt', where), 'years': d['years']}
        if 'record_interval' in d:
            values['record_interval'] = check_number(d['record_interval'], 'record_interval', where)
        try:
            return cls(**values)
        except ValueError as e:
            raise ValueError(f'{where}: {e}') from None


@dataclass(frozen=True, kw_only=True)
class SimulationInput:
    """計算の入力：壁、外気、室内、計算の設定。"""

    # 壁の入力
    wall: WallInput

    # 外気の入力
    outdoor: OutdoorInput

    # 室内の条件
    indoor: CosineIndoor

    # 計算の設定
    simulation: SimulationSettings

    @classmethod
    def from_dict(cls, d: dict, base_dir: str = '.') -> Self:
        """計算の入力の dict から作る。

        Args:
            d: 計算の入力の dict
            base_dir: 気象データのファイルの場所の起点（入力のファイルがあるフォルダ）
        """

        if not isinstance(d, dict):
            raise ValueError(f'計算の入力は dict（JSON のオブジェクト）で書いてください（{type(d).__name__}）。')
        check_keys(d, _SIMULATION_INPUT_KEYS, where='計算の入力')
        for key in _SIMULATION_INPUT_KEYS:
            if not isinstance(d[key], dict):
                raise ValueError(f'計算の入力の {key} は dict（JSON のオブジェクト）で書いてください（{type(d[key]).__name__}）。')

        return cls(
            wall=WallInput.from_dict(d['wall']),
            outdoor=OutdoorInput.from_dict(d['outdoor'], base_dir=base_dir),
            indoor=_indoor_from_dict(d['indoor']),
            simulation=SimulationSettings.from_dict(d['simulation']),
        )

    @classmethod
    def from_json(cls, path: str) -> Self:
        """入力のファイル（JSON、UTF-8）から作る。気象データのファイルの場所は、入力のファイルがあるフォルダから読む。

        Args:
            path: 入力のファイルの場所
        """

        with open(path, encoding='utf-8') as f:
            try:
                d = json.load(f)
            except json.JSONDecodeError as e:
                raise ValueError(f'入力のファイルが JSON として読めません（{path}）：{e}') from None
        return cls.from_dict(d, base_dir=os.path.dirname(os.path.abspath(path)))
