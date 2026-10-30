"""太陽位置：各時刻の太陽高度と太陽方位角（計算方法の文書の評価方法「太陽位置」）。

各値のコメントは、1 行目に「名前 記号, 単位」（無次元は -）、2 行目に「式 (k)」と補足を書く。
コメントの「式 (k)」は、計算方法の文書の太陽位置の節（docs/source/contents/03_06_eval_solar_position.md）の式の番号。
変数の名前は、その節の記号に合わせる（例：時角 τ は tau）。どれも 1 つの時刻の値なので、時刻の添え字 n は付けない。
"""

import math

from modules.time_axis import SECONDS_PER_DAY, SECONDS_PER_HOUR, SECONDS_PER_YEAR


def get_solar_position(t: float, phi_lat: float, lambda_lon: float) -> tuple[float, float, float, float]:
    """太陽高度と太陽方位角の正弦と余弦を求める。

    太陽が天頂にある場合（cos h_sun = 0）は、太陽方位角は定義されない。日本の緯度では、太陽が天頂に来ることはない。

    Args:
        t: 時刻（1 月 1 日 0 時からの経過時間）, s
        phi_lat: 緯度（北緯を正とする）, deg
        lambda_lon: 経度（東経を正とする）, deg

    Returns:
        太陽高度の正弦 sin h_sun, -、余弦 cos h_sun, -、太陽方位角（南を 0、西回りを正とする）の正弦 sin a_sun, -、余弦 cos a_sun, -
    """

    # 通算日 d（1 月 1 日を 1 とする）, d
    # 式 (1)。1 年は 365 日とし、2 年目以降も同じ暦を繰り返す
    d = int((t % SECONDS_PER_YEAR) // SECONDS_PER_DAY) + 1

    # その日の標準時 t_m（日本の標準時）, h
    # 式 (2)
    t_m = (t % SECONDS_PER_DAY) / SECONDS_PER_HOUR

    # 日角 ω, rad
    # 式 (3)
    omega = 2.0 * math.pi * d / 366.0

    # 赤緯 δ, rad
    # 式 (4)
    delta = (0.006322
             - 0.405748 * math.cos(omega + 0.153231)
             - 0.005880 * math.cos(2.0 * omega + 0.207099)
             - 0.003233 * math.cos(3.0 * omega + 0.620129))

    # 均時差 e, h
    # 式 (5)
    e = (-0.000279
         + 0.122772 * math.cos(omega + 1.498311)
         - 0.165458 * math.cos(2.0 * omega - 1.261546)
         - 0.005354 * math.cos(3.0 * omega - 1.157100))

    # 時角 τ, rad
    # 式 (6)。経度 135 deg（日本の標準時の子午線）との差の分だけ、南中の時刻がずれる
    tau = math.pi / 12.0 * ((t_m + e - 12.0) + (lambda_lon - 135.0) / 15.0)

    # 緯度 φ_lat, rad
    # 三角関数に入れるときだけ rad に直す
    phi = math.radians(phi_lat)

    # 太陽高度の正弦 sin h_sun と余弦 cos h_sun, -
    # 式 (7)・(8)。太陽高度は -π/2〜π/2 とするので、余弦は 0 以上
    sin_h_sun = math.sin(phi) * math.sin(delta) + math.cos(phi) * math.cos(delta) * math.cos(tau)
    cos_h_sun = math.sqrt(1.0 - sin_h_sun ** 2)

    # 太陽方位角の正弦 sin a_sun と余弦 cos a_sun, -
    # 式 (9)・(10)
    sin_a_sun = math.cos(delta) * math.sin(tau) / cos_h_sun
    cos_a_sun = (math.cos(delta) * math.sin(phi) * math.cos(tau) - math.sin(delta) * math.cos(phi)) / cos_h_sun

    return sin_h_sun, cos_h_sun, sin_a_sun, cos_a_sun
