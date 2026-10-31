"""太陽位置：各時刻の太陽高度と太陽方位角（計算方法の文書の評価方法「太陽位置」）。

各値のコメントは、1 行目に「名前 記号, 単位」（無次元は -）、2 行目に「式 (k)」と補足を書く。
コメントの「式 (k)」は、計算方法の文書の太陽位置の節（docs/source/contents/03_06_eval_solar_position.md）の式の番号。
変数の名前は、その節の記号に合わせる（例：時角 τ は tau）。ステップ n の値には、記号と同じく _n を付ける（例：τ^n は tau_n）。時刻によらない値（緯度 φ_lat）には付けない。
"""

import math

from modules.time_axis import SECONDS_PER_DAY, SECONDS_PER_HOUR, SECONDS_PER_YEAR


def get_solar_position(t_n: float, phi_lat: float, lambda_lon: float) -> tuple[float, float, float, float]:
    """ステップ n の、太陽高度と太陽方位角の正弦と余弦を求める。

    太陽が天頂にある場合（cos h_sun^n = 0）は、太陽方位角は定義されない。日本の緯度では、太陽が天頂に来ることはない。

    Args:
        t_n: ステップ n の時刻 t^n（1 月 1 日 0 時からの経過時間）, s
        phi_lat: 緯度（北緯を正とする）, deg
        lambda_lon: 経度（東経を正とする）, deg

    Returns:
        ステップ n の、太陽高度の正弦 sin h_sun^n, -、余弦 cos h_sun^n, -、太陽方位角（南を 0、西回りを正とする）の正弦 sin a_sun^n, -、余弦 cos a_sun^n, -
    """

    # 通算日 d^n（1 月 1 日を 1 とする）, d
    # 式 (1)。1 年は 365 日とし、2 年目以降も同じ暦を繰り返す
    d_n = int((t_n % SECONDS_PER_YEAR) // SECONDS_PER_DAY) + 1

    # その日の標準時 t_m^n（日本の標準時）, h
    # 式 (2)
    t_m_n = (t_n % SECONDS_PER_DAY) / SECONDS_PER_HOUR

    # 日角 ω^n, rad
    # 式 (3)
    omega_n = 2.0 * math.pi * d_n / 366.0

    # 赤緯 δ^n, rad
    # 式 (4)
    delta_n = (0.006322
               - 0.405748 * math.cos(omega_n + 0.153231)
               - 0.005880 * math.cos(2.0 * omega_n + 0.207099)
               - 0.003233 * math.cos(3.0 * omega_n + 0.620129))

    # 均時差 e^n, h
    # 式 (5)
    e_n = (-0.000279
           + 0.122772 * math.cos(omega_n + 1.498311)
           - 0.165458 * math.cos(2.0 * omega_n - 1.261546)
           - 0.005354 * math.cos(3.0 * omega_n - 1.157100))

    # 時角 τ^n, rad
    # 式 (6)。経度 135 deg（日本の標準時の子午線）との差の分だけ、南中の時刻がずれる
    tau_n = math.pi / 12.0 * ((t_m_n + e_n - 12.0) + (lambda_lon - 135.0) / 15.0)

    # 緯度 φ_lat, rad
    # 三角関数に入れるときだけ rad に直す
    phi = math.radians(phi_lat)

    # 太陽高度の正弦 sin h_sun^n と余弦 cos h_sun^n, -
    # 式 (7)・(8)。太陽高度は -π/2〜π/2 とするので、余弦は 0 以上
    sin_h_sun_n = math.sin(phi) * math.sin(delta_n) + math.cos(phi) * math.cos(delta_n) * math.cos(tau_n)
    cos_h_sun_n = math.sqrt(1.0 - sin_h_sun_n ** 2)

    # 太陽方位角の正弦 sin a_sun^n と余弦 cos a_sun^n, -
    # 式 (9)・(10)
    sin_a_sun_n = math.cos(delta_n) * math.sin(tau_n) / cos_h_sun_n
    cos_a_sun_n = (math.cos(delta_n) * math.sin(phi) * math.cos(tau_n) - math.sin(delta_n) * math.cos(phi)) / cos_h_sun_n

    return sin_h_sun_n, cos_h_sun_n, sin_a_sun_n, cos_a_sun_n
