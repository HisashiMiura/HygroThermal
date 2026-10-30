"""時間の軸：経過時間の起点と、1 日・1 年の長さ（段階 6）。

経過時間の起点は 1 月 1 日 0 時。1 年は 365 日（うるう年なし）で、2 年目以降も同じ周期で続く。
外気・室内の境界条件は、この経過時間, s から値を作る。
"""

# 1 時間, s
SECONDS_PER_HOUR = 3600.0

# 1 日, s
SECONDS_PER_DAY = 24 * SECONDS_PER_HOUR

# 1 年（365 日）, s
SECONDS_PER_YEAR = 365 * SECONDS_PER_DAY

# 1 年の時間の数（1 時間ごとの気象データの行の数）
HOURS_PER_YEAR = 365 * 24
