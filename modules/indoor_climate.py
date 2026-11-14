"""室内の条件：ステップ n の室内の境界条件を作る（計算方法の文書の評価方法「室内の条件」）。

各値のコメントは、1 行目に「名前 記号, 単位」（無次元は -）、2 行目に「式 (k)」と補足を書く。
コメントの「式 (k)」は、計算方法の文書の「室内の条件」（docs/source/contents/03_08_eval_indoor.md）の式の番号。
変数の名前は、その節の記号に合わせる。ステップ n の値には、記号と同じく _n を付ける（例：θ_in,n は theta_in_n）。
経過時間の起点は 1 月 1 日 0 時、1 年は 365 日（time_axis.py）。
"""

from dataclasses import dataclass

import numpy as np

from modules.boundary import IndoorBoundary
from modules.time_axis import SECONDS_PER_DAY, SECONDS_PER_YEAR


@dataclass(frozen=True, kw_only=True)
class CosineIndoor:
    """室温が 1 年周期の余弦で変わり、相対湿度は一定の室内。

    いちばん暑い時刻 t_peak は、いちばん暑い日（通算日 d_peak）の 12 時。
    Fortran の SUBROUTINE ROOM は 1 日の中で一定（日ごとの段）で、その日の値が極大になる。ここでは滑らかに変え、日の真ん中を極大にする。
    """

    # 年平均の室温 θ_in,ave, deg.C
    theta_in_ave: float

    # 室温の振幅 θ_in,amp, K
    theta_in_amp: float

    # 室内の相対湿度 φ_in, %（一定）
    rh_in: float

    # 室温がいちばん高い日の通算日 d_peak, d（1 月 1 日が 1）
    d_peak: int = 212

    def __post_init__(self):
        """作るときに値を確かめる。"""

        # 振幅が負ならエラー（負の振幅は、極大と極小が入れ替わるだけで紛らわしい）
        if self.theta_in_amp < 0.0:
            raise ValueError(f'室温の振幅 theta_in_amp は 0 以上にしてください（{self.theta_in_amp}）。')

        # 相対湿度が 0 より大きく 100 以下でなければエラー
        if not 0.0 < self.rh_in <= 100.0:
            raise ValueError(f'相対湿度 rh_in は 0 より大きく 100 以下にしてください（{self.rh_in}）。')

        # いちばん暑い日が 1〜365 の整数でなければエラー
        if not isinstance(self.d_peak, int) or not 1 <= self.d_peak <= 365:
            raise ValueError(f'いちばん暑い日 d_peak は 1〜365 の整数にしてください（{self.d_peak}）。')

    def get_indoor_boundary(self, t_n: float) -> IndoorBoundary:
        """ステップ n の室内の境界条件。

        Args:
            t_n: ステップ n の時刻 t_n（1 月 1 日 0 時からの経過時間）, s

        Returns:
            ステップ n の室内の境界条件
        """

        # 室温がいちばん高い時刻 t_peak, s
        # 式 (1)。通算日 d_peak の 12 時
        t_peak = (self.d_peak - 1 + 0.5) * SECONDS_PER_DAY

        # 室温 θ_in,n, deg.C
        # 式 (2)。余弦は 1 年で元に戻るので、t_n を 1 年で割った余りはとらない
        theta_in_n = self.theta_in_ave + self.theta_in_amp * np.cos(2.0 * np.pi * (t_n - t_peak) / SECONDS_PER_YEAR)

        # 室内の相対湿度 φ_in,n, %
        # 式 (3)
        rh_in_n = self.rh_in

        return IndoorBoundary.from_theta_rh(theta=theta_in_n, rh=rh_in_n)
