"""時間を進めるループ（段階 6、設計の文書 3 章）。

境界条件は、経過時間（秒）を渡すとその時刻の値を返す関数で受け取る（Δt と切り離すため）。
1 ステップを進める部分は solver.step。木材腐朽を求める層があれば、毎ステップ、腐朽の状態も進める（decay）。

計算方法の文書の評価方法との対応：「時間の進め方」（docs/source/contents/03_13_eval_time.md）。
コメントの「式 (k)」は、とくに断りがなければその節の式の番号。
"""

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np

from modules import decay
from modules.air_water_properties import T_0
from modules.balance import drainage_is
from modules.boundary import Boundaries
from modules.newton import NewtonNotConverged
from modules.solver import PoresFilled, step
from modules.state import State
from modules.wall_geometry import FilmNode, WallGeometry
from modules.wall_state import WallState

# 境界条件の関数：計算を始めてからの経過時間, s → その時刻の外気と室内の境界条件
BoundariesAt = Callable[[float], Boundaries]


@dataclass(frozen=True, eq=False)
class SimulationResult:
    """時間を進めた結果。

    記録は最初の状態（経過時間 0）と、record_every ステップごとの状態（添え字 n は記録の番号）。
    質点ごと・記録ごとの配列は、数式の T_i,n と同じ並び（i が先、n があと）で、形は [i, n]。
    """

    # 記録した時刻（計算を始めてからの経過時間）, s, [n]
    time_ns: np.ndarray

    # 各質点の、記録した時刻の温度, K, [i, n]
    t_is_ns: np.ndarray

    # 各質点の、記録した時刻の水分化学ポテンシャル, J/kg, [i, n]
    mu_is_ns: np.ndarray

    # 各ステップのニュートン法の修正の回数, [ステップ]
    iterations_ns: np.ndarray

    # 最後のステップのあとの状態（記録の間隔によらず、続きの計算に使える）
    state_last: WallState

    # 各質点の、記録した時刻の木材腐朽の質量減少率, -, [i, n]（木材腐朽を求めない質点は 0）
    mass_loss_is_ns: np.ndarray | None = None

    # 最後のステップのあとの木材腐朽の状態（木材腐朽を求める層がなければ None）
    decay_last: decay.DecayState | None = None

    # 各質点の、記録した時刻までに流れ落ちた水の累計 D_i,n, kg/m2, [i, n]（水膜の質点のほかは 0）
    drained_is_ns: np.ndarray | None = None

    def state(self, n: int) -> WallState:
        """記録 n の壁の状態。"""
        return WallState(t_is=self.t_is_ns[:, n], mu_is=self.mu_is_ns[:, n])


def simulate(
        geometry: WallGeometry,
        state_0: WallState,
        dt: float,
        n_steps: int,
        boundaries_at: BoundariesAt,
        record_every: int = 1
    ) -> SimulationResult:
    """境界条件に沿って、n_steps ステップ進める。

    ステップ n（0 〜 n_steps - 1）は、時刻 n Δt の状態から、ステップの終わりの時刻 (n + 1) Δt の境界条件
    boundaries_at((n + 1) Δt) を使って（後退差分）、時刻 (n + 1) Δt の状態を求める。

    Args:
        geometry: 壁の形
        state_0: 最初（経過時間 0）の状態
        dt: 時間刻み, s
        n_steps: 進めるステップの数
        boundaries_at: 経過時間, s → その時刻の境界条件
        record_every: 何ステップごとに記録するか（1 なら毎ステップ）

    Returns:
        時間を進めた結果

    Raises:
        Condensation: 通気層の空気の質点で結露して、相対湿度 100% 以下では収支が合わないとき（solver.step）
        PoresFilled: 空隙率を与えた材料の質点の水分量が、空隙を水で満たす水分量を超えたとき（solver.step）
        NewtonNotConverged: それ以外の理由で収束しなかったとき（solver.step）
        どの例外にも、止まるまでの結果を属性 partial_result（SimulationResult。記録は止まる前の最後の記録まで、
        修正の回数と state_last は解けた最後のステップまで）に、解けなかったステップの終わりの時刻 t_n+1（s）を
        属性 time_failed に付けて出す。

    Notes:
        計算方法の文書の「時間の進め方」に対応する。
    """

    # 時間刻みが正でなければエラー
    if not dt > 0.0:
        raise ValueError(f'時間刻み dt は正にしてください（{dt}）。')

    # ステップの数が 0 以上の整数でなければエラー（0 なら最初の状態だけを返す）
    if not isinstance(n_steps, int) or n_steps < 0:
        raise ValueError(f'ステップの数 n_steps は 0 以上の整数にしてください（{n_steps}）。')

    # 記録の間隔が 1 以上の整数でなければエラー
    if not isinstance(record_every, int) or record_every < 1:
        raise ValueError(f'記録の間隔 record_every は 1 以上の整数にしてください（{record_every}）。')

    # 記録の数（最初の状態 ＋ record_every ステップごと）
    n_records = 1 + n_steps // record_every

    time_ns = np.empty(n_records)
    t_is_ns = np.empty((len(state_0), n_records))
    mu_is_ns = np.empty((len(state_0), n_records))
    mass_loss_is_ns = np.zeros((len(state_0), n_records))
    drained_is_ns = np.zeros((len(state_0), n_records))
    iterations_ns = np.empty(n_steps, dtype=int)

    # 最初の状態を記録する
    time_ns[0], t_is_ns[:, 0], mu_is_ns[:, 0] = 0.0, state_0.t_is, state_0.mu_is

    state = state_0

    # 木材腐朽を求める質点と、ステップ n と n - 1 の木材腐朽の状態（はじめは腐朽していない）
    decay_nodes = decay.DecayNodes.from_geometry(geometry)
    decay_n = decay_nmns = decay.DecayState.initial(len(state_0))

    # ステップ 0 の境界条件（通気層の換気量に使う）
    b_n = boundaries_at(0.0)

    # 水膜の質点があるか（なければ、流れ落ちる水はない）と、ステップ n までに流れ落ちた水の累計 D_i,n, kg/m2, [i]（はじめは 0）
    has_film = any(isinstance(node, FilmNode) for node in geometry.nodes)
    d_is_n = np.zeros(len(state_0))

    def result(n_done: int) -> SimulationResult:
        """ステップ 0 〜 n_done - 1 を解いたところまでの結果。"""
        r_done = 1 + n_done // record_every
        return SimulationResult(
            time_ns=time_ns[:r_done],
            t_is_ns=t_is_ns[:, :r_done],
            mu_is_ns=mu_is_ns[:, :r_done],
            iterations_ns=iterations_ns[:n_done],
            state_last=state,
            mass_loss_is_ns=mass_loss_is_ns[:, :r_done],
            decay_last=decay_n if decay_nodes.any else None,
            drained_is_ns=drained_is_ns[:, :r_done],
        )

    # 解けなかったら、止まるまでの結果と、解けなかったステップの終わりの時刻を例外に付けて出し直す
    try:
        for n in range(n_steps):

            # ステップ n + 1 の時刻 t_n+1, s
            # 式 (1)
            t_npls = (n + 1) * dt

            # ステップ n + 1 の境界条件で、1 ステップ進める（換気量だけは、ステップ n の外気で求める）
            # 4 節
            b_npls = boundaries_at(t_npls)

            # 木材腐朽を求める層がなければ、腐朽で生まれる水はなし
            if not decay_nodes.any:
                state, iterations_ns[n] = step(geometry, state, dt, b_n.outdoor, b_npls.outdoor, b_npls.indoor)

            # 木材腐朽を求める層があれば、ステップ n - 1 から n までの腐朽で生まれた水を足して解き、
            # 求めたステップ n + 1 の温度・相対湿度で、腐朽の状態を進める
            # 4 節
            else:
                j_gen_is_npls = decay.water_generation_is(decay_nodes, decay_n, decay_nmns, dt)
                state, iterations_ns[n] = step(geometry, state, dt, b_n.outdoor, b_npls.outdoor, b_npls.indoor, j_gen_is_npls)
                rh_is_npls = np.array([State(t=t, mu=min(mu, 0.0)).rh for t, mu in zip(state.t_is, state.mu_is)])
                decay_nmns, decay_n = decay_n, decay.advance(decay_n, decay_nodes, state.t_is - T_0, rh_is_npls, dt)

            # ステップ n + 1 までに流れ落ちた水の累計 D_i,n+1, kg/m2, [i]
            # 式 (2)。流れ落ちる水は、収支と同じくステップ n + 1 の状態で求める
            if has_film:
                d_is_n = d_is_n + drainage_is(geometry, state.t_is, state.mu_is) * dt

            # 次のステップの始めの境界条件は、このステップの終わりの境界条件
            b_n = b_npls

            # record_every ステップごとに記録する
            # 5 節
            if (n + 1) % record_every == 0:
                r = (n + 1) // record_every
                time_ns[r], t_is_ns[:, r], mu_is_ns[:, r] = t_npls, state.t_is, state.mu_is
                mass_loss_is_ns[:, r] = decay_n.mass_loss_is
                drained_is_ns[:, r] = d_is_n
    except (NewtonNotConverged, PoresFilled) as e:
        e.partial_result = result(n)
        e.time_failed = t_npls
        raise

    return result(n_steps)
