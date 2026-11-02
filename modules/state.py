from dataclasses import dataclass
import numpy as np

from modules.air_water_properties import T_0, get_mu, get_p_vs, get_rh


@dataclass
class State:

    # 絶対温度, K
    t: float

    # 水分化学ポテンシャル, J/kg
    mu: float

    @classmethod
    def init(cls, t: float, rh: float):
        """Stateクラスの初期化を行うクラスメソッド

        Args:
            t: 絶対温度, K
            rh: 相対湿度, %
        
        Returns:
            Stateクラスのインスタンス
        """

        # 水分化学ポテンシャルの初期値を計算
        wp = get_mu(rh=rh, t=t)

        return State(t=t, mu=wp)

    @property
    def theta(self) -> float:
        return self.t - T_0
    
    @theta.setter
    def theta(self, value):
        self.t = value + T_0
    
    @property
    def p_v_sat(self) -> float:
        """飽和水蒸気圧, Pa（「湿り空気と水の性質」の式 (1)）"""

        return get_p_vs(self.t)
    
    @property
    def rh(self) -> float:
        """相対湿度, %（「湿り空気と水の性質」の式 (3)）"""

        return get_rh(mu=self.mu, t=self.t)
    
    @property
    def p_v(self) -> float:
        """水蒸気圧, Pa（「湿り空気と水の性質」の式 (4)）"""

        return self.rh / 100.0 * self.p_v_sat
    
