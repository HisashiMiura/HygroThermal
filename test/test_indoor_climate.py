import pytest

from modules.boundary import IndoorBoundary
from modules.indoor_climate import CosineIndoor
from modules.time_axis import SECONDS_PER_DAY, SECONDS_PER_YEAR

# 今の Python 版の既定値（平均 22.5 deg.C、振幅 4.5 K、相対湿度 60%）
_INDOOR = CosineIndoor(theta_in_ave=22.5, theta_in_amp=4.5, rh_in=60.0)

# 7 月 31 日（通算日 212）12 時：1 月 1 日 0 時から 211.5 日後
_PEAK = 211.5 * SECONDS_PER_DAY


class TestCosineIndoor:

    def test_peak(self):
        # 7 月 31 日 12 時に、平均 + 振幅
        assert _INDOOR.get_indoor_boundary(_PEAK).state.theta == pytest.approx(27.0)

    def test_bottom(self):
        # 極大の半年後に、平均 - 振幅
        assert _INDOOR.get_indoor_boundary(_PEAK + SECONDS_PER_YEAR / 2).state.theta == pytest.approx(18.0)

    def test_average(self):
        # 極大の 1/4 年前と 1/4 年後は、年平均
        assert _INDOOR.get_indoor_boundary(_PEAK - SECONDS_PER_YEAR / 4).state.theta == pytest.approx(22.5)
        assert _INDOOR.get_indoor_boundary(_PEAK + SECONDS_PER_YEAR / 4).state.theta == pytest.approx(22.5)

    def test_period(self):
        # 365 日後（2 年目の同じ時刻）は同じ値
        time = 40.3 * SECONDS_PER_DAY
        assert _INDOOR.get_indoor_boundary(time + SECONDS_PER_YEAR).state.theta == pytest.approx(_INDOOR.get_indoor_boundary(time).state.theta)

    def test_rh(self):
        # 相対湿度はいつでも一定。返すのは室内の境界条件
        for time in (0.0, _PEAK, 300.0 * SECONDS_PER_DAY):
            indoor = _INDOOR.get_indoor_boundary(time)
            assert isinstance(indoor, IndoorBoundary)
            assert indoor.state.rh == pytest.approx(60.0)

    def test_peak_day(self):
        # いちばん暑い日を変えると、その日の 12 時が極大
        indoor = CosineIndoor(theta_in_ave=20.0, theta_in_amp=5.0, rh_in=50.0, d_peak=1)
        assert indoor.get_indoor_boundary(0.5 * SECONDS_PER_DAY).state.theta == pytest.approx(25.0)

    @pytest.mark.parametrize('kwargs, match', [
        ({'theta_in_amp': -1.0}, 'theta_in_amp'),
        ({'rh_in': 0.0}, 'rh_in'),
        ({'rh_in': 101.0}, 'rh_in'),
        ({'d_peak': 0}, 'd_peak'),
        ({'d_peak': 366}, 'd_peak'),
        ({'d_peak': 212.5}, 'd_peak'),
    ])
    def test_invalid(self, kwargs, match):
        args = {'theta_in_ave': 22.5, 'theta_in_amp': 4.5, 'rh_in': 60.0} | kwargs
        with pytest.raises(ValueError, match=match):
            CosineIndoor(**args)
