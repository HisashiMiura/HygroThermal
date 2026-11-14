import math

import pytest

from modules import air_water_properties as aw


class TestPVs:

    def test_20c(self):
        # 20 deg.C の飽和水蒸気圧は、およそ 2339 Pa
        assert aw.get_p_vs(20.0 + aw.T_0) == pytest.approx(2339.0, rel=1e-3)


class TestMuRh:

    def test_saturation(self):
        # 相対湿度 100% では μ = 0
        assert aw.get_mu(100.0, 293.15) == 0.0

    def test_round_trip(self):
        # 式 (2) と式 (3) は互いに逆
        for rh in (10.0, 50.0, 99.0):
            assert aw.get_rh(aw.get_mu(rh, 293.15), 293.15) == pytest.approx(rh, rel=1e-12)

    def test_value(self):
        # μ = R_v T ln(φ / 100)
        assert aw.get_mu(50.0, 293.15) == pytest.approx(461.5 * 293.15 * math.log(0.5))

    def test_rh_upper_limit(self):
        # μ > 0 でも相対湿度は 100% を超えない
        assert aw.get_rh(1000.0, 293.15) == 100.0


class TestRhoA:

    def test_value(self):
        # ρ_a = 353 / T
        assert aw.get_rho_a(293.15) == pytest.approx(353.0 / 293.15)
