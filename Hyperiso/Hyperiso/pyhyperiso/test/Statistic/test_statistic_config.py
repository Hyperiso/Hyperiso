import pytest

from pyhyperiso.core.Statistic.StatisticConfig import AdvancedStatisticConfig


def test_fit_parameter_sensitivity_defaults_to_cpp():
    cfg = AdvancedStatisticConfig()
    cpp = cfg.to_cpp()

    assert cpp.fit_parameter_sensitivity_check is True
    assert cpp.fit_parameter_sensitivity_probe_fraction == pytest.approx(0.05)
    assert cpp.fit_parameter_sensitivity_rel_cutoff == pytest.approx(1e-10)
    assert cpp.fit_parameter_sensitivity_abs_cutoff == pytest.approx(1e-12)
    assert cpp.fit_parameter_sensitivity_keep_on_failure is True


def test_fit_parameter_sensitivity_custom_values_to_cpp():
    cfg = AdvancedStatisticConfig(
        fit_parameter_sensitivity_check=False,
        fit_parameter_sensitivity_probe_fraction=0.1,
        fit_parameter_sensitivity_rel_cutoff=2e-9,
        fit_parameter_sensitivity_abs_cutoff=3e-11,
        fit_parameter_sensitivity_keep_on_failure=False,
    )
    cpp = cfg.to_cpp()

    assert cpp.fit_parameter_sensitivity_check is False
    assert cpp.fit_parameter_sensitivity_probe_fraction == pytest.approx(0.1)
    assert cpp.fit_parameter_sensitivity_rel_cutoff == pytest.approx(2e-9)
    assert cpp.fit_parameter_sensitivity_abs_cutoff == pytest.approx(3e-11)
    assert cpp.fit_parameter_sensitivity_keep_on_failure is False
