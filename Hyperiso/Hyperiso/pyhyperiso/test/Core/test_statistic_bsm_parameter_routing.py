from pathlib import Path


ROOT = Path(__file__).resolve().parents[5]


def source(relative: str) -> str:
    return (ROOT / relative).read_text(encoding="utf-8")


def test_statistic_parameter_mutations_preserve_bsm_scope():
    proxy_cpp = source(
        "Hyperiso/Hyperiso/core/src/Statistic/adapters/"
        "StatParamOptimizerProxy.cpp"
    )
    proxy_h = source(
        "Hyperiso/Hyperiso/core/src/Statistic/adapters/"
        "StatParamOptimizerProxy.h"
    )
    observable_proxy = source(
        "Hyperiso/Hyperiso/core/src/Statistic/adapters/"
        "ObservableInterfaceProxy.cpp"
    )
    port = source(
        "Hyperiso/Hyperiso/core/src/Statistic/ports/"
        "IStatParamOptimizerProxy.h"
    )

    assert "poa_bsm({ParameterType::BSM})" in proxy_cpp
    assert "pid.type.value() == ParameterType::BSM" in proxy_cpp
    assert "void set_value(const ParamId& pid" in proxy_h
    assert "virtual void set_value(const ParamId& pid" in port
    assert "spop_->set_value(pid, value);" in observable_proxy

    assert "spop_->set_value(s.block, s.code" not in observable_proxy
