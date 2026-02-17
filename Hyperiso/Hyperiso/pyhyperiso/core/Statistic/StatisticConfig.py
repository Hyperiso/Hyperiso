from typing import Dict, List
from dataclasses import dataclass, field
from pyhyperiso.core.Common.GeneralEnum import QCDOrder, Observables
from pyhyperiso.core.Common.ParamId import ParamId
from pyhyperiso.core.Common.Mapper import ObservableMapper
from pyhyperiso.phyperiso.pyhyperiso.statistic import StatisticConfig as _CppStatisticConfig

@dataclass
class StatisticConfig:
    obss : Dict[Observables, QCDOrder]= field(default_factory=lambda: {})
    p_specs : List[ParamId] = field(default_factory=list)
    MC_draws : int = 100
    skew_abs_threshold : float = 0.2
    
    MLE_max_iter : int = 100
    MLE_tol : float = 1e-6
    
    def to_cpp(self) -> _CppStatisticConfig:
        """Converts this Python wrapper to a C++ Config object."""
        cpp = _CppStatisticConfig()
        cpp.obss = {ObservableMapper.to_id(x)._to_cpp():y.value for x,y in self.obss.items()}
        cpp.p_specs = {p._cpp_obj for p in self.p_specs}
        cpp.MC_draws = self.MC_draws
        cpp.skew_abs_threshold = self.skew_abs_threshold
        cpp.MLE_max_iter = self.MLE_max_iter
        cpp.MLE_tol = self.MLE_tol

        return cpp
    
    def __repr__(self):
        return f"StatisticConfig(obss={self.obss}, p_specs={self.p_specs}, MC_draws={self.MC_draws}, skew_abs_threshold = {self.skew_abs_threshold}, MLE_max_iter = {self.MLE_max_iter}, MLE_tot = {self.MLE_tol})"
    
    
if __name__ == "__main__":
    a = StatisticConfig()
    print(a)
    
    a.to_cpp()
    
    import pyhyperiso.phyperiso.pyhyperiso.common as cA
    import pyhyperiso.phyperiso.pyhyperiso.common as cB

    print("A:", cA.__file__)
    print("B:", cB.__file__)
    print("same QCDOrder object?", cA.QCDOrder is cB.QCDOrder)
    print("same LO?", cA.QCDOrder.LO == cB.QCDOrder.LO)