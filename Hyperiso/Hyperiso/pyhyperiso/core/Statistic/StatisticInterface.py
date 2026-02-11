from __future__ import annotations

from typing import Any, List

from pyhyperiso.core.Statistic.StatisticConfig import StatisticConfig
from pyhyperiso.phyperiso.pyhyperiso.statistic import StatisticInterface as _CppStatisticInterface

from pyhyperiso.core.Statistic.MCResult import MCResult
from pyhyperiso.core.Statistic.GaussianSummary import GaussianSummary


class StatisticInterface:
    """High-level Python interface to the underlying C++ statistical engine.

    This wrapper hides pybind11-bound objects and returns Python-native wrappers.
    """

    def __init__(self, config: StatisticConfig):
        """Create a StatisticInterface.

        Args:
            config: High-level Python statistic configuration.
        """
        self._cpp = _CppStatisticInterface(config.to_cpp())

    def compute_uncertainties(self) -> Dict[Any, GaussianSummary]:
        """Compute uncertainties and return Gaussian summaries.

        Returns:
            List[GaussianSummary]: Summary objects (Python wrappers).
        """
        cpp_res = self._cpp.compute_uncertainties()

        return{x:GaussianSummary.from_cpp(y) for x,y in cpp_res.items()}

    def compute_uncertainties_and_sampling(self) -> MCResult:
        """Compute uncertainties and also return a Monte Carlo realization.

        Returns:
            MCResult: Python wrapper containing the MC realization and summaries.
        """
        cpp_res = self._cpp.compute_uncertainties_and_sampling()
        return MCResult.from_cpp(cpp_res)
    

if __name__ == "__main__":
    
    from pyhyperiso.core.Core.HyperisoMaster import PyHyperisoMaster
    from pyhyperiso.core.Core.HyperisoConfig import PyHyperisoConfig, ExternalFlag
    from pyhyperiso.core.Common.GeneralEnum import Model, Observables, QCDOrder, ParameterType
    from pyhyperiso.core.Common.ParamId import ParamId
    from pathlib import Path
    config = PyHyperisoConfig(
        flags={
            ExternalFlag.IS_LHA_SPECTRUM: False,
            ExternalFlag.HAS_WILSON_INPUT: False,
            ExternalFlag.HAS_TH_OBSERVABLE_INPUT: False,
            ExternalFlag.HYP_AS_SM_MARTY: True
        },
        model=Model.SM,
        mty_model_name="MSSM_UFO",
        mty_model_path=Path("/my/custom/marty/path")
    )
    
    hyp = PyHyperisoMaster()
    lha_file_path = "/home/theo/hyperiso/Assets/lha/si_input.flha" 
    
    hyp.init(lha_file=lha_file_path, config=config)
    
    obs = {Observables.BR_BS_MUMU : QCDOrder.LO, Observables.BR_BS_MUMU_UNTAG : QCDOrder.LO, Observables.BR_BD_MUMU : QCDOrder.LO}
    config_stat = StatisticConfig(obss = obs)
    config_stat.p_specs = [ParamId(ParameterType.DECAY, "B_ll", 1)]
    config_stat.MC_draws = 10000
    si = StatisticInterface(config_stat)
    
    print(si.compute_uncertainties())
    
    # print(si.compute_uncertainties_and_sampling())
    
    # print(len(si.compute_uncertainties_and_sampling().mc_real.sampled_obss))
    # print(len(si.compute_uncertainties_and_sampling().mc_real.sampled_params))
    
    
    
    import matplotlib.pyplot as plt
    from collections import defaultdict

    data = si.compute_uncertainties_and_sampling().mc_real.sampled_obss

    values_by_obs = defaultdict(list)

    for d in data:
        for obs, val in d.items():
            values_by_obs[str(obs)].append(val)  

    for obs, values in values_by_obs.items():
        plt.figure()
        plt.hist(values, bins=50)
        plt.xlabel(obs)
        plt.ylabel("Counts")
        plt.title(f"Distribution of {obs}")
        plt.tight_layout()
        plt.show()