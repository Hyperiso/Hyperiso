from pyhyperiso.core.Statistic.StatisticConfig import StatisticConfig
from pyhyperiso.phyperiso.pyhyperiso.statistic import StatisticInterface as _CppStatisticInterface
from pyhyperiso.core.Common.General import PyParamId
from pyhyperiso.core.Common.GeneralEnum import ParameterType

class StatisticInterface:
    def __init__(self, config : StatisticConfig):
        self.si = _CppStatisticInterface(config.to_cpp())
    
    def compute_uncertainties(self):
        return self.si.compute_uncertainties()
    
    def compute_MLE(self):
        return self.si.compute_MLE()
    

if __name__ == "__main__":
    
    from pyhyperiso.core.Core.HyperisoMaster import PyHyperisoMaster
    from pyhyperiso.core.Core.HyperisoConfig import PyHyperisoConfig, ExternalFlag
    from pyhyperiso.core.Common.GeneralEnum import Model, Observables, QCDOrder
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
    config_stat.p_specs = [PyParamId(ParameterType.DECAY, "B_ll", 1)]

    si = StatisticInterface(config_stat)
    
    print(si.compute_uncertainties())
    
    