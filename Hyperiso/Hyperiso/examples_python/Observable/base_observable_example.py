from pathlib import Path

from pyhyperiso.Common import (
    Model,
    Observables,
    QCDOrder,
    ParameterType,
    Decays
)

from pyhyperiso.Core import HyperisoConfig, HyperisoMaster, ExternalFlag, ParameterProvider

from pyhyperiso.Observable import ObservableInterface

if __name__ == "__main__":
    print("🔧 Initializing PyHyperisoMaster with custom PyHyperisoConfig...")

    config = HyperisoConfig(
        model=Model.SM,
    )
    
    hyp = HyperisoMaster()
    lha_file_path = "lha/si_input.flha" 

    print("\n🚀 Calling init with config...")
    hyp.init(lha_file=lha_file_path, config=config)
    

    from pyhyperiso.core.BusinessLogic.DecayConfig import BKstarllConfig
    
    kstar_conf = BKstarllConfig()
    
    interface = ObservableInterface()
    interface.add_observable(Observables.BR_B_XS_GAMMA, QCDOrder.NNLO, True)
    interface.add_observable(Observables.BR_BS_MUMU, QCDOrder.NNLO, True)
    interface.add_observable(Observables.BR_BD_MUMU, QCDOrder.NNLO, True)
    interface.add_observable(Observables.BR_B0__D_TAU_NU, QCDOrder.NNLO, True)
    interface.add_observable(Observables.BR_B__KSTAR_GAMMA, QCDOrder.NNLO, True)
    interface.set_decay_config(Decays.B__Kstar_l_l, kstar_conf)
    print(interface.compute_observable(Observables.BR_B_XS_GAMMA))  # Scalar(...)
    print(interface.compute_observable(Observables.BR_BS_MUMU))  # Scalar(...)
    print(interface.compute_all())
    
    test_values = []
    
    from pyhyperiso.core.Core.ParameterSetter import ParameterSetter, ParamId
    py_set = ParameterSetter()
    for i in range(50, 90, 10):
        py_set.mutate(ParamId(ParameterType.SM, "MASS", 24), i)
        print("aah : ", ParameterProvider().get_by_pid(ParamId(ParameterType.SM, "MASS", 24)))
        # interface.enable_obs()
        print("aah : ", ParameterProvider().get_by_pid(ParamId(ParameterType.SM, "MASS", 24)))
        test_values.append(interface.compute_observable(Observables.BR_B_XS_GAMMA)[0].value)
    
    import matplotlib
    matplotlib.use("TkAgg")
    import matplotlib.pyplot as plt
    
    plt.scatter(range(50, 90, 10), test_values)
    
    plt.show()
    # test_values = []
    # from pyhyperiso.core.Core.ParameterSetter import PyParameterSetter, ParamId, ParameterType
    # py_set = PyParameterSetter()
    # for i in range(1, 81):
    #     py_set.mutate(ParamId(ParameterType.WILSON, "B_SCALE", 1), i)
    #     test_values.append(interface.get_R(WilsonRequest(WGroup.B, WCoeff.C9, QCDOrder.LO, ContributionType.TOTAL)))
    
    # import matplotlib
    # matplotlib.use("TkAgg")
    # import matplotlib.pyplot as plt
    
    # plt.scatter(range(1, 81), test_values)
    
    # plt.show()