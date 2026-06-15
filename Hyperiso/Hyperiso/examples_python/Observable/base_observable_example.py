from pathlib import Path

from pyhyperiso.Common import (
    Model,
    Observables,
    QCDOrder,
    ParameterType,
    Decays,
    BinnedObservableId
)

from pyhyperiso.Core import HyperisoConfig, HyperisoMaster, ExternalFlag, ParameterProvider

from pyhyperiso.Observable import ObservableInterface
from pyhyperiso.Observable import BKstarllConfig, BV_FF_Src

if __name__ == "__main__":

    config = HyperisoConfig(
        model=Model.SM,
    )
    hyp = HyperisoMaster()
    lha_file_path = "lha/si_input.flha" 
    hyp.init(lha_file=lha_file_path, config=config)
    

    #Creation of the Observable interface. This interface need to be created after the HyperisoMaster, and will register all observables needed for the calculation.
    interface = ObservableInterface()
    
    #For some decay, you can fine-tuned the config if you want to use a specific set of form factor or change some parameters for the calculation
    #For heavy calculation (e.g Kstarll) with integration, you can use the "kstar_conf.n_threads" option to run on multiple threads.
    #This is also available within the interface with : interface.set_bkstarll_threads(nb_of_thread) (same for bkll and bsphi).
    kstar_conf = BKstarllConfig()
    
    #Use the add_observable to add a new observable, with its order in QCD for the calculation.
    interface.add_observable(Observables.BR_B_XS_GAMMA, QCDOrder.NNLO)
    
    #With bin, use the BinnedObservableId class to add the bin. You can add multiple bin, which will be treated as differents observables.
    interface.add_binned_observable(BinnedObservableId(Observables.F_L_B0__KSTAR0_MU_MU, {1.1, 6.0}), QCDOrder.NNLO)
    
    #If you want to add all the observables of a decay, use the add_observables_from_decay(Decays, QCDOrder) API.
    interface.add_observables_from_decay(Decays.B__l_l, QCDOrder.NNLO)
    
    #Compute the observable. Here bins are not necessary, the API will return a vector of Observablue, which contain id, bin and value.
    print(interface.compute_observable(Observables.F_L_B0__KSTAR0_MU_MU))  # List(ObservableValue(...))
    
    #This allows to change the form factor choice of a particular decay.
    kstar_conf.ff_src = BV_FF_Src.GRvDV 
    
    #Set the new config inside the interface.
    interface.set_decay_config(Decays.B__Kstar_l_l, kstar_conf)
    
    print(interface.compute_observable(Observables.F_L_B0__KSTAR0_MU_MU))  # List(ObservableValue(...))
    print(interface.compute_observable(Observables.BR_BS_MUMU))  # List(ObservableValue(...))
    
    #This API allows to calculate all observables present in the interface (all add through the differents API).
    print(interface.compute_all())
    