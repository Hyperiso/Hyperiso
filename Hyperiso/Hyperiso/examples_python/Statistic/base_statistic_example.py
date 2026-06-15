from pathlib import Path

from pyhyperiso.Common import (
    Model,
    Observables,
    QCDOrder,
    ParameterType,
    Decays,
    BinnedObservableId,
    ParamId
)

from pyhyperiso.Core import HyperisoConfig, HyperisoMaster, ExternalFlag, ParameterProvider

from pyhyperiso.Observable import ObservableInterface

from pyhyperiso.Statistic import StatisticInterface, StatisticConfig, ContourOptions, ContourAlgorithm

if __name__ == "__main__":

    config = HyperisoConfig(
        model=Model.SM,
    )
    hyp = HyperisoMaster()
    lha_file_path = "lha/si_input.flha" 
    hyp.init(lha_file=lha_file_path, config=config)
    
    oi = ObservableInterface()
    
    oi.add_observables_from_decay(Decays.B__l_l, QCDOrder.NNLO)
    
    #Configuration for the stastitical interfaces.
    #A lot of options can be fine-tuned but we advice the user only to change few of them.
    sc = StatisticConfig()
    
    #The most useful parameter to change is the number of draw for the MonteCarlo generator.
    sc.MC_draws = 100
    
    #Main interface for uncertainty calculation, MLE and making 2D contours.
    si = StatisticInterface(sc, oi)
    
    #This API calculate symmetric and asymetric uncertaities of the observables given in the ObservableInterface.
    #It calculate the skew as well to know if the distribution is symmetric or not.
    result_u = si.compute_uncertainties()
    
    print(result_u)
    print("\n\n")
    #One can also fit over one or multiple parameter through the following API. It return the best values for the p_specs and the nuisances.
    result_MLE = si.compute_MLE([ParamId(ParameterType.FLAVOR, "FCONST", [511, 1]), ParamId(ParameterType.FLAVOR, "FCONST", [531, 1])])
    
    print(result_MLE)
    print("\n\n")
    
    co = ContourOptions()
    
    result_contour = si.compute_confidence_contour(ParamId(ParameterType.FLAVOR, "FCONST", [511, 1]), ParamId(ParameterType.FLAVOR, "FCONST", [531, 1]), 1, [-0.5, 0.5, -0.5, 0.5], co)
    
    print(result_contour)
    
    

