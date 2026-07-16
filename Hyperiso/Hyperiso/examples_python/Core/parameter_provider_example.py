from pyhyperiso.Common import ParameterType, ParamId, DataType, CorrelationType, Observables

from pyhyperiso.Core import (
    HyperisoConfig,
    HyperisoMaster,
    ParameterProvider,
    CorrelationProvider,
    APIAdapter,
    APIPath,
)

if __name__ == "__main__":
    # All informations for the creation of the HyperisoMaster are provided in the base_core_example.py example.
    config = HyperisoConfig()
    hyp = HyperisoMaster()
    lha_file_path = "lha/si_input.flha"
    hyp.init(lha_file=lha_file_path, config=config)

    # This create a parameterprovider for SM parameters. This allows to get parameters without specify each time the parameterType
    sm_param_provider = ParameterProvider(ParameterType.SM)

    # If you want parameters from differents sections, use a no-typed ParameterProvider:
    no_pt_param_provider = ParameterProvider()

    # Using the sm_param_provider, one can test if a parameter exist using the following API:
    print(sm_param_provider.exists_by_block("SMINPUTS", 2))

    # Using the no_pt_param_provider, one can test if a parameter exist using the following API:
    print(no_pt_param_provider.exists_by_pid(ParamId(ParameterType.SM, "SMINPUTS", 2)))

    # One can get the value of a Parameter through the following API. If one want to get the uncertainty (STAT, SYS or COMBINE), one can change the DataType.
    # Return a Scalar (real,img)
    print(sm_param_provider.get_by_block("MASS", 25, DataType.VALUE))

    # Same with a no pt provider.
    # Return a Scalar (real,img)
    print(no_pt_param_provider.get_by_pid(ParamId(ParameterType.SM, "MASS", 25), DataType.VALUE))

    # One can get a copy of a Parameter (class) with all informations (value, id, uncertainties) through :
    param = no_pt_param_provider.get_parameter(ParamId(ParameterType.WILSON, "B_SCALE", 1))
    print(param.value)

    # One might also want to get correlation between parameters or observables.
    # This is useful for the statistic part
    cp = CorrelationProvider()

    # The following API can be use to get correlation between parameters.
    print(
        cp.correlation_from_paramid(
            ParamId(ParameterType.DECAY, "B_K", "1_1_0"),
            ParamId(ParameterType.DECAY, "B_K", "1_1_1"),
            CorrelationType.STAT,
        )
    )

    # For observables, One can use (The last argument is the name of the experiment, default value is "DEFAULT" when its neither CMS or LHCB, Belle, etc.):
    print(
        cp.correlation_from_observable(
            Observables.F_L_B0__KSTAR0_MU_MU,
            Observables.P_1_B0__KSTAR0_MU_MU,
            CorrelationType.STAT,
            "CMS",
        )
    )

    # If one doesn't know in which section a block is in, one can use the following Proxy :

    apiadapt = APIAdapter()

    # This will return a vector of possible ParameterType the block is in.
    # Usually a block is in one section only, but for special cases, like MASS, the parameters can be either in BSM or SM section.
    print(apiadapt.get_type_of_block("SMINPUTS"))

    # This proxy is also useful to get blocks as list of parameters :
    print(apiadapt.get_block_infos("MASS", ParameterType.SM))

    # It can also be use to get informations on paths :
    print(apiadapt.get_path(APIPath.LHA_PATH))
