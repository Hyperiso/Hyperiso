from pathlib import Path

from pyhyperiso.Common import (
    Model,
    WGroup,
    WCoeff,
    QCDOrder,
    ParameterType,
    ContributionType,
    ScaleType,
    WilsonBasis
    )

from pyhyperiso.Core import (
    HyperisoConfig, 
    HyperisoMaster, 
    ExternalFlag, 
    ParameterProvider,
    BlockLogger
)

from pyhyperiso.Wilson import WilsonBuildConfig, WilsonInterface, WilsonRequest

if __name__ == "__main__":

    #See base_core_example for informations about these classes.
    config = HyperisoConfig(
        flags={
            ExternalFlag.IS_LHA_SPECTRUM: False, 
            ExternalFlag.HAS_WILSON_INPUT: False,
            ExternalFlag.HAS_TH_OBSERVABLE_INPUT: False,
            ExternalFlag.HYP_AS_SM_MARTY: True,
        },
        model=Model.SM,
        mty_model_name="marty_model_name",
        mty_model_path=Path("/my/custom/marty/path")
    )

    hyp = HyperisoMaster()
    lha_file_path = "lha/si_input.flha" 

    hyp.init(lha_file=lha_file_path, config=config)
    
    #Build configuration for wilson coefficients, with matching/running scales information, name of the Wilson group and order in QCD.
    config = WilsonBuildConfig(
        groups={WGroup.B, WGroup.BScalar}, #Groups to be build.
        matching_scale=81.0, #Matching scale. This scale is global, meaning all group will have the same matching scale.
        hadronic_scale=2.0, #Running scale. This scale is global like the matching scale. Specials scales for K mesons/D mesons coefficients or others need to be modified within Hyperiso with the ParameterSetter API (D_SCALE, K_SCCALE).
        order=QCDOrder.LO #Order in QCD. Careful when using this with MARTY. SM contribution can go up to NNLO but the BSM part will stay at LO.
    )


    interface = WilsonInterface() #Wilson Interface for building and requesting wilson coefficients.
    
    interface.build(config) #API to build the wilson coefficients specifed within the configuration.


    config.groups = [WGroup.BPrime]
    interface.add_wilson_group(config) #You can add a new group after the build. The same QCDOrder and scales will be used.
     
    #Requests to the interface for wilson coefficients.
    #group : name of the group the coefficient is in.
    #coefficient : name of the coefficient (need to be inside the group).
    #order : QCD order at which the coefficient is requested. If above the QCDOrder used to build, the coefficient will be 0.
    #contribution : SM, BSM or TOTAL = SM + BSM.
    #scale_type : If HADRONIC the interface will return the running coefficients, if MATCHING, then it will be the matching coefficient.
    #wilson_basis : Basis used for the running of the coefficients. Use STANDARD by default, only TRADITIONAL for B group.
    req = WilsonRequest(
        group=WGroup.B,
        coefficient=WCoeff.C9,
        order=QCDOrder.NNLO,
        contribution=ContributionType.TOTAL,
        scale_type=ScaleType.HADRONIC,
        wilson_basis=WilsonBasis.STANDARD
    )
    coefs = {WCoeff.C1, WCoeff.C2, WCoeff.C3, WCoeff.C4, WCoeff.C5, WCoeff.C6, WCoeff.C7, WCoeff.C8, WCoeff.C9, WCoeff.C10}
    coefs_primes = {WCoeff.CP1, WCoeff.CP2, WCoeff.CP3, WCoeff.CP4, WCoeff.CP5, WCoeff.CP6, WCoeff.CP7, WCoeff.CP8, WCoeff.CP9, WCoeff.CP10, WCoeff.CPQ1, WCoeff.CPQ2}
    coefs_scalar = {WCoeff.CQ1, WCoeff.CQ2}
    
    for coef in coefs:
        req.coefficient = coef
        print(req.coefficient, " : ", interface.get_sep_order_matching(req.group, req.coefficient, req.contribution)) #Get a dictionnary with all requested order and their contributions.
    
    print("\n\n\n")

    req.group = WGroup.BPrime
    
    for coef in coefs_primes:
        req.coefficient = coef
        print(req.coefficient, " : ", interface.get_FR(req)) #Get full running, with sum of all order up to the requested order.

    print("\n\n\n")

    print(interface.get_all_matching(WGroup.BScalar, QCDOrder.LO, ContributionType.TOTAL))

    