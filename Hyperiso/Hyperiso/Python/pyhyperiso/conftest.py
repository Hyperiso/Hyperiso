from Python.Phyperiso import MemoryManager
from Python.Phyperiso import Parameters, QCDHelper
from Python.Phyperiso import WilsonInterface
from Python.Phyperiso import Model, WGroup, WCoeff, ParameterType, QCDOrder

def main():
    mm = MemoryManager()
    mm.init("Test/InputFiles/testInput.flha", Model.SM) # Initialize program manager with LHA file

    wi = WilsonInterface() # Initialize interface and build the required groups
    wi.build(
        [WGroup.B.value, WGroup.BPrime.value],                     # Coefficient groups
        2 * Parameters(ParameterType.SM)("MASS", 24),  # Matching scale
        QCDHelper.mass_b_1S() / 2,                     # Hadronic scale
        QCDOrder.NNLO                                  # QCD Order
    )
    
    # Retrieve coefficient values
    with open("out.dat", "w") as f:
        f.write("C7 matching (LO+NLO+NNLO) : "  + str(wi.get_full_matching_coefficient(WGroup.B, WCoeff.C7, QCDOrder.NNLO).real) + "\n")
        f.write("C7 hardronic (LO+NLO+NNLO) : " + str(wi.get_full_run_coefficient(WGroup.B, WCoeff.C7, QCDOrder.NNLO).real) + "\n")
        f.write("CP7 matching (LO) : "          + str(wi.get_full_matching_coefficient(WGroup.BPrime, WCoeff.CP7, QCDOrder.LO).real) + "\n")
        f.write("CP7 hardronic (LO) : "         + str(wi.get_full_run_coefficient(WGroup.BPrime, WCoeff.CP7, QCDOrder.LO).real) + "\n")

    return 0

if __name__ == "__main__":
    main()