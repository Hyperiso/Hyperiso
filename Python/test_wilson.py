
# import ctypes

import phyperiso.core as core
import phyperiso.wilson as wilson
import matplotlib.pyplot as plt
import numpy as np

def test_coefficient_manager():
    model_name = "SM"
    # mm = core.MemoryManager.get_instance("Test/InputFiles/testInput.flha", [0])
    # mm.init()
    # param = core.Parameters.get_instance(0)
    # param = core.Parameters.get_instance(0)
    # param = core.Parameters.get_instance(0)
    # print(param.alpha_s(81))
    range2 = np.arange(50,250, 1)
    manager = wilson.CoefficientManager.get_instance(model_name)
    # param = core.Parameters.get_instance(0)
    wilson.CoefficientManager.initialize("Test/InputFiles/testInput.flha", [0])
    print(model_name)
    group = wilson.BCoefficientGroup()

    group_name = "BCoefficientGroup"
    manager.register_coefficient_group(group_name, group)

    Q_match = 81
    Q = 42
    manager.set_q_match(group_name, Q_match)
    # print("mtop : ", manager.get_params("MASS", 6))
    # manager.set_params("BCoefficientGroup", "MASS", 6, 150)
    # print("new mtop : ", manager.get_params("MASS", 6))
    liste = []
    liste_NLO = []
    liste_NNLO = []
    liste_alpha = []
    order = "LO"
    coeff_name = "C7"
    for mt in range2:
        order = "LO"
        manager.set_params("BCoefficientGroup", "MASS", 6, mt)
        manager.set_matching_coefficient(group_name, order)
        manager.set_group_scale(group_name, Q)
        manager.set_run_coefficient(group_name, order)
        liste_alpha.append(manager.get_alpha_s("BCoefficientGroup"))
        liste.append(manager.get_matching_coefficient(group_name, coeff_name, order).real)

        order = "NLO"
        manager.set_matching_coefficient(group_name, order)
        liste_NLO.append(manager.get_matching_coefficient(group_name, coeff_name, order).real)

        order = "NNLO"
        manager.set_matching_coefficient(group_name, order)
        liste_NNLO.append(manager.get_matching_coefficient(group_name, coeff_name, order).real)



    LO = liste[1:]
    NLO = np.array(liste_NLO[1:]) *np.array( liste_alpha[1:]) / (2* np.pi)
    NNLO =  liste_NNLO[1:] *np.array( liste_alpha[1:])*np.array( liste_alpha[1:]) / (2* np.pi)**2
    plt.plot(range2[1:], LO, color = "gold", label = "LO")
    plt.plot(range2[1:], LO +NLO, color = "orange", label = "LO + NLO")
    plt.plot(range2[1:], LO + NLO +NNLO, color = "red", label = "LO + NLO + NNLO")
    plt.grid(True)
    plt.legend()
    plt.xlabel(r"$m_{t,pole}$ [GeV]")
    plt.ylabel(r"C7")
    plt.tick_params(direction="in", left = True, right = True, top = True, bottom = True)
    plt.xlim(50,250)
    plt.title(r"C7 Coefficient as a function of $m_{t,pole}$")
    plt.savefig("C7_coeff", dpi = 800)
    plt.show()
    # manager.set_group_scale(group_name, Q)
    # manager.set_run_coefficient(group_name, "LO")
    
    
    # try:
    #     matching_value = manager.get_matching_coefficient(group_name, coeff_name, order)
    #     run_value = manager.get_run_coefficient(group_name, coeff_name, order)
    #     print(f"Matching value for {coeff_name} at Q_match = {Q_match}: {matching_value}")
    #     print(f"Running value for {coeff_name} at Q = {Q}: {run_value}")
    # except Exception as e:
    #     print(f"Error during coefficient retrieval: {e}")

if __name__ == "__main__":
    test_coefficient_manager()