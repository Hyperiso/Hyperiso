
# import ctypes

import phyperiso.core as core
import phyperiso.wilson as wilson

def test_coefficient_manager():
    model_name = "SM"
    # mm = core.MemoryManager.get_instance("Test/InputFiles/testInput.flha", [0])
    # mm.init()
    # param = core.Parameters.get_instance(0)
    # param = core.Parameters.get_instance(0)
    # param = core.Parameters.get_instance(0)
    # print(param.alpha_s(81))
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
    print("mtop : ", manager.get_params("MASS", 6))
    manager.set_params("BCoefficientGroup", "MASS", 6, 150)
    print("new mtop : ", manager.get_params("MASS", 6))
    liste = []

    order = "LO"
    coeff_name = "C7"
    for mt in range(50,200):
        manager.set_params("BCoefficientGroup", "MASS", 6, mt)
        manager.set_matching_coefficient(group_name, "LO")
        liste.append(manager.get_matching_coefficient(group_name, coeff_name, order))

    print(liste)
    manager.set_group_scale(group_name, Q)
    manager.set_run_coefficient(group_name, "LO")
    
    
    try:
        matching_value = manager.get_matching_coefficient(group_name, coeff_name, order)
        run_value = manager.get_run_coefficient(group_name, coeff_name, order)
        print(f"Matching value for {coeff_name} at Q_match = {Q_match}: {matching_value}")
        print(f"Running value for {coeff_name} at Q = {Q}: {run_value}")
    except Exception as e:
        print(f"Error during coefficient retrieval: {e}")

if __name__ == "__main__":
    test_coefficient_manager()