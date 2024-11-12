import phyperiso.core as core
import phyperiso.wilson as wilson

def test_coefficient_manager():
    model_name = "SM"
    # mm = core.MemoryManager.get_instance("Test/InputFiles/testInput.flha", [0])
    # mm.init()
    # param = core.Parameters.get_instance(0)
    # print(param.alpha_s(81))
    manager = wilson.CoefficientManager.get_instance(model_name)
    wilson.CoefficientManager.initialize("Test/InputFiles/testInput.flha", [0])
    print(model_name)
    group = wilson.BCoefficientGroup()
    print("wow1")
    # group_name = "BCoefficientGroup"
    # # group = wilson.BPrimeCoefficientGroup()
    # print("wow")
    # manager.register_coefficient_group(group_name, group)

    # Q_match = 81
    # Q = 42
    # manager.set_q_match(group_name, Q_match)
    # manager.set_group_scale(group_name, Q)
    # order = "LO"
    # coeff_name = "C2"
    
    # try:
    #     matching_value = manager.get_matching_coefficient(group_name, coeff_name, order)
    #     run_value = manager.get_run_coefficient(group_name, coeff_name, order)
    #     print(f"Matching value for {coeff_name} at Q_match = {Q_match}: {matching_value}")
    #     print(f"Running value for {coeff_name} at Q = {Q}: {run_value}")
    # except Exception as e:
    #     print(f"Error during coefficient retrieval: {e}")

if __name__ == "__main__":
    test_coefficient_manager()