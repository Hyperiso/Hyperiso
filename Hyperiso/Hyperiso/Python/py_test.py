from phyperiso.pyhyperiso import core, wilson, observable


def test_core_module():
    print("Testing Core Module...")
    
    # Test MemoryManager
    memory_manager = core.MemoryManager.get_instance()
    memory_manager.init("Test/InputFiles/testInput.flha", core.Model.SM, is_spectrum=True, has_wilsons=True, has_obs=True)
    # input_path = memory_manager.get_input_lha_path()
    # print(f"Input LHA Path: {input_path}")

    # Test Parameters
    parameters = core.Parameters.get_instance(core.ParameterType.SM)
    alpha_s = parameters.alpha_s(91.2)  # Example scale
    print(f"Alpha_s at 91.2 GeV: {alpha_s}")

    # running_mass = parameters.running_mass(4.18, 4.18, 91.2)  # Example for a quark
    # print(f"Running mass: {running_mass}")


def test_wilson_module():
    print("\nTesting Wilson Module...")
    
    # Test Wilson Parameters
    wilson_params = wilson.wilson_parameters.WilsonParameters.get_instance()
    wilson_params.set_mu(1000.0)
    wilson_params.set_mu_w(160.0)
    print("Wilson Parameters set successfully.")

    # Test Coefficient Manager
    coeff_manager = wilson.coefficient_manager.CoefficientManager.get_instance("SM")
    coeff_manager.initialize("Test/InputFiles/testInput.flha", core.Model.SM)
    # state = coeff_manager.get_state()
    # print(f"Coefficient Manager State: {state}")


def test_observable_module():
    print("\nTesting Observable Module...")
    
    # Test Observable Interface
    obs_interface = observable.ObservableInterface()
    obs_value = obs_interface.compute_observable(observable.Observables.ISOSPIN_ASYMMETRY_B_KSTAR_GAMMA)  # Example observable
    print(f"Computed Observable Value: {obs_value}")

    obs_interface.set_param("SMINPUTS", 6, 173.0, core.ParameterType.SM)  # Example of setting a parameter
    print("Parameter set successfully in Observable Interface.")


if __name__ == "__main__":
    try:
        print("a")
        test_core_module()
        test_wilson_module()
        test_observable_module()
        print("\nAll modules tested successfully!")
    except Exception as e:
        print(f"An error occurred during testing: {e}")