from pyhyperiso.phyperiso.pyhyperiso.core import HyperisoMaster as _CppHyperisoMaster
from pyhyperiso.core.Common.GeneralEnum import Model
from pyhyperiso.core.Core.Config import PyConfig, ExternalFlag

class PyHyperisoMaster:
    """High-level Python wrapper for the C++ HyperisoMaster class.

    Provides a Pythonic interface to initialize and monitor the Hyperiso system.
    """

    def __init__(self):
        """Initializes a new Hyperiso controller."""
        self._cpp_obj = _CppHyperisoMaster()

    def init(self, lha_file: str, config: PyConfig = None):
        """Initializes Hyperiso with an LHA file and an optional config.

        Args:
            lha_file (str): Path to the LHA input file.
            config (_CppConfig, optional): Native Config object. If not provided,
                a default config will be used.
        """
        if config is not None:
            self._cpp_obj.init(lha_file, config.to_cpp())
        else:
            self._cpp_obj.init(lha_file)

    def check_flag(self, flag: ExternalFlag) -> bool:
        """Checks whether a specific external flag is active.

        Args:
            flag (ExternalFlag): The Python enum flag to query.

        Returns:
            bool: True if the flag is set, False otherwise.
        """
        return self._cpp_obj.check_flag(flag.value)

    @property
    def model(self) -> Model:
        """Returns the physics model currently used.

        Returns:
            Model: Active model enumeration (e.g., Model.SM, Model.SUSY).
        """
        return Model(self._cpp_obj.get_model())

    def __repr__(self) -> str:
        """String representation of the PyHyperisoMaster instance.

        Returns:
            str: Descriptive string including the model.
        """
        return f"<PyHyperisoMaster model={self.model.name}>"
    
    
if __name__ == "__main__":
    from pathlib import Path

    print("🔧 Initializing PyHyperisoMaster with custom PyConfig...")

    # Création du config avec un flag activé
    config = PyConfig(
        flags={
            ExternalFlag.IS_LHA_SPECTRUM: True,
            ExternalFlag.HAS_WILSON_INPUT: False,
            ExternalFlag.HAS_TH_OBSERVABLE_INPUT: False,
            ExternalFlag.USE_MARTY: False
        },
        model=Model.SM,
        mty_model_name="MSSM_UFO",
        mty_model_path=Path("/my/custom/marty/path")
    )

    print("🔧 PyConfig content:")
    print(config)

    # Initialisation de HyperisoMaster
    hyp = PyHyperisoMaster()
    lha_file_path = "lha/camilia.flha"  # adapte ce chemin à ton repo local

    print("\n🚀 Calling init with config...")
    hyp.init(lha_file=lha_file_path, config=config)

    print("✅ Current model:", hyp.model.name)
    print("✅ Flag IS_LHA_SPECTRUM:", hyp.check_flag(ExternalFlag.IS_LHA_SPECTRUM))
    print("✅ Flag USE_MARTY:", hyp.check_flag(ExternalFlag.USE_MARTY))
