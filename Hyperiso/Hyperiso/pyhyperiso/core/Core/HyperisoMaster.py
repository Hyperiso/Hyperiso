from pyhyperiso.phyperiso.pyhyperiso.core import HyperisoMaster as _CppHyperisoMaster
from pyhyperiso.core.Common.GeneralEnum import Model
from Hyperiso.Hyperiso.pyhyperiso.core.Core.HyperisoConfig import PyHyperisoConfig, ExternalFlag
import os

class PyHyperisoMaster:
    """High-level Python wrapper for the C++ HyperisoMaster class.

    Provides a Pythonic interface to initialize and monitor the Hyperiso system.
    """

    def __init__(self):
        """Initializes a new Hyperiso controller."""
        self._cpp_obj = _CppHyperisoMaster()
        self._config = None
        
    def init(self, lha_file: str, config: PyHyperisoConfig = None):
        """Initializes Hyperiso with an LHA file and an optional config.

        Args:
            lha_file (str): Path to the LHA input file.
            config (PyHyperisoConfig, optional): Config object with Hyperiso input flags. If not provided,
                a default config will be used.
        """
        print("hein ?")
        if not os.path.isabs(lha_file):
            lha_file = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", "..", "..", "Assets", lha_file))
            print(lha_file)
        if config is not None:
            self._cpp_obj.init(lha_file, config.to_cpp())
            self.config = config
        else:
            self._cpp_obj.init(lha_file)
            self.config = PyHyperisoConfig()

    def switch_lha(self, lha_file: str, config: PyHyperisoConfig = None):
        """Initializes Hyperiso with an LHA file and an optional config.

        Args:
            lha_file (str): Path to the LHA input file.
            config (PyHyperisoConfig, optional): Basic Config with hyperiso flags inputs. If not provided,
                a default config will be used.
        """
        if config is not None:
            self._cpp_obj.switch_lha(lha_file, config.to_cpp())
            self.config = config
        else:
            self._cpp_obj.init(lha_file, self.config.to_cpp())
            
            
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

    print("🔧 Initializing PyHyperisoMaster with custom PyHyperisoConfig...")

    config = PyHyperisoConfig(
        flags={
            ExternalFlag.IS_LHA_SPECTRUM: True,
            ExternalFlag.HAS_WILSON_INPUT: False,
            ExternalFlag.HAS_TH_OBSERVABLE_INPUT: False,
            # ExternalFlag.USE_MARTY: False
        },
        model=Model.SM,
        mty_model_name="MSSM_UFO",
        mty_model_path=Path("/my/custom/marty/path")
    )

    print("🔧 PyHyperisoConfig content:")
    print(config)

    hyp = PyHyperisoMaster()
    lha_file_path = "lha/camilia.flha"

    print("\n🚀 Calling init with config...")
    hyp.init(lha_file=lha_file_path, config=config)

    print("✅ Current model:", hyp.model.name)
    print("✅ Flag IS_LHA_SPECTRUM:", hyp.check_flag(ExternalFlag.IS_LHA_SPECTRUM))

    hyp.switch_lha("lha/testinput_thdm.lha",config)