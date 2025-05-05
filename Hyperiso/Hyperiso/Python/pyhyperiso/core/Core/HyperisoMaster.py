from pyhyperiso.core.Core import HyperisoMaster as _CppHyperisoMaster
from pyhyperiso.phyperiso.pyhyperiso.common import Config as _CppConfig
from pyhyperiso.core.Common.GeneralEnum import ExternalFlag, Model


class PyHyperisoMaster:
    """High-level Python wrapper for the C++ HyperisoMaster class.

    Provides a Pythonic interface to initialize and monitor the Hyperiso system.
    """

    def __init__(self):
        """Initializes a new Hyperiso controller."""
        self._cpp_obj = _CppHyperisoMaster()

    def init(self, lha_file: str, config: _CppConfig = None):
        """Initializes Hyperiso with an LHA file and an optional config.

        Args:
            lha_file (str): Path to the LHA input file.
            config (_CppConfig, optional): Native Config object. If not provided,
                a default config will be used.
        """
        if config is not None:
            self._cpp_obj.init(lha_file, config)
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
