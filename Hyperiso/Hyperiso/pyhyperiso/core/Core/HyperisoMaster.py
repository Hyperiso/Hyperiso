"""High-level controller for initializing and switching Hyperiso sessions."""

from __future__ import annotations

import os
from typing import Optional

from pyhyperiso.phyperiso.pyhyperiso.core import HyperisoMaster as _CppHyperisoMaster
from pyhyperiso.core.Common.GeneralEnum import Model
from pyhyperiso.core.Core.HyperisoConfig import ExternalFlag, HyperisoConfig


class HyperisoMaster:
    """High-level Python wrapper around the C++ ``HyperisoMaster``.

    ``HyperisoMaster`` owns the initialization lifecycle of the global C++
    runtime. It loads an LHA file, applies a :class:`HyperisoConfig`, and makes
    the model/flags available to the rest of the Python wrappers.

    Examples:
        >>> hyp = HyperisoMaster()
        >>> hyp.init("lha/si_input.flha", HyperisoConfig())
        >>> hyp.model
        <Model.SM: ...>
    """

    def __init__(self) -> None:
        """Create an uninitialized Hyperiso controller."""
        self._cpp_obj = _CppHyperisoMaster()
        self.config: Optional[HyperisoConfig] = None

    @staticmethod
    def _resolve_lha_path(lha_file: str) -> str:
        """Resolve project-relative LHA paths used by the historical wrapper.

        Args:
            lha_file: Absolute path or path relative to the project ``Assets``
                directory.

        Returns:
            Absolute path passed to the C++ layer.
        """
        if os.path.isabs(lha_file):
            return lha_file
        return os.path.abspath(
            os.path.join(
                os.path.dirname(__file__),
                "..",
                "..",
                "..",
                "..",
                "..",
                "Assets",
                lha_file,
            )
        )

    def init(self, lha_file: str, config: Optional[HyperisoConfig] = None) -> None:
        """Initialize Hyperiso with an LHA file.

        Args:
            lha_file: Absolute LHA path or project-relative path under
                ``Assets``.
            config: Optional initialization config. When omitted, the C++
                default configuration is used and ``self.config`` is set to a
                default :class:`HyperisoConfig` instance.
        """
        resolved_lha = self._resolve_lha_path(lha_file)
        if config is not None:
            self._cpp_obj.init(resolved_lha, config.to_cpp())
            self.config = config
        else:
            self._cpp_obj.init(resolved_lha)
            self.config = HyperisoConfig()

    def switch_lha(self, lha_file: str, config: Optional[HyperisoConfig] = None) -> None:
        """Switch the active LHA input file.

        Args:
            lha_file: Path to the new LHA file. Unlike ``init()``, this method
                forwards the path exactly as provided, matching the historical
                C++ wrapper behavior.
            config: Optional config to apply during the switch. If omitted, the
                previously stored config is reused when available.

        Raises:
            RuntimeError: If no config is available for a config-less switch.
        """
        if config is not None:
            self._cpp_obj.switch_lha(lha_file, config.to_cpp())
            self.config = config
            return

        if self.config is None:
            raise RuntimeError("HyperisoMaster.switch_lha() requires a config before init() has been called.")

        self._cpp_obj.init(lha_file, self.config.to_cpp())

    def check_flag(self, flag: ExternalFlag) -> bool:
        """Return whether an external flag is active.

        Args:
            flag: External flag to query.

        Returns:
            ``True`` when the flag is active in the current C++ runtime.
        """
        return bool(self._cpp_obj.check_flag(flag.value))

    @property
    def model(self) -> Model:
        """Return the physics model currently active in C++."""
        return Model(self._cpp_obj.get_model())

    def __repr__(self) -> str:
        """Return a compact representation including the active model."""
        return f"<PyHyperisoMaster model={self.model.name}>"


__all__ = ["HyperisoMaster"]
    
    
if __name__ == "__main__":
    from pathlib import Path

    print("🔧 Initializing PyHyperisoMaster with custom PyHyperisoConfig...")

    config = HyperisoConfig(
        flags={
            ExternalFlag.IS_LHA_SPECTRUM: True,
            ExternalFlag.HAS_WILSON_INPUT: False,
            ExternalFlag.HAS_TH_OBSERVABLE_INPUT: False,
            ExternalFlag.HYP_AS_SM_MARTY : True
        },
        model=Model.SM,
        mty_model_name="MSSM_UFO",
        mty_model_path=Path("/my/custom/marty/path")
    )

    print("🔧 PyHyperisoConfig content:")
    print(config)

    hyp = HyperisoMaster()
    lha_file_path = "lha/camilia.flha"

    print("\n🚀 Calling init with config...")
    hyp.init(lha_file=lha_file_path, config=config)

    print("✅ Current model:", hyp.model.name)
    print("✅ Flag IS_LHA_SPECTRUM:", hyp.check_flag(ExternalFlag.IS_LHA_SPECTRUM))

    hyp.switch_lha("lha/testinput_thdm.lha",config)