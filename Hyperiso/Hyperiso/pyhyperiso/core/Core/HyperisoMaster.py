"""High-level controller for initializing and switching Hyperiso sessions."""

from __future__ import annotations

import os
import tempfile
from enum import Enum
from importlib import resources
from pathlib import Path
from typing import Any, Mapping, Optional, Union

from pyhyperiso.phyperiso.pyhyperiso.core import APIPath as _CppAPIPath
from pyhyperiso.phyperiso.pyhyperiso.core import HyperisoMaster as _CppHyperisoMaster
from pyhyperiso.core.Common.GeneralEnum import Model
from pyhyperiso.core.Core.HyperisoConfig import ExternalFlag, HyperisoConfig

PathLike = Union[str, os.PathLike[str]]


class APIPath(Enum):
    """Filesystem path keys accepted by ``HyperisoMaster.pre_init_set_paths``.

    ``LHA_PATH`` is exposed for read-only diagnostics through ``APIAdapter`` but
    is intentionally rejected by ``pre_init_set_paths`` because the active LHA
    file is provided through ``init(...)`` or ``switch_lha(...)``.
    """

    LHA_PATH = _CppAPIPath.LHA_PATH
    ASSETS_ROOT = _CppAPIPath.ASSETS_ROOT

    DEFAULT_PARAM_VALUES = _CppAPIPath.DEFAULT_PARAM_VALUES
    DEFAULT_OBS_VALUES = _CppAPIPath.DEFAULT_OBS_VALUES
    DEFAULT_PARAM_CORR = _CppAPIPath.DEFAULT_PARAM_CORR
    DEFAULT_OBS_CORR = _CppAPIPath.DEFAULT_OBS_CORR
    DEFAULT_NUISANCES = _CppAPIPath.DEFAULT_NUISANCES

    USER_SM_PARAMS = _CppAPIPath.USER_SM_PARAMS
    USER_FLAVOR_PARAMS = _CppAPIPath.USER_FLAVOR_PARAMS
    USER_DECAY_PARAMS = _CppAPIPath.USER_DECAY_PARAMS
    USER_OBS_VALUES = _CppAPIPath.USER_OBS_VALUES
    USER_PARAM_CORR = _CppAPIPath.USER_PARAM_CORR
    USER_OBS_CORR = _CppAPIPath.USER_OBS_CORR
    USER_NUISANCES = _CppAPIPath.USER_NUISANCES

    PARAM_MAPPING_DIR = _CppAPIPath.PARAM_MAPPING_DIR
    TEMPLATE_DIR = _CppAPIPath.TEMPLATE_DIR
    SPECTRUM_DIR = _CppAPIPath.SPECTRUM_DIR
    MARTY_TEMP_DIR = _CppAPIPath.MARTY_TEMP_DIR


class HyperisoMaster:
    """High-level Python wrapper around the C++ ``HyperisoMaster``.

    The wrapper configures package-friendly defaults before initialization:
    read-only assets are looked up under ``pyhyperiso/assets`` when available,
    and writable caches are placed under the user cache directory unless the
    caller provides explicit directories.
    """

    def __init__(
        self,
        *,
        configure_default_paths: bool = True,
        assets_root: Optional[PathLike] = None,
        cache_root: Optional[PathLike] = None,
    ) -> None:
        """Create an uninitialized Hyperiso controller.

        Args:
            configure_default_paths: When ``True``, configure packaged assets
                and cache directories before any call to ``init``.
            assets_root: Optional replacement for the packaged read-only
                assets directory.
            cache_root: Optional root used to create ``MartyTemp`` and
                ``Spectrum`` writable cache directories.
        """
        self._cpp_obj = _CppHyperisoMaster()
        self.config: Optional[HyperisoConfig] = None

        if configure_default_paths:
            self.configure_default_paths(assets_root=assets_root, cache_root=cache_root)

    @staticmethod
    def _packaged_assets_root() -> Optional[Path]:
        """Return ``pyhyperiso/assets`` when it is available as a filesystem path."""
        try:
            assets = resources.files("pyhyperiso").joinpath("assets")
        except Exception:
            return None

        if assets.is_dir():
            return Path(str(assets)).resolve()
        return None

    @staticmethod
    def _legacy_source_assets_root() -> Optional[Path]:
        """Return the historical source-tree ``Assets`` directory when present."""
        assets = (
            Path(__file__).resolve().parent
            / ".."
            / ".."
            / ".."
            / ".."
            / ".."
            / "Assets"
        ).resolve()
        return assets if assets.is_dir() else None

    @staticmethod
    def _default_cache_root() -> Path:
        """Return a writable cache root without adding a hard platformdirs dependency."""
        env_root = os.environ.get("HYPERISO_CACHE_ROOT")
        if env_root:
            return Path(env_root).expanduser().resolve()

        xdg_cache = os.environ.get("XDG_CACHE_HOME")
        if xdg_cache:
            return (Path(xdg_cache).expanduser() / "pyhyperiso").resolve()

        home = os.environ.get("HOME")
        if home:
            return (Path(home).expanduser() / ".cache" / "pyhyperiso").resolve()

        return (Path(tempfile.gettempdir()) / "pyhyperiso").resolve()

    @staticmethod
    def _to_cpp_api_path(path_key: Any) -> Any:
        """Convert a Python APIPath-like value to the bound C++ enum value."""
        cpp_value = getattr(path_key, "value", path_key)
        if isinstance(cpp_value, str):
            try:
                return getattr(_CppAPIPath, cpp_value)
            except AttributeError as exc:
                raise ValueError(f"Unknown APIPath name: {cpp_value}") from exc
        return cpp_value

    def configure_default_paths(
        self,
        *,
        assets_root: Optional[PathLike] = None,
        cache_root: Optional[PathLike] = None,
    ) -> None:
        """Configure package assets and writable cache directories before init.

        Args:
            assets_root: Optional read-only assets directory. When omitted, the
                wrapper first tries ``pyhyperiso/assets`` and then the historical
                source-tree ``Assets`` directory.
            cache_root: Optional writable root. ``MartyTemp`` and ``Spectrum``
                are created below it. When omitted, a standard user cache root is
                used.
        """
        resolved_assets = (
            Path(os.fspath(assets_root)).expanduser().resolve()
            if assets_root is not None
            else self._packaged_assets_root() or self._legacy_source_assets_root()
        )

        if resolved_assets is not None:
            self.pre_init_set_paths({APIPath.ASSETS_ROOT: resolved_assets})

        resolved_cache_root = (
            Path(os.fspath(cache_root)).expanduser().resolve()
            if cache_root is not None
            else self._default_cache_root()
        )
        self.pre_init_set_marty_cache_dir(resolved_cache_root / "MartyTemp")
        self.pre_init_set_spectrum_cache_dir(resolved_cache_root / "Spectrum")

    @staticmethod
    def _resolve_lha_path(lha_file: PathLike) -> str:
        """Prepare an LHA path for the C++ layer.

        Absolute paths are forwarded as absolute paths. Relative paths are kept
        relative so that the C++ ``MemoryManager`` resolves them under the active
        ``ASSETS_ROOT``. This is important for wheels where ``Assets`` no longer
        lives next to the source checkout.
        """
        path = os.fspath(lha_file)
        return os.path.abspath(path) if os.path.isabs(path) else path

    def init(self, lha_file: PathLike, config: Optional[HyperisoConfig] = None) -> None:
        """Initialize Hyperiso with an LHA file.

        Args:
            lha_file: Absolute LHA path, or path relative to the active
                ``ASSETS_ROOT``.
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

    def pre_init_add_block(
        self,
        block_name: str,
        item_count: int = 2,
        value_idx: int = 1,
        scale_idx: int = -1,
        rg_idx: int = -1,
        bin_idx: int = -1,
        global_scale: bool = False,
    ) -> None:
        """Register an additional LHA block prototype before initialization."""
        self._cpp_obj.pre_init_add_block(
            block_name,
            item_count,
            value_idx,
            scale_idx,
            rg_idx,
            bin_idx,
            global_scale,
        )

    def pre_init_set_marty_path(self, marty_path: PathLike) -> None:
        """Register an existing MARTY installation before initialization."""
        self._cpp_obj.pre_init_set_marty_path(os.path.abspath(os.fspath(marty_path)))

    def pre_init_set_paths(self, path_overrides: Mapping[Any, PathLike]) -> None:
        """Override selected Hyperiso filesystem paths before initialization.

        Args:
            path_overrides: Mapping from ``APIPath`` keys to replacement paths.
                Keys may be this wrapper's ``APIPath`` enum, the low-level C++
                ``APIPath`` enum, another wrapper enum exposing a ``.value``
                C++ enum, or a string matching an ``APIPath`` name. Values may
                be strings or ``os.PathLike`` objects.
        """
        cpp_overrides = {
            self._to_cpp_api_path(path_key): os.path.abspath(os.fspath(path_value))
            for path_key, path_value in path_overrides.items()
        }
        self._cpp_obj.pre_init_set_paths(cpp_overrides)

    def pre_init_set_marty_cache_dir(self, cache_dir: PathLike) -> None:
        """Set the writable MARTY generated-code/cache directory before init."""
        self._cpp_obj.pre_init_set_marty_cache_dir(os.path.abspath(os.fspath(cache_dir)))

    def pre_init_set_spectrum_cache_dir(self, cache_dir: PathLike) -> None:
        """Set the writable spectrum cache directory before init."""
        self._cpp_obj.pre_init_set_spectrum_cache_dir(os.path.abspath(os.fspath(cache_dir)))

    def switch_lha(self, lha_file: PathLike, config: Optional[HyperisoConfig] = None) -> None:
        """Switch the active LHA input file."""
        resolved_lha = self._resolve_lha_path(lha_file)
        if config is not None:
            self._cpp_obj.switch_lha(resolved_lha, config.to_cpp())
            self.config = config
            return

        if self.config is None:
            raise RuntimeError("HyperisoMaster.switch_lha() requires a config before init() has been called.")

        self._cpp_obj.switch_lha(resolved_lha, self.config.to_cpp())

    def check_flag(self, flag: ExternalFlag) -> bool:
        """Return whether an external flag is active."""
        return bool(self._cpp_obj.check_flag(flag.value))

    @property
    def model(self) -> Model:
        """Return the physics model currently active in C++."""
        return Model(self._cpp_obj.get_model())

    def __repr__(self) -> str:
        """Return a compact representation including the active model."""
        return f"<PyHyperisoMaster model={self.model.name}>"


__all__ = ["HyperisoMaster", "APIPath"]
    
    
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
    lha_file_path = "lha/zprime_input.flha"

    print("\n🚀 Calling init with config...")
    hyp.init(lha_file=lha_file_path, config=config)

    print("✅ Current model:", hyp.model.name)
    print("✅ Flag IS_LHA_SPECTRUM:", hyp.check_flag(ExternalFlag.IS_LHA_SPECTRUM))

    hyp.switch_lha("lha/testinput_thdm.lha",config)