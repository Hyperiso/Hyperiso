"""Python wrappers for decay-configuration objects.

This module mirrors the decay configuration structs exposed by the C++
observable binding. The Python classes are lightweight value wrappers: they
store Python enum values and build the corresponding pybind11 C++ config object
through ``to_cpp()`` when passed to
``ObservableInterface.set_decay_config(...)``.

Example:
    >>> from pyhyperiso.core.BusinessLogic.DecayConfig import BKllConfig, BPFFSource
    >>> from pyhyperiso.core.BusinessLogic.ObservableInterface import ObservableInterface
    >>> from pyhyperiso.core.Common.GeneralEnum import Decays
    >>> oi = ObservableInterface()
    >>> cfg = BKllConfig(ff_src=BPFFSource.GKvD_SR_LAT, n_threads=8)
    >>> oi.set_decay_config(Decays.B__K_l_l, cfg)
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Any, Type, TypeVar

from pyhyperiso.phyperiso.pyhyperiso.observable import (
    B_FF_Type as _CppB_FF_Type,
    BP_FF_Src as _CppBP_FF_Src,
    BV_FF_Src as _CppBV_FF_Src,
    LbL_FF_Src as _CppLbL_FF_Src,
    BDlnuConfig as _CppBDlnuConfig,
    BDstarlnuConfig as _CppBDstarlnuConfig,
    BKllConfig as _CppBKllConfig,
    BKstarGammaConfig as _CppBKstarGammaConfig,
    BKstarllConfig as _CppBKstarllConfig,
    BsPhiConfig as _CppBsPhiConfig,
    BXsllConfig as _CppBXsllConfig,
    DecayConfig as _CppDecayConfig,
    KllDecayConfig as _CppKllDecayConfig,
    LbLllConfig as _CppLbLllConfig,
)


class BFFType(Enum):
    """Form-factor treatment for exclusive B decays.

    Attributes:
        FULL: Use the full form-factor basis.
        SOFT: Use the soft form-factor approximation.
    """

    FULL = _CppB_FF_Type.FULL
    SOFT = _CppB_FF_Type.SOFT


class BPFFSource(Enum):
    """Source of B -> pseudoscalar form factors.

    These values are used by :class:`BKllConfig`.
    """

    AS = _CppBP_FF_Src.AS
    GRvDV = _CppBP_FF_Src.GRvDV
    GKvD_SR_LAT = _CppBP_FF_Src.GKvD_SR_LAT
    GKvD_SR = _CppBP_FF_Src.GKvD_SR
    FLAG24 = _CppBP_FF_Src.FLAG24
    HPQCD22 = _CppBP_FF_Src.HPQCD22


class BVFFSource(Enum):
    """Source of B -> vector form factors.

    These values are used by :class:`BKstarllConfig`,
    :class:`BKstarGammaConfig`, and :class:`BsPhiConfig`.
    """

    BSZ_SR_LAT = _CppBV_FF_Src.BSZ_SR_LAT
    BSZ_SR = _CppBV_FF_Src.BSZ_SR
    GRvDV = _CppBV_FF_Src.GRvDV
    GKvD_SR_LAT = _CppBV_FF_Src.GKvD_SR_LAT
    GKvD_SR = _CppBV_FF_Src.GKvD_SR
    HLMW = _CppBV_FF_Src.HLMW


class LbLFFSource(Enum):
    """Source of Lambda_b -> Lambda form factors.

    Currently only the values exposed by the available C++ binding are listed.
    """

    DM = _CppLbL_FF_Src.DM


class BDlnuBCharge(Enum):
    """B-meson charge choice for :class:`BDlnuConfig`."""

    B_0 = _CppBDlnuConfig.BCharge.B_0
    B_PLUS = _CppBDlnuConfig.BCharge.B_PLUS


class BDstarlnuBCharge(Enum):
    """B-meson charge choice for :class:`BDstarlnuConfig`."""

    B_0 = _CppBDstarlnuConfig.BCharge.B_0
    B_PLUS = _CppBDstarlnuConfig.BCharge.B_PLUS


class BKllBCharge(Enum):
    """B-meson charge choice for :class:`BKllConfig`."""

    B_0 = _CppBKllConfig.BCharge.B_0
    B_PLUS = _CppBKllConfig.BCharge.B_PLUS


class BKllLepton(Enum):
    """Lepton-generation choice for :class:`BKllConfig`."""

    E = _CppBKllConfig.Lepton.E
    MU = _CppBKllConfig.Lepton.MU
    TAU = _CppBKllConfig.Lepton.TAU


class BKstarllPowerCorrectionsImpl(Enum):
    """Power-correction prescription for :class:`BKstarllConfig`."""

    BFS = _CppBKstarllConfig.PowerCorrectionsImpl.BFS
    BCvDV = _CppBKstarllConfig.PowerCorrectionsImpl.BCvDV
    KMPW = _CppBKstarllConfig.PowerCorrectionsImpl.KMPW


class BKstarllBCharge(Enum):
    """B-meson charge choice for :class:`BKstarllConfig`."""

    B_0 = _CppBKstarllConfig.BCharge.B_0
    B_PLUS = _CppBKstarllConfig.BCharge.B_PLUS


class BKstarllLepton(Enum):
    """Lepton-generation choice for :class:`BKstarllConfig`."""

    E = _CppBKstarllConfig.Lepton.E
    MU = _CppBKstarllConfig.Lepton.MU
    TAU = _CppBKstarllConfig.Lepton.TAU


class BKstarGammaBCharge(Enum):
    """B-meson charge choice for :class:`BKstarGammaConfig`."""

    B_0 = _CppBKstarGammaConfig.BCharge.B_0
    B_PLUS = _CppBKstarGammaConfig.BCharge.B_PLUS


class BsPhiLepton(Enum):
    """Lepton-generation choice for :class:`BsPhiConfig`."""

    E = _CppBsPhiConfig.Lepton.E
    MU = _CppBsPhiConfig.Lepton.MU
    TAU = _CppBsPhiConfig.Lepton.TAU


class BXsllLepton(Enum):
    """Lepton-generation choice for :class:`BXsllConfig`."""

    E = _CppBXsllConfig.Lepton.E
    MU = _CppBXsllConfig.Lepton.MU
    TAU = _CppBXsllConfig.Lepton.TAU


class LbLllLepton(Enum):
    """Lepton-generation choice for :class:`LbLllConfig`."""

    E = _CppLbLllConfig.Lepton.E
    MU = _CppLbLllConfig.Lepton.MU
    TAU = _CppLbLllConfig.Lepton.TAU


EnumT = TypeVar("EnumT", bound=Enum)


def _as_cpp_enum(value: EnumT | Any, enum_type: Type[EnumT], name: str) -> Any:
    """Return the C++ enum value stored by a Python wrapper enum.

    Args:
        value: Python enum value or raw pybind11 enum value.
        enum_type: Expected Python enum class.
        name: Field name used in error messages.

    Returns:
        The pybind11 enum value to assign on the C++ config object.

    Raises:
        TypeError: If ``value`` is a different Python enum type.
    """
    if isinstance(value, enum_type):
        return value.value
    if isinstance(value, Enum):
        raise TypeError(f"{name} must be {enum_type.__name__}, got {type(value).__name__}.")
    return value


class DecayConfig:
    """Base decay-configuration wrapper.

    This class is useful for decay engines whose concrete configuration is the
    empty C++ ``DecayConfig`` marker. It is also used as the common Python base
    class for legacy C++ config structs that do not inherit from ``DecayConfig``
    but are still accepted by ``ObservableInterface.set_decay_config`` through
    typed pybind overloads.
    """

    @classmethod
    def from_cpp(cls, cpp_obj: _CppDecayConfig) -> "DecayConfig":
        """Wrap an existing C++ ``DecayConfig`` object.

        Args:
            cpp_obj: Bound C++ config object.

        Returns:
            A base Python configuration wrapper.
        """
        inst = cls()
        inst._cpp_obj = cpp_obj
        return inst

    def to_cpp(self) -> _CppDecayConfig:
        """Build a C++ base decay configuration.

        Returns:
            Bound C++ ``DecayConfig`` instance.
        """
        return _CppDecayConfig()

    def _to_cpp(self) -> _CppDecayConfig:
        """Alias for ``to_cpp()`` used by other wrappers."""
        return self.to_cpp()


@dataclass(slots=True)
class BDlnuConfig(DecayConfig):
    """Configuration for ``B -> D l nu`` decays.

    Args:
        charge: B-meson charge convention. Defaults to ``B_PLUS``.
    """

    charge: BDlnuBCharge = BDlnuBCharge.B_PLUS

    @classmethod
    def from_cpp(cls, cpp_obj: _CppBDlnuConfig) -> "BDlnuConfig":
        """Create a Python wrapper from a C++ ``BDlnuConfig``."""
        return cls(charge=BDlnuBCharge(cpp_obj.charge))

    def to_cpp(self) -> _CppBDlnuConfig:
        """Build the C++ ``BDlnuConfig`` object."""
        cfg = _CppBDlnuConfig()
        cfg.charge = _as_cpp_enum(self.charge, BDlnuBCharge, "charge")
        return cfg


@dataclass(slots=True)
class BDstarlnuConfig(DecayConfig):
    """Configuration for ``B -> D* l nu`` decays.

    Args:
        charge: B-meson charge convention. Defaults to ``B_PLUS``.
    """

    charge: BDstarlnuBCharge = BDstarlnuBCharge.B_PLUS

    @classmethod
    def from_cpp(cls, cpp_obj: _CppBDstarlnuConfig) -> "BDstarlnuConfig":
        """Create a Python wrapper from a C++ ``BDstarlnuConfig``."""
        return cls(charge=BDstarlnuBCharge(cpp_obj.charge))

    def to_cpp(self) -> _CppBDstarlnuConfig:
        """Build the C++ ``BDstarlnuConfig`` object."""
        cfg = _CppBDstarlnuConfig()
        cfg.charge = _as_cpp_enum(self.charge, BDstarlnuBCharge, "charge")
        return cfg


@dataclass(slots=True)
class BKllConfig(DecayConfig):
    """Configuration for exclusive ``B -> K l+ l-`` decays.

    Args:
        ff_src: Source of B -> K form factors.
        ff_type: Form-factor treatment, full or soft.
        charge: B-meson charge convention.
        gen: Lepton generation.
        n_threads: Requested number of worker threads. ``0`` may be interpreted
            by the C++ decay as "use hardware concurrency" when supported.
    """

    ff_src: BPFFSource = BPFFSource.AS
    ff_type: BFFType = BFFType.FULL
    charge: BKllBCharge = BKllBCharge.B_PLUS
    gen: BKllLepton = BKllLepton.MU
    n_threads: int = 1

    @classmethod
    def from_cpp(cls, cpp_obj: _CppBKllConfig) -> "BKllConfig":
        """Create a Python wrapper from a C++ ``BKllConfig``."""
        return cls(
            ff_src=BPFFSource(cpp_obj.ff_src),
            ff_type=BFFType(cpp_obj.ff_type),
            charge=BKllBCharge(cpp_obj.charge),
            gen=BKllLepton(cpp_obj.gen),
            n_threads=int(cpp_obj.n_threads),
        )

    def to_cpp(self) -> _CppBKllConfig:
        """Build the C++ ``BKllConfig`` object."""
        cfg = _CppBKllConfig()
        cfg.ff_src = _as_cpp_enum(self.ff_src, BPFFSource, "ff_src")
        cfg.ff_type = _as_cpp_enum(self.ff_type, BFFType, "ff_type")
        cfg.charge = _as_cpp_enum(self.charge, BKllBCharge, "charge")
        cfg.gen = _as_cpp_enum(self.gen, BKllLepton, "gen")
        cfg.n_threads = int(self.n_threads)
        return cfg


@dataclass(slots=True)
class BKstarllConfig(DecayConfig):
    """Configuration for exclusive ``B -> K* l+ l-`` decays.

    Args:
        ff_src: Source of B -> K* form factors.
        ff_type: Form-factor treatment, full or soft.
        power_corr_impl: Non-factorisable power-correction prescription.
        charge: B-meson charge convention.
        gen: Lepton generation.
        n_threads: Requested number of worker threads. ``0`` may be interpreted
            by the C++ decay as "use hardware concurrency" when supported.
    """

    ff_src: BVFFSource = BVFFSource.BSZ_SR_LAT
    ff_type: BFFType = BFFType.FULL
    power_corr_impl: BKstarllPowerCorrectionsImpl = BKstarllPowerCorrectionsImpl.BFS
    charge: BKstarllBCharge = BKstarllBCharge.B_PLUS
    gen: BKstarllLepton = BKstarllLepton.MU
    n_threads: int = 1

    @classmethod
    def from_cpp(cls, cpp_obj: _CppBKstarllConfig) -> "BKstarllConfig":
        """Create a Python wrapper from a C++ ``BKstarllConfig``."""
        return cls(
            ff_src=BVFFSource(cpp_obj.ff_src),
            ff_type=BFFType(cpp_obj.ff_type),
            power_corr_impl=BKstarllPowerCorrectionsImpl(cpp_obj.power_corr_impl),
            charge=BKstarllBCharge(cpp_obj.charge),
            gen=BKstarllLepton(cpp_obj.gen),
            n_threads=int(cpp_obj.n_threads),
        )

    def to_cpp(self) -> _CppBKstarllConfig:
        """Build the C++ ``BKstarllConfig`` object."""
        cfg = _CppBKstarllConfig()
        cfg.ff_src = _as_cpp_enum(self.ff_src, BVFFSource, "ff_src")
        cfg.ff_type = _as_cpp_enum(self.ff_type, BFFType, "ff_type")
        cfg.power_corr_impl = _as_cpp_enum(self.power_corr_impl, BKstarllPowerCorrectionsImpl, "power_corr_impl")
        cfg.charge = _as_cpp_enum(self.charge, BKstarllBCharge, "charge")
        cfg.gen = _as_cpp_enum(self.gen, BKstarllLepton, "gen")
        cfg.n_threads = int(self.n_threads)
        return cfg


@dataclass(slots=True)
class BKstarGammaConfig(DecayConfig):
    """Configuration for exclusive ``B -> K* gamma`` decays.

    Args:
        ff_src: Source of B -> K* form factors.
        charge: B-meson charge convention.
    """

    ff_src: BVFFSource = BVFFSource.BSZ_SR_LAT
    charge: BKstarGammaBCharge = BKstarGammaBCharge.B_PLUS

    @classmethod
    def from_cpp(cls, cpp_obj: _CppBKstarGammaConfig) -> "BKstarGammaConfig":
        """Create a Python wrapper from a C++ ``BKstarGammaConfig``."""
        return cls(
            ff_src=BVFFSource(cpp_obj.ff_src),
            charge=BKstarGammaBCharge(cpp_obj.charge),
        )

    def to_cpp(self) -> _CppBKstarGammaConfig:
        """Build the C++ ``BKstarGammaConfig`` object."""
        cfg = _CppBKstarGammaConfig()
        cfg.ff_src = _as_cpp_enum(self.ff_src, BVFFSource, "ff_src")
        cfg.charge = _as_cpp_enum(self.charge, BKstarGammaBCharge, "charge")
        return cfg


@dataclass(slots=True)
class BsPhiConfig(DecayConfig):
    """Configuration for exclusive ``Bs -> phi l+ l-`` decays.

    Args:
        ff_src: Source of Bs -> phi form factors.
        ff_type: Form-factor treatment, full or soft.
        gen: Lepton generation.
        n_threads: Requested number of worker threads.
    """

    ff_src: BVFFSource = BVFFSource.BSZ_SR_LAT
    ff_type: BFFType = BFFType.FULL
    gen: BsPhiLepton = BsPhiLepton.MU
    n_threads: int = 1

    @classmethod
    def from_cpp(cls, cpp_obj: _CppBsPhiConfig) -> "BsPhiConfig":
        """Create a Python wrapper from a C++ ``BsPhiConfig``."""
        return cls(
            ff_src=BVFFSource(cpp_obj.ff_src),
            ff_type=BFFType(cpp_obj.ff_type),
            gen=BsPhiLepton(cpp_obj.gen),
            n_threads=int(cpp_obj.n_threads),
        )

    def to_cpp(self) -> _CppBsPhiConfig:
        """Build the C++ ``BsPhiConfig`` object."""
        cfg = _CppBsPhiConfig()
        cfg.ff_src = _as_cpp_enum(self.ff_src, BVFFSource, "ff_src")
        cfg.ff_type = _as_cpp_enum(self.ff_type, BFFType, "ff_type")
        cfg.gen = _as_cpp_enum(self.gen, BsPhiLepton, "gen")
        cfg.n_threads = int(self.n_threads)
        return cfg


@dataclass(slots=True)
class BXsllConfig(DecayConfig):
    """Configuration for inclusive ``B -> X_s l+ l-`` decays.

    Args:
        gen: Lepton generation. Defaults to ``MU``.
    """

    gen: BXsllLepton = BXsllLepton.MU

    @classmethod
    def from_cpp(cls, cpp_obj: _CppBXsllConfig) -> "BXsllConfig":
        """Create a Python wrapper from a C++ ``BXsllConfig``."""
        return cls(gen=BXsllLepton(cpp_obj.gen))

    def to_cpp(self) -> _CppBXsllConfig:
        """Build the C++ ``BXsllConfig`` object."""
        cfg = _CppBXsllConfig()
        cfg.gen = _as_cpp_enum(self.gen, BXsllLepton, "gen")
        return cfg


@dataclass(slots=True)
class KllDecayConfig(DecayConfig):
    """Configuration for ``K_L,S -> l+ l-`` decays.

    Args:
        N_L_sign: Sign convention for the long-distance ``K_L -> gamma gamma``
            term. Defaults to ``1``.
        gen: Lepton generation using the integer convention expected by the C++
            backend. Defaults to ``2``.
    """

    N_L_sign: int = 1
    gen: int = 2

    @classmethod
    def from_cpp(cls, cpp_obj: _CppKllDecayConfig) -> "KllDecayConfig":
        """Create a Python wrapper from a C++ ``KllDecayConfig``."""
        return cls(N_L_sign=int(cpp_obj.N_L_sign), gen=int(cpp_obj.gen))

    def to_cpp(self) -> _CppKllDecayConfig:
        """Build the C++ ``KllDecayConfig`` object."""
        cfg = _CppKllDecayConfig()
        cfg.N_L_sign = int(self.N_L_sign)
        cfg.gen = int(self.gen)
        return cfg


@dataclass(slots=True)
class LbLllConfig(DecayConfig):
    """Configuration for ``Lambda_b -> Lambda l+ l-`` decays.

    Args:
        ff_src: Source of Lambda_b -> Lambda form factors.
        gen: Lepton generation. Defaults to ``MU``.
        n_threads: Requested number of worker threads. ``0`` may be interpreted
            by the C++ layer as ``std::thread::hardware_concurrency``.
    """

    ff_src: LbLFFSource = LbLFFSource.DM
    gen: LbLllLepton = LbLllLepton.MU
    n_threads: int = 1

    @classmethod
    def from_cpp(cls, cpp_obj: _CppLbLllConfig) -> "LbLllConfig":
        """Create a Python wrapper from a C++ ``LbLllConfig``."""
        return cls(
            ff_src=LbLFFSource(cpp_obj.ff_src),
            gen=LbLllLepton(cpp_obj.gen),
            n_threads=int(cpp_obj.n_threads),
        )

    def to_cpp(self) -> _CppLbLllConfig:
        """Build the C++ ``LbLllConfig`` object."""
        cfg = _CppLbLllConfig()
        cfg.ff_src = _as_cpp_enum(self.ff_src, LbLFFSource, "ff_src")
        cfg.gen = _as_cpp_enum(self.gen, LbLllLepton, "gen")
        cfg.n_threads = int(self.n_threads)
        return cfg


# Backward/low-level aliases matching the C++ naming style.
B_FF_Type = BFFType
BP_FF_Src = BPFFSource
BV_FF_Src = BVFFSource
LbL_FF_Src = LbLFFSource


__all__ = [
    "BFFType",
    "BPFFSource",
    "BVFFSource",
    "LbLFFSource",
    "B_FF_Type",
    "BP_FF_Src",
    "BV_FF_Src",
    "LbL_FF_Src",
    "BDlnuBCharge",
    "BDstarlnuBCharge",
    "BKllBCharge",
    "BKllLepton",
    "BKstarllPowerCorrectionsImpl",
    "BKstarllBCharge",
    "BKstarllLepton",
    "BKstarGammaBCharge",
    "BsPhiLepton",
    "BXsllLepton",
    "LbLllLepton",
    "DecayConfig",
    "BDlnuConfig",
    "BDstarlnuConfig",
    "BKllConfig",
    "BKstarllConfig",
    "BKstarGammaConfig",
    "BsPhiConfig",
    "BXsllConfig",
    "KllDecayConfig",
    "LbLllConfig",
]
