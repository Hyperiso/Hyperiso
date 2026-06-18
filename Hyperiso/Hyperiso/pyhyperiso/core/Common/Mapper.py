"""Python wrappers around C++ enum/name mappers.

The C++ core exposes a family of mapper classes that convert between enum
values, canonical names, and sometimes domain-specific identifiers. This module
keeps those conversions available from Python while returning Python wrapper
objects where appropriate, for example :class:`ObservableId` and
:class:`LhaID`.
"""

from typing import AnyStr, Iterable, Optional, Sequence, Union

from pyhyperiso.phyperiso.pyhyperiso.common import ContributionTypeMapper as _CppContributionTypeMapper
from pyhyperiso.phyperiso.pyhyperiso.common import DecayMapper as _CppDecayMapper
from pyhyperiso.phyperiso.pyhyperiso.common import CustomObservableSpec as _CppCustomObservableSpec
from pyhyperiso.phyperiso.pyhyperiso.common import _CppDecayId
from pyhyperiso.phyperiso.pyhyperiso.common import GroupMapper as _CppGroupMapper
from pyhyperiso.phyperiso.pyhyperiso.common import MassTypeMapper as _CppMassTypeMapper
from pyhyperiso.phyperiso.pyhyperiso.common import ModelMapper as _CppModelMapper
from pyhyperiso.phyperiso.pyhyperiso.common import ObservableMapper as _CppObservableMapper
from pyhyperiso.phyperiso.pyhyperiso.common import OrderMapper as _CppOrderMapper
from pyhyperiso.phyperiso.pyhyperiso.common import ParameterTypeMapper as _CppParameterTypeMapper
from pyhyperiso.phyperiso.pyhyperiso.common import ScaleTypeMapper as _CppScaleTypeMapper
from pyhyperiso.phyperiso.pyhyperiso.common import WCoefMapper as _CppWCoefMapper
from pyhyperiso.phyperiso.pyhyperiso.common import WilsonBasisMapper as _CppWilsonBasisMapper

from pyhyperiso.core.Common.GeneralEnum import (
    ContributionType,
    Decays,
    MassType,
    Model,
    Observables,
    ParameterType,
    QCDOrder,
    ScaleType,
    WCoeff,
    WGroup,
    WilsonBasis,
)
from pyhyperiso.core.Common.LhaID import LhaID
try:
    from pyhyperiso.core.Common.SymbolId import ObservableId, DecayId, WGroupId, WCoefId, _unwrap_optional
except ImportError:
    from pyhyperiso.core.Common.SymbolId import ObservableId, WGroupId, WCoefId, _unwrap_optional

    class DecayId:
        """Python fallback wrapper around the bound C++ ``DecayId`` type.

        Add the same class to ``pyhyperiso.core.Common.SymbolId`` if you prefer
        keeping all symbol-id wrappers in one module.
        """

        def __init__(self, value: Union[str, "_CppDecayId"]):
            if isinstance(value, _CppDecayId):
                self._cpp_obj = value
            else:
                self._cpp_obj = _CppDecayId(str(value))

        @property
        def name(self) -> str:
            """Return the canonical string stored by the C++ id."""
            return str(self._cpp_obj)

        def _to_cpp(self):
            """Return the bound C++ object."""
            return self._cpp_obj

        def __str__(self) -> str:
            return self.name

        def __repr__(self) -> str:
            return f"DecayId({self.name!r})"

        def __eq__(self, other) -> bool:
            if isinstance(other, DecayId):
                return self.name.lower() == other.name.lower()
            if isinstance(other, str):
                return self.name.lower() == other.lower()
            return False

        def __hash__(self) -> int:
            return hash(self.name.lower())


def _cpp_lhaid_or_none(ext: Optional[LhaID]):
    """Return the C++ ``LhaID`` object or ``None`` for optional mapper APIs."""
    if ext is None:
        return None
    if isinstance(ext, LhaID):
        return ext._cpp_obj
    return ext


def _wrap_observable_id(cpp_id) -> ObservableId:
    """Wrap a bound C++ ``ObservableId`` as a Python ``ObservableId``."""
    return ObservableId(str(cpp_id))


def _wrap_decay_id(cpp_id) -> DecayId:
    """Wrap a bound C++ ``DecayId`` as a Python ``DecayId``."""
    return DecayId(str(cpp_id))


def _wrap_wgroup_id(cpp_id) -> WGroupId:
    """Wrap a bound C++ ``WGroupId`` as a Python ``WGroupId``."""
    return WGroupId(str(cpp_id))


def _wrap_wcoef_id(cpp_id) -> WCoefId:
    """Wrap a bound C++ ``WCoefId`` as a Python ``WCoefId``."""
    return WCoefId(str(cpp_id))


def _unwrap_optional_or_none(cpp_optional, label: str):
    """Return the value stored in a pybind optional, or ``None`` if empty."""
    if cpp_optional is None:
        return None
    try:
        return _unwrap_optional(cpp_optional, label)
    except KeyError:
        return None


class CustomObservableSpec:
    """Python value object describing a custom observable registration.

    Args:
        canonical: Canonical observable name to register.
        aliases: Optional aliases.
        ext: Optional FLHA id.
    """

    def __init__(
        self,
        canonical: str,
        aliases: Optional[Sequence[str]] = None,
        ext: Optional[LhaID] = None,
    ):
        self.canonical = canonical
        self.aliases = list(aliases or [])
        self.ext = ext

    def _to_cpp(self):
        """Convert this spec to the bound C++ ``CustomObservableSpec``."""
        return _CppCustomObservableSpec(
            self.canonical,
            self.aliases,
            _cpp_lhaid_or_none(self.ext),
        )

    @classmethod
    def from_any(cls, value: Union["CustomObservableSpec", dict, tuple]):
        """Build a spec from a spec object, dict, or tuple.

        Accepted tuple shapes are ``(canonical,)``, ``(canonical, aliases)``,
        and ``(canonical, aliases, ext)``.
        """
        if isinstance(value, cls):
            return value

        if isinstance(value, dict):
            return cls(
                value["canonical"],
                value.get("aliases"),
                value.get("ext"),
            )

        if isinstance(value, tuple):
            if len(value) == 1:
                return cls(value[0])
            if len(value) == 2:
                return cls(value[0], value[1])
            if len(value) == 3:
                return cls(value[0], value[1], value[2])

        raise TypeError(
            "Custom observable spec must be CustomObservableSpec, dict, "
            "or tuple(canonical[, aliases[, ext]])"
        )



class OrderMapper:
    """Mapper for :class:`QCDOrder` values."""

    def __init__(self):
        pass

    def str(self, obs_id: QCDOrder):
        """Return the canonical C++ string for a QCD order.

        Args:
            obs_id: Python ``QCDOrder`` enum value.

        Returns:
            str: Canonical name stored by the C++ mapper.
        """
        return _CppOrderMapper.str(obs_id.value)

    def id_of(self, obs_id: QCDOrder):
        """Round-trip a QCD order through its canonical string.

        Args:
            obs_id: Python ``QCDOrder`` enum value.

        Returns:
            Any: Bound C++ enum/id returned by the mapper.
        """
        return _CppOrderMapper.id_of(self.str(obs_id))

    def get_str(self):
        """Return the C++ mapper's primary string table."""
        return _CppOrderMapper.get_str()

    def get_str_all(self):
        """Return all string aliases known to the C++ mapper."""
        return _CppOrderMapper.get_str_all()

    def get_enum(self):
        """Return enum values known to the C++ mapper."""
        return [x for x in _CppOrderMapper.get_enum()]


class ParameterTypeMapper:
    """Mapper for :class:`ParameterType` namespaces."""

    def __init__(self):
        pass

    def str(self, obs_id: ParameterType):
        """Return the canonical string for a parameter namespace."""
        return _CppParameterTypeMapper.str(obs_id.value)

    def id_of(self, obs_id: ParameterType):
        """Return the C++ id associated with a parameter namespace."""
        return _CppParameterTypeMapper.id_of(self.str(obs_id))

    def get_str(self):
        """Return the C++ mapper's primary string table."""
        return _CppParameterTypeMapper.get_str()

    def get_str_all(self):
        """Return all string aliases known to the C++ mapper."""
        return _CppParameterTypeMapper.get_str_all()

    def get_enum(self):
        """Return enum values known to the C++ mapper."""
        return [x for x in _CppParameterTypeMapper.get_enum()]


class ModelMapper:
    """Mapper for physics :class:`Model` values."""

    def __init__(self):
        pass

    def str(self, obs_id: Model):
        """Return the canonical model name."""
        return _CppModelMapper.str(obs_id.value)

    def id_of(self, obs_id: Model):
        """Return the C++ id associated with a model."""
        return _CppModelMapper.id_of(self.str(obs_id))

    def get_str(self):
        """Return the C++ mapper's primary string table."""
        return _CppModelMapper.get_str()

    def get_str_all(self):
        """Return all string aliases known to the C++ mapper."""
        return _CppModelMapper.get_str_all()

    def get_enum(self):
        """Return enum values known to the C++ mapper."""
        return [x for x in _CppModelMapper.get_enum()]


class WilsonBasisMapper:
    """Mapper for Wilson-coefficient basis conventions."""

    def __init__(self):
        pass

    def str(self, obs_id: WilsonBasis):
        """Return the canonical string for a Wilson basis."""
        return _CppWilsonBasisMapper.str(obs_id.value)

    def id_of(self, obs_id: WilsonBasis):
        """Return the C++ id associated with a Wilson basis."""
        return _CppWilsonBasisMapper.id_of(self.str(obs_id))

    def get_str(self):
        """Return the C++ mapper's primary string table."""
        return _CppWilsonBasisMapper.get_str()

    def get_str_all(self):
        """Return all string aliases known to the C++ mapper."""
        return _CppWilsonBasisMapper.get_str_all()

    def get_enum(self):
        """Return enum values known to the C++ mapper."""
        return [x for x in _CppWilsonBasisMapper.get_enum()]


class ContributionTypeMapper:
    """Mapper for Wilson-coefficient contribution components.

    Contribution types typically distinguish SM-only, BSM-only, and total
    coefficients in the Wilson pipeline.
    """

    def __init__(self):
        pass

    def str(self, obs_id: ContributionType):
        """Return the canonical contribution-type string."""
        return _CppContributionTypeMapper.str(obs_id.value)

    def id_of(self, obs_id: ContributionType):
        """Return the C++ id associated with a contribution type."""
        return _CppContributionTypeMapper.id_of(self.str(obs_id))

    def get_str(self):
        """Return the C++ mapper's primary string table."""
        return _CppContributionTypeMapper.get_str()

    def get_str_all(self):
        """Return all string aliases known to the C++ mapper."""
        return _CppContributionTypeMapper.get_str_all()

    def get_enum(self):
        """Return enum values known to the C++ mapper."""
        return [x for x in _CppContributionTypeMapper.get_enum()]


class MassTypeMapper:
    """Mapper for QCD mass conventions."""

    def __init__(self):
        pass

    def str(self, obs_id: MassType):
        """Return the canonical mass-type string."""
        return _CppMassTypeMapper.str(obs_id.value)

    def id_of(self, obs_id: MassType):
        """Return the C++ id associated with a mass convention."""
        return _CppMassTypeMapper.id_of(self.str(obs_id))

    def get_str(self):
        """Return the C++ mapper's primary string table."""
        return _CppMassTypeMapper.get_str()

    def get_str_all(self):
        """Return all string aliases known to the C++ mapper."""
        return _CppMassTypeMapper.get_str_all()

    def get_enum(self):
        """Return enum values known to the C++ mapper."""
        return [x for x in _CppMassTypeMapper.get_enum()]


class ScaleTypeMapper:
    """Mapper for scale categories, such as matching and hadronic scales."""

    def __init__(self):
        pass

    def str(self, obs_id: ScaleType):
        """Return the canonical scale-type string."""
        return _CppScaleTypeMapper.str(obs_id.value)

    def id_of(self, obs_id: ScaleType):
        """Return the C++ id associated with a scale category."""
        return _CppScaleTypeMapper.id_of(self.str(obs_id))

    def get_str(self):
        """Return the C++ mapper's primary string table."""
        return _CppScaleTypeMapper.get_str()

    def get_str_all(self):
        """Return all string aliases known to the C++ mapper."""
        return _CppScaleTypeMapper.get_str_all()

    def get_enum(self):
        """Return enum values known to the C++ mapper."""
        return [x for x in _CppScaleTypeMapper.get_enum()]


class GroupMapper:
    """Mapper for Wilson-coefficient groups.

    The wrapper accepts both legacy :class:`WGroup` enums and runtime strings.
    String lookups return :class:`WGroupId`, which is the preferred API for
    custom groups and config/CLI driven workflows.
    """

    def __init__(self):
        pass

    @staticmethod
    def str(group) -> str:
        """Return the canonical string for a Wilson group enum or dynamic id."""
        if isinstance(group, WGroup):
            return _CppGroupMapper.str(group.value)
        if isinstance(group, WGroupId):
            return _CppGroupMapper.canonical(group._to_cpp())
        raise TypeError("str() expects WGroup or WGroupId")

    @staticmethod
    def id_of(group) -> WGroupId:
        """Resolve a group enum/name/alias to a dynamic :class:`WGroupId`."""
        if isinstance(group, WGroupId):
            return group
        if isinstance(group, WGroup):
            return _wrap_wgroup_id(_CppGroupMapper.to_id(group.value))
        return _wrap_wgroup_id(_CppGroupMapper.id_of(str(group)))

    @staticmethod
    def to_id(group: WGroup) -> WGroupId:
        """Convert a builtin Wilson-group enum to :class:`WGroupId`."""
        return _wrap_wgroup_id(_CppGroupMapper.to_id(group.value))

    @staticmethod
    def canonical(group_id: WGroupId) -> str:
        """Return the canonical name of a dynamic Wilson group id."""
        return _CppGroupMapper.canonical(group_id._to_cpp())

    @staticmethod
    def block_name(group, scale: ScaleType, basis: WilsonBasis = WilsonBasis.STANDARD) -> str:
        """Return the scale/basis block name used by Wilson parameter blocks."""
        if isinstance(group, WGroup):
            return _CppGroupMapper.str(group.value, scale.value, basis.value)
        group_id = GroupMapper.id_of(group)
        return _CppGroupMapper.str_id(group_id._to_cpp(), scale.value, basis.value)

    @staticmethod
    def register_custom(canonical: str, aliases=None, external: Optional[str] = None) -> bool:
        """Register a custom Wilson group in the dynamic mapper.

        Args:
            canonical: Canonical group name.
            aliases: Optional aliases accepted by :meth:`id_of`.
            external: Optional external key used by C++ mapper extensions.
        """
        return _CppGroupMapper.register_custom(canonical, list(aliases or []), external)

    def get_str(self):
        """Return the C++ mapper's primary string table."""
        return _CppGroupMapper.get_str()

    def get_str_all(self):
        """Return all string aliases known to the C++ mapper."""
        return _CppGroupMapper.get_str_all()

    def get_enum(self):
        """Return builtin enum values known to the C++ mapper."""
        return [WGroup(x) for x in _CppGroupMapper.get_enum()]


class WCoefMapper:
    """Mapper for Wilson-coefficient identifiers.

    Runtime lookups return :class:`WCoefId`, so custom coefficients do not need
    a corresponding static :class:`WCoeff` enum value.
    """

    def __init__(self):
        pass

    @staticmethod
    def str(coef) -> str:
        """Return the canonical string for a Wilson coefficient."""
        if isinstance(coef, WCoeff):
            return _CppWCoefMapper.str(coef.value)
        if isinstance(coef, WCoefId):
            return _CppWCoefMapper.canonical(coef._to_cpp())
        raise TypeError("str() expects WCoeff or WCoefId")

    @staticmethod
    def id_of(coef) -> WCoefId:
        """Resolve a coefficient enum/name/alias to :class:`WCoefId`."""
        if isinstance(coef, WCoefId):
            return coef
        if isinstance(coef, WCoeff):
            return _wrap_wcoef_id(_CppWCoefMapper.to_id(coef.value))
        return _wrap_wcoef_id(_CppWCoefMapper.id_of(str(coef)))

    @staticmethod
    def to_id(coef: WCoeff) -> WCoefId:
        """Convert a builtin coefficient enum to :class:`WCoefId`."""
        return _wrap_wcoef_id(_CppWCoefMapper.to_id(coef.value))

    @staticmethod
    def canonical(coef_id: WCoefId) -> str:
        """Return the canonical name of a dynamic coefficient id."""
        return _CppWCoefMapper.canonical(coef_id._to_cpp())

    @staticmethod
    def register_custom(canonical: str, aliases=None, flha=(0, 0)) -> bool:
        """Register a custom Wilson coefficient.

        Args:
            canonical: Canonical coefficient name.
            aliases: Optional aliases accepted by :meth:`id_of`.
            flha: External FLHA pair used by Wilson input/output conventions.
        """
        return _CppWCoefMapper.register_custom(canonical, list(aliases or []), tuple(flha))

    @staticmethod
    def flha_base(coef) -> tuple[int, int]:
        """Return the external FLHA base pair for a coefficient enum or id."""
        if isinstance(coef, WCoeff):
            return tuple(_CppWCoefMapper.flha_base(coef.value))
        coef_id = WCoefMapper.id_of(coef)
        return tuple(_CppWCoefMapper.flha_base(coef_id._to_cpp()))

    def get_str(self):
        """Return the C++ mapper's primary string table."""
        return _CppWCoefMapper.get_str()

    def get_str_all(self):
        """Return all string aliases known to the C++ mapper."""
        return _CppWCoefMapper.get_str_all()

    def get_enum(self):
        """Return builtin enum values known to the C++ mapper."""
        return [WCoeff(x) for x in _CppWCoefMapper.get_enum()]


class ObservableMapper:
    """Mapper for builtin and custom observables.

    Dynamic/user-facing code should use :class:`ObservableId`; the static
    :class:`Observables` enum is kept for builtin legacy paths only.
    """

    def __init__(self):
        pass

    @staticmethod
    def str(obs: Observables) -> str:
        """Return the canonical C++ name for a builtin observable enum."""
        return _CppObservableMapper.str(obs.value)

    @staticmethod
    def enum_elt_legacy(name: AnyStr) -> Observables:
        """Resolve a builtin observable name to the static enum.

        Custom observables cannot be represented by :class:`Observables`; use
        :meth:`id_of` for dynamic lookup.
        """
        return Observables(_CppObservableMapper.enum_elt_legacy(name))

    @staticmethod
    def id_of(name: AnyStr) -> ObservableId:
        """Resolve a canonical name or alias to an :class:`ObservableId`."""
        return _wrap_observable_id(_CppObservableMapper.id_of(name))

    @staticmethod
    def to_id(obs: Observables) -> ObservableId:
        """Convert a builtin observable enum to an :class:`ObservableId`."""
        return _wrap_observable_id(_CppObservableMapper.to_id(obs.value))

    @staticmethod
    def canonical(obs_id: ObservableId) -> str:
        """Return the canonical name of an internal observable id."""
        return _CppObservableMapper.canonical(obs_id._to_cpp())

    @staticmethod
    def get_str():
        """Return builtin observable names."""
        return _CppObservableMapper.get_str()

    @staticmethod
    def get_str_all():
        """Return builtin and custom observable names."""
        return _CppObservableMapper.get_str_all()

    @staticmethod
    def get_enum():
        """Return builtin observable enum values known to the mapper."""
        return [Observables(x) for x in _CppObservableMapper.get_enum()]

    @staticmethod
    def from_flha(flha_id: LhaID) -> ObservableId:
        """Convert an FLHA code to an observable id.

        Raises:
            KeyError: If no observable is associated with the FLHA code.
        """
        cpp_opt = _CppObservableMapper.from_flha(flha_id._cpp_obj)
        cpp_id = _unwrap_optional(cpp_opt, "ObservableId from FLHA")
        return _wrap_observable_id(cpp_id)

    @staticmethod
    def flha(obs) -> LhaID:
        """Return the FLHA code associated with an observable.

        Args:
            obs: Either a builtin :class:`Observables` or an
                :class:`ObservableId`.
        """
        if isinstance(obs, Observables):
            return LhaID(_CppObservableMapper.flha(obs.value))
        if isinstance(obs, ObservableId):
            return LhaID(_CppObservableMapper.flha(obs._to_cpp()))
        raise TypeError("flha() expects Observables or ObservableId")

    @staticmethod
    def register_custom(
        canonical: str,
        parent_decay: Union[str, Decays, DecayId],
        aliases: Optional[Sequence[str]] = None,
        ext: Optional[LhaID] = None,
    ) -> bool:
        """Register a custom observable and attach it to a parent decay.

        Args:
            canonical: Canonical observable name.
            parent_decay: Parent decay as name/alias, :class:`Decays`, or
                :class:`DecayId`.
            aliases: Optional observable aliases.
            ext: Optional FLHA id.

        Returns:
            bool: ``True`` if the C++ registry accepted the symbol.
        """
        aliases = list(aliases or [])
        cpp_ext = _cpp_lhaid_or_none(ext)

        if isinstance(parent_decay, Decays):
            return _CppObservableMapper.register_custom_with_decay_enum(
                canonical,
                parent_decay.value,
                aliases,
                cpp_ext,
            )

        if isinstance(parent_decay, DecayId):
            return _CppObservableMapper.register_custom_with_decay_id(
                canonical,
                parent_decay._to_cpp(),
                aliases,
                cpp_ext,
            )

        return _CppObservableMapper.register_custom(
            canonical,
            str(parent_decay),
            aliases,
            cpp_ext,
        )


class DecayMapper:
    """Mapper for builtin and custom decays."""

    def __init__(self):
        pass

    @staticmethod
    def str(decay: Union[Decays, DecayId]) -> str:
        """Return the canonical string for a builtin enum or dynamic id."""
        if isinstance(decay, Decays):
            return _CppDecayMapper.str(decay.value)
        if isinstance(decay, DecayId):
            return _CppDecayMapper.canonical(decay._to_cpp())
        raise TypeError("str() expects Decays or DecayId")

    @staticmethod
    def enum_elt_legacy(name: AnyStr) -> Decays:
        """Resolve a builtin decay name to the static enum."""
        return Decays(_CppDecayMapper.enum_elt_legacy(name))

    @staticmethod
    def id_of(decay: Union[AnyStr, Decays]) -> DecayId:
        """Resolve a decay name/alias or enum to :class:`DecayId`."""
        if isinstance(decay, Decays):
            return DecayMapper.to_id(decay)
        return _wrap_decay_id(_CppDecayMapper.id_of(decay))

    @staticmethod
    def to_id(decay: Decays) -> DecayId:
        """Convert a builtin decay enum to a dynamic :class:`DecayId`."""
        return _wrap_decay_id(_CppDecayMapper.to_id(decay.value))

    @staticmethod
    def canonical(decay_id: DecayId) -> str:
        """Return the canonical name of a dynamic decay id."""
        return _CppDecayMapper.canonical(decay_id._to_cpp())

    @staticmethod
    def get_str():
        """Return builtin decay names."""
        return _CppDecayMapper.get_str()

    @staticmethod
    def get_str_all():
        """Return builtin and custom decay names."""
        return _CppDecayMapper.get_str_all()

    @staticmethod
    def get_enum():
        """Return builtin decay enum values."""
        return [Decays(x) for x in _CppDecayMapper.get_enum()]

    @staticmethod
    def get_observables(decay: Union[Decays, DecayId, str]):
        """Return observables attached to a decay.

        Returns:
            list[Observables] for builtin :class:`Decays` input, and
            list[ObservableId] for :class:`DecayId` or string input.
        """
        if isinstance(decay, Decays):
            return [Observables(o) for o in _CppDecayMapper.get_observables(decay.value)]

        if isinstance(decay, DecayId):
            return [
                _wrap_observable_id(o)
                for o in _CppDecayMapper.get_observables(decay._to_cpp())
            ]

        return [
            _wrap_observable_id(o)
            for o in _CppDecayMapper.get_observables_by_name(str(decay))
        ]

    @staticmethod
    def get_observable_ids(decay: Decays):
        """Return builtin decay observables converted to :class:`ObservableId`."""
        return [
            _wrap_observable_id(o)
            for o in _CppDecayMapper.get_observable_ids(decay.value)
        ]

    @staticmethod
    def get_decay(obs: Observables) -> Decays:
        """Return the builtin decay enum associated with a builtin observable."""
        return Decays(_CppDecayMapper.get_decay(obs.value))

    @staticmethod
    def get_decay_id(obs: Union[Observables, ObservableId]) -> Optional[DecayId]:
        """Return the dynamic parent decay id of an observable, if known."""
        if isinstance(obs, Observables):
            cpp_opt = _CppDecayMapper.get_decay_id(obs.value)
        elif isinstance(obs, ObservableId):
            cpp_opt = _CppDecayMapper.get_decay_id(obs._to_cpp())
        else:
            raise TypeError("get_decay_id() expects Observables or ObservableId")

        cpp_id = _unwrap_optional_or_none(cpp_opt, "DecayId")
        if cpp_id is None:
            return None
        return _wrap_decay_id(cpp_id)

    @staticmethod
    def get_decay_id_or_throw(obs: ObservableId) -> DecayId:
        """Return the dynamic parent decay id or raise from C++ if absent."""
        return _wrap_decay_id(_CppDecayMapper.get_decay_id_or_throw(obs._to_cpp()))

    @staticmethod
    def has_observables(decay: Union[Decays, DecayId, str]) -> bool:
        """Return whether a decay has at least one observable."""
        if isinstance(decay, Decays):
            decay = DecayMapper.to_id(decay)
        elif isinstance(decay, str):
            decay = DecayMapper.id_of(decay)
        if not isinstance(decay, DecayId):
            raise TypeError("has_observables() expects Decays, DecayId, or str")
        return _CppDecayMapper.has_observables(decay._to_cpp())

    @staticmethod
    def register_custom_with_observables(
        canonical: str,
        observables: Sequence[Union[CustomObservableSpec, dict, tuple]],
        aliases: Optional[Sequence[str]] = None,
    ) -> bool:
        """Register a custom decay together with at least one observable."""
        cpp_specs = [CustomObservableSpec.from_any(o)._to_cpp() for o in observables]
        return _CppDecayMapper.register_custom_with_observables(
            canonical,
            cpp_specs,
            list(aliases or []),
        )

__all__ = [
    "OrderMapper",
    "ParameterTypeMapper",
    "ModelMapper",
    "WilsonBasisMapper",
    "ContributionTypeMapper",
    "MassTypeMapper",
    "ScaleTypeMapper",
    "GroupMapper",
    "WCoefMapper",
    "WGroupId",
    "WCoefId",
    "ObservableMapper",
    "DecayMapper",
    "CustomObservableSpec",
    "DecayId",
]

if __name__ == "__main__":
    
    ordermapper = OrderMapper()
    
    print(ordermapper.get_str())
    print(ordermapper.get_str_all())
    print(ordermapper.get_enum())
    
    gm = GroupMapper()
    
    print(gm.str(WGroup.B))
    
    print(gm.id_of(WGroup.B))
    print(type(gm.id_of(WGroup.B)))
    
    om = ObservableMapper()
    
    print(om.str(Observables.A_FB_B__D_TAU_NU))
    
    print(om.flha(Observables.BR_BU_TAU_NU))
    print(type(om.flha(Observables.BR_BU_TAU_NU)))
    
    decaymapper = DecayMapper()
    

    print(decaymapper.get_observables(Decays.B__D_l_nu))
    print(decaymapper.get_decay(Observables.BR_D__MU_NU))
    
    wcoefmapper = WCoefMapper()
    
