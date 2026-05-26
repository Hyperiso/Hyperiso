"""Python wrappers around C++ enum/name mappers.

The C++ core exposes a family of mapper classes that convert between enum
values, canonical names, and sometimes domain-specific identifiers. This module
keeps those conversions available from Python while returning Python wrapper
objects where appropriate, for example :class:`ObservableId` and
:class:`LhaID`.
"""

from typing import AnyStr

from pyhyperiso.phyperiso.pyhyperiso.common import ContributionTypeMapper as _CppContributionTypeMapper
from pyhyperiso.phyperiso.pyhyperiso.common import DecayMapper as _CppDecayMapper
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
from pyhyperiso.core.Common.SymbolId import ObservableId, _unwrap_optional


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
    """Mapper for Wilson-coefficient groups."""

    def __init__(self):
        pass

    def str(self, group_id: WGroup):
        """Return the canonical string for a Wilson group."""
        return _CppGroupMapper.str(group_id.value)

    def id_of(self, group_id: WGroup):
        """Return the C++ id associated with a Wilson group."""
        return _CppGroupMapper.id_of(self.str(group_id))

    def get_str(self):
        """Return the C++ mapper's primary string table."""
        return _CppGroupMapper.get_str()

    def get_str_all(self):
        """Return all string aliases known to the C++ mapper."""
        return _CppGroupMapper.get_str_all()

    def get_enum(self):
        """Return enum values known to the C++ mapper."""
        return [x for x in _CppGroupMapper.get_enum()]


class WCoefMapper:
    """Mapper for Wilson-coefficient identifiers."""

    def __init__(self):
        pass

    def str(self, obs_id: WCoeff):
        """Return the canonical string for a Wilson coefficient."""
        return _CppWCoefMapper.str(obs_id.value)

    def id_of(self, obs_id: WCoeff):
        """Return the C++ id associated with a Wilson coefficient."""
        return _CppWCoefMapper.id_of(self.str(obs_id))

    def get_str(self):
        """Return the C++ mapper's primary string table."""
        return _CppWCoefMapper.get_str()

    def get_str_all(self):
        """Return all string aliases known to the C++ mapper."""
        return _CppWCoefMapper.get_str_all()

    def get_enum(self):
        """Return enum values known to the C++ mapper."""
        return [x for x in _CppWCoefMapper.get_enum()]


class ObservableMapper:
    """Static helpers converting observable enums, ids, names, and FLHA codes.

    Example:
        >>> oid = ObservableMapper.to_id(Observables.BR_BU_TAU_NU)
        >>> name = ObservableMapper.canonical(oid)
    """

    def __init__(self):
        pass

    @staticmethod
    def str(obs_id: Observables):
        """Return the canonical C++ name for an observable enum."""
        return _CppObservableMapper.str(obs_id.value)

    @staticmethod
    def id_of(name: AnyStr) -> ObservableId:
        """Return the internal observable id associated with a name.

        Args:
            name: Canonical name or alias understood by the C++ mapper.

        Returns:
            ObservableId: Python wrapper around the internal observable id.
        """
        cpp_id = _CppObservableMapper.id_of(name)
        return ObservableId(str(cpp_id))

    @staticmethod
    def to_id(obs: Observables) -> ObservableId:
        """Convert a public observable enum to an internal observable id."""
        cpp_id = _CppObservableMapper.to_id(obs.value)
        return ObservableId(str(cpp_id))

    @staticmethod
    def canonical(obs_id: ObservableId) -> str:
        """Return the canonical name of an internal observable id."""
        return _CppObservableMapper.canonical(obs_id._to_cpp())

    @staticmethod
    def get_str():
        """Return the C++ mapper's primary string table."""
        return _CppObservableMapper.get_str()

    @staticmethod
    def get_str_all():
        """Return all string aliases known to the C++ mapper."""
        return _CppObservableMapper.get_str_all()

    @staticmethod
    def get_enum():
        """Return all public observable enum values known to the mapper."""
        return [Observables(x) for x in _CppObservableMapper.get_enum()]

    @staticmethod
    def from_flha(flha_id: LhaID) -> ObservableId:
        """Convert an FLHA code to an observable id.

        Args:
            flha_id: FLHA identifier wrapper.

        Returns:
            ObservableId: Matching observable id.

        Raises:
            KeyError: If no observable is associated with the FLHA code.
        """
        cpp_opt = _CppObservableMapper.from_flha(flha_id._cpp_obj)
        cpp_id = _unwrap_optional(cpp_opt, "ObservableId from FLHA")
        return ObservableId(str(cpp_id))

    @staticmethod
    def flha(obs) -> LhaID:
        """Return the FLHA code associated with an observable.

        Args:
            obs: Either a public ``Observables`` enum or an internal
                ``ObservableId``.

        Returns:
            LhaID: FLHA code used by the C++ observable database.

        Raises:
            TypeError: If ``obs`` is neither ``Observables`` nor
                ``ObservableId``.
        """
        if isinstance(obs, Observables):
            cpp_lha = _CppObservableMapper.flha(obs.value)
            return LhaID(cpp_lha)
        elif isinstance(obs, ObservableId):
            cpp_lha = _CppObservableMapper.flha(obs._to_cpp())
            return LhaID(cpp_lha)
        else:
            raise TypeError("flha() expects Observables or ObservableId")


class DecayMapper:
    """Mapper linking decay categories and observables."""

    def __init__(self):
        pass

    def str(self, obs_id: Decays):
        """Return the canonical string for a decay enum."""
        return _CppDecayMapper.str(obs_id.value)

    def id_of(self, obs_id: Decays):
        """Return the C++ id associated with a decay enum."""
        return _CppDecayMapper.id_of(self.str(obs_id))

    def get_str(self):
        """Return the C++ mapper's primary string table."""
        return _CppDecayMapper.get_str()

    def get_str_all(self):
        """Return all string aliases known to the C++ mapper."""
        return _CppDecayMapper.get_str_all()

    def get_enum(self):
        """Return enum values known to the C++ mapper."""
        return [x for x in _CppDecayMapper.get_enum()]

    def get_observables(self, decay: Decays):
        """Return public observables attached to a decay category.

        Args:
            decay: Decay enum.

        Returns:
            list[Observables]: Observables associated with the decay.
        """
        cpp_list = _CppDecayMapper.get_observables(decay.value)
        return [Observables(o) for o in cpp_list]

    def get_decay(self, obs: Observables):
        """Return the decay category associated with an observable."""
        cpp_decay = _CppDecayMapper.get_decay(obs.value)
        return Decays(cpp_decay)


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
    "ObservableMapper",
    "DecayMapper",
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
    
