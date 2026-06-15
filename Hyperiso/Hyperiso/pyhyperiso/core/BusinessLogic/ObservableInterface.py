"""Python wrapper for the observable computation interface.

This module exposes a Pythonic facade over the C++ ``ObservableInterface``.
It lets users select observables, compute theory predictions, inspect
experimental inputs, manage observable dependencies, and update parameters
used by the observable layer.

The public API accepts and returns Python wrapper objects only. Conversion to
and from the pybind11 objects is intentionally kept inside this module.

Example:
    >>> from pyhyperiso.core.BusinessLogic.ObservableInterface import ObservableInterface
    >>> from pyhyperiso.core.Common.GeneralEnum import Observables, QCDOrder
    >>> oi = ObservableInterface()
    >>> oi.add_observable(Observables.BR_BS_MUMU, QCDOrder.NNLO, add_dependencies=True)
    >>> values = oi.compute_observable(Observables.BR_BS_MUMU)
    >>> central = oi.compute_observable_central(Observables.BR_BS_MUMU)
"""

from __future__ import annotations

from typing import Dict, List, Mapping, Set

from pyhyperiso.phyperiso.pyhyperiso.observable import ObservableInterface as _CppObservableInterface
from pyhyperiso.core.BusinessLogic.ObservableValue import ObservableValue
from pyhyperiso.core.Common.BinnedObservableId import BinnedObservableId
from pyhyperiso.core.Common.GeneralEnum import Decays, Observables, ParameterType, QCDOrder, UncertaintyType
from pyhyperiso.core.Common.LhaID import LhaID
from pyhyperiso.core.Common.ParamId import ParamId
from pyhyperiso.core.Common.SymbolId import ObservableId
from pyhyperiso.core.BusinessLogic.DecayConfig import DecayConfig

def _require(value, typ, name: str):
    if not isinstance(value, typ):
        raise TypeError(f"{name} doit être {typ.__name__}, reçu {type(value)!r}.")
    return value


def _cpp_observable_enum(obs: Observables):
    return _require(obs, Observables, "obs").value


def _cpp_observable_id(obs: ObservableId):
    return _require(obs, ObservableId, "obs")._to_cpp()


def _cpp_binned_observable_id(obs: BinnedObservableId):
    return _require(obs, BinnedObservableId, "obs").to_cpp()


def _cpp_qcd_order(order: QCDOrder):
    return _require(order, QCDOrder, "qcd_order").value


def _cpp_decay(decay: Decays):
    return _require(decay, Decays, "decay").value


def _cpp_uncertainty_type(u_type: UncertaintyType):
    return _require(u_type, UncertaintyType, "u_type").value


def _cpp_param_id(pid: ParamId):
    return _require(pid, ParamId, "pid").to_cpp()


def _param_from_cpp(cpp_obj) -> ParamId:
    return ParamId.from_cpp(cpp_obj)


def _single_lha_code(code: LhaID) -> int:
    _require(code, LhaID, "pid.code")
    parts = code.get_parts()
    if len(parts) != 1:
        raise ValueError("set_param/get_param via ObservableInterface exige un LhaID à une seule entrée.")
    return int(parts[0])


class ObservableInterface:
    """High-level API used to configure and compute observables.

    ``ObservableInterface`` mirrors the C++ class of the same name. It owns an
    underlying C++ observable manager and exposes a chainable Python API for
    common workflows:

    * register observables by enum, by internal ``ObservableId`` or by binned id;
    * compute one observable or all registered observables;
    * inspect experimental central values and uncertainties;
    * add explicit parameter dependencies for likelihood/statistical workflows;
    * set or read model parameters through the global parameter store.

    Notes:
        ``add_dependencies=True`` asks the C++ observable manager to attach the
        full dependency allow-list known for the selected observable. This is
        important for subsequent statistical workflows, because those
        dependencies are used to decide which nuisance parameters should be
        propagated or profiled.

        The QCD order is passed to the C++ decay computation. Depending on the
        backend and observable, the actual Wilson-coefficient order may be
        limited by backend capabilities.

    Example:
        >>> oi = ObservableInterface()
        >>> oi.add_observable(Observables.BR_BS_MUMU, QCDOrder.NNLO, True)
        >>> oi.get_current_observables()
        [...]
        >>> oi.compute_observable_central(Observables.BR_BS_MUMU)
        3.6e-09
    """

    def __init__(self) -> None:
        """Create a new observable interface backed by a fresh C++ manager."""
        self._cpp_obj = _CppObservableInterface()

    @classmethod
    def from_cpp(cls, cpp_obj) -> "ObservableInterface":
        """Wrap an existing C++ ``ObservableInterface`` instance.

        Args:
            cpp_obj: Bound C++ object produced by pybind11.

        Returns:
            A Python ``ObservableInterface`` sharing ownership of ``cpp_obj``.
        """
        inst = cls.__new__(cls)
        inst._cpp_obj = cpp_obj
        return inst

    def _to_cpp(self):
        """Return the underlying pybind11 object.

        Returns:
            The wrapped C++ ``ObservableInterface`` object.
        """
        return self._cpp_obj

    def add_observable(
        self,
        obs: Observables,
        qcd_order: QCDOrder,
        add_dependencies: bool = True,
    ) -> "ObservableInterface":
        """Register an observable using the public observable enum.

        Args:
            obs: Observable enum value to register.
            qcd_order: Maximum QCD order requested for the corresponding decay.
            add_dependencies: Whether to also attach the C++ dependency
                allow-list for the observable.

        Returns:
            ``self`` to support fluent chaining.
        """
        self._cpp_obj.add_observable(_cpp_observable_enum(obs), _cpp_qcd_order(qcd_order), bool(add_dependencies))
        return self

    def add_observable_id(
        self,
        obs: ObservableId,
        qcd_order: QCDOrder,
        add_dependencies: bool = True,
    ) -> "ObservableInterface":
        """Register an observable using an internal ``ObservableId``.

        Args:
            obs: Internal observable identifier.
            qcd_order: Maximum QCD order requested for the corresponding decay.
            add_dependencies: Whether to also attach known parameter
                dependencies for this observable.

        Returns:
            ``self`` to support fluent chaining.
        """
        self._cpp_obj.add_observable(_cpp_observable_id(obs), _cpp_qcd_order(qcd_order), bool(add_dependencies))
        return self

    def add_binned_observable(
        self,
        obs: BinnedObservableId,
        qcd_order: QCDOrder,
        add_dependencies: bool = True,
    ) -> "ObservableInterface":
        """Register a single binned observable.

        Args:
            obs: Binned observable identifier, including the observable id and
                bin range.
            qcd_order: Maximum QCD order requested for the corresponding decay.
            add_dependencies: Whether to also attach known parameter
                dependencies for this observable.

        Returns:
            ``self`` to support fluent chaining.
        """
        self._cpp_obj.add_observable(_cpp_binned_observable_id(obs), _cpp_qcd_order(qcd_order), bool(add_dependencies))
        return self

    def add_observables(
        self,
        obs_names: Mapping[Observables, QCDOrder],
        add_dependencies: bool = True,
    ) -> "ObservableInterface":
        """Register several enum-based observables at once.

        Args:
            obs_names: Mapping from observable enum values to requested QCD
                orders.
            add_dependencies: Whether to attach dependencies for every
                observable in the mapping.

        Returns:
            ``self`` to support fluent chaining.

        Example:
            >>> oi.add_observables({
            ...     Observables.BR_BS_MUMU: QCDOrder.NNLO,
            ...     Observables.BR_BD_MUMU: QCDOrder.NNLO,
            ... }, add_dependencies=True)
        """
        cpp_map = {_cpp_observable_enum(k): _cpp_qcd_order(v) for k, v in obs_names.items()}
        self._cpp_obj.add_observables(cpp_map, bool(add_dependencies))
        return self

    def add_observable_ids(
        self,
        obs_names: Mapping[ObservableId, QCDOrder],
        add_dependencies: bool = True,
    ) -> "ObservableInterface":
        """Register several id-based observables at once.

        Args:
            obs_names: Mapping from internal observable ids to requested QCD
                orders.
            add_dependencies: Whether to attach dependencies for every
                observable in the mapping.

        Returns:
            ``self`` to support fluent chaining.
        """
        cpp_map = {_cpp_observable_id(k): _cpp_qcd_order(v) for k, v in obs_names.items()}
        self._cpp_obj.add_observables(cpp_map, bool(add_dependencies))
        return self

    def add_observables_from_decay(
        self,
        decay: Decays,
        qcd_order: QCDOrder,
        add_dependencies: bool = True,
    ) -> "ObservableInterface":
        """Register every observable attached to a decay channel.

        Args:
            decay: Decay family whose observables should be registered.
            qcd_order: Maximum QCD order requested for all observables in the
                decay family.
            add_dependencies: Whether to attach dependencies for the selected
                observables.

        Returns:
            ``self`` to support fluent chaining.
        """
        self._cpp_obj.add_observables(_cpp_decay(decay), _cpp_qcd_order(qcd_order), bool(add_dependencies))
        return self

    def add_observable_parameter(self, obs: Observables, pid: ParamId) -> "ObservableInterface":
        """Add one explicit parameter dependency to an enum observable.

        Args:
            obs: Observable enum value.
            pid: Parameter identifier to mark as a dependency.

        Returns:
            ``self`` to support fluent chaining.
        """
        self._cpp_obj.add_observable_parameter(_cpp_observable_enum(obs), _cpp_param_id(pid))
        return self

    def add_observable_id_parameter(self, obs: ObservableId, pid: ParamId) -> "ObservableInterface":
        """Add one explicit parameter dependency to an id-based observable.

        Args:
            obs: Internal observable identifier.
            pid: Parameter identifier to mark as a dependency.

        Returns:
            ``self`` to support fluent chaining.
        """
        self._cpp_obj.add_observable_parameter(_cpp_observable_id(obs), _cpp_param_id(pid))
        return self

    def add_observable_parameters(self, obs: Observables, pids: Set[ParamId]) -> "ObservableInterface":
        """Add several dependencies to an enum observable.

        Args:
            obs: Observable enum value.
            pids: Set of parameter identifiers to attach.

        Returns:
            ``self`` to support fluent chaining.
        """
        self._cpp_obj.add_observable_parameters(_cpp_observable_enum(obs), {_cpp_param_id(p) for p in pids})
        return self

    def add_observable_id_parameters(self, obs: ObservableId, pids: Set[ParamId]) -> "ObservableInterface":
        """Add several dependencies to an id-based observable.

        Args:
            obs: Internal observable identifier.
            pids: Set of parameter identifiers to attach.

        Returns:
            ``self`` to support fluent chaining.
        """
        self._cpp_obj.add_observable_parameters(_cpp_observable_id(obs), {_cpp_param_id(p) for p in pids})
        return self

    def remove_observable(self, obs: Observables) -> "ObservableInterface":
        """Remove one enum observable from the current selection.

        Args:
            obs: Observable enum value to remove.

        Returns:
            ``self`` to support fluent chaining.
        """
        self._cpp_obj.remove_observable(_cpp_observable_enum(obs))
        return self

    def remove_observable_id(self, obs: ObservableId) -> "ObservableInterface":
        """Remove one id-based observable from the current selection.

        Args:
            obs: Internal observable identifier to remove.

        Returns:
            ``self`` to support fluent chaining.
        """
        self._cpp_obj.remove_observable(_cpp_observable_id(obs))
        return self

    def remove_observables(self, ids: Set[Observables]) -> "ObservableInterface":
        """Remove several enum observables from the current selection.

        Args:
            ids: Set of observable enum values to remove.

        Returns:
            ``self`` to support fluent chaining.
        """
        self._cpp_obj.remove_observables({_cpp_observable_enum(x) for x in ids})
        return self

    def remove_observables_from_decay(self, decay: Decays) -> "ObservableInterface":
        """Remove every selected observable associated with a decay family.

        Args:
            decay: Decay family whose observables should be removed.

        Returns:
            ``self`` to support fluent chaining.
        """
        self._cpp_obj.remove_observables(_cpp_decay(decay))
        return self

    def compute_observable(self, obs: Observables) -> List[ObservableValue]:
        """Compute the theory prediction for one enum observable.

        Args:
            obs: Observable enum value to compute.

        Returns:
            List of observable values. Unbinned observables usually return a
            single entry. Binned observables return one entry per bin.
        """
        return [ObservableValue.from_cpp(v) for v in self._cpp_obj.compute_observable(_cpp_observable_enum(obs))]

    def compute_observable_id(self, obs: ObservableId) -> List[ObservableValue]:
        """Compute the theory prediction for one id-based observable.

        Args:
            obs: Internal observable identifier to compute.

        Returns:
            List of observable values. Binned observables can contain multiple
            entries with non-null ``ObservableValue.bin``.
        """
        return [ObservableValue.from_cpp(v) for v in self._cpp_obj.compute_observable(_cpp_observable_id(obs))]

    def compute_observable_central(self, obs: Observables) -> float:
        """Compute a scalar central value for an unbinned enum observable.

        Args:
            obs: Observable enum value to compute.

        Returns:
            The single predicted value.

        Raises:
            RuntimeError: If the C++ layer returns no value or more than one
                value, which indicates that the observable is binned and should
                be read with :meth:`compute_observable`.
        """
        vals = self.compute_observable(obs)
        if not vals:
            raise RuntimeError("compute_observable a renvoyé une liste vide.")
        if len(vals) != 1:
            raise RuntimeError("Observable binné : utilisez compute_observable() pour récupérer toutes les bins.")
        return vals[0].value

    def compute_observable_id_central(self, obs: ObservableId) -> float:
        """Compute a scalar central value for an unbinned id-based observable.

        Args:
            obs: Internal observable identifier to compute.

        Returns:
            The single predicted value.

        Raises:
            RuntimeError: If the C++ layer returns no value or more than one
                value, which indicates that the observable is binned and should
                be read with :meth:`compute_observable_id`.
        """
        vals = self.compute_observable_id(obs)
        if not vals:
            raise RuntimeError("compute_observable_id a renvoyé une liste vide.")
        if len(vals) != 1:
            raise RuntimeError("Observable binné : utilisez compute_observable_id() pour récupérer toutes les bins.")
        return vals[0].value

    def compute_all(self) -> Dict[ObservableId, List[ObservableValue]]:
        """Compute all currently registered observables.

        Returns:
            Mapping from observable id to a list of predicted values. Each list
            may contain several bins.
        """
        cpp_map = self._cpp_obj.compute_all()
        return {
            ObservableId(str(cpp_id)): [ObservableValue.from_cpp(v) for v in values]
            for cpp_id, values in cpp_map.items()
        }

    def get_exp_value(self, obs: Observables) -> float:
        """Return the experimental central value for an enum observable.

        Args:
            obs: Observable enum value.

        Returns:
            Experimental central value stored in the C++ observable database.
        """
        return float(self._cpp_obj.get_exp_value(_cpp_observable_enum(obs)))

    def get_exp_value_id(self, obs: ObservableId) -> float:
        """Return the experimental central value for an id-based observable.

        Args:
            obs: Internal observable identifier.

        Returns:
            Experimental central value stored in the C++ observable database.
        """
        return float(self._cpp_obj.get_exp_value(_cpp_observable_id(obs)))

    def get_exp_uncertainty(
        self,
        obs: Observables,
        u_type: UncertaintyType = UncertaintyType.COMBINED,
    ) -> float:
        """Return an experimental uncertainty for an enum observable.

        Args:
            obs: Observable enum value.
            u_type: Type of uncertainty to retrieve. The default is the combined
                uncertainty.

        Returns:
            Requested experimental uncertainty.
        """
        return float(self._cpp_obj.get_exp_uncertainty(_cpp_observable_enum(obs), _cpp_uncertainty_type(u_type)))

    def get_exp_uncertainty_id(
        self,
        obs: ObservableId,
        u_type: UncertaintyType = UncertaintyType.COMBINED,
    ) -> float:
        """Return an experimental uncertainty for an id-based observable.

        Args:
            obs: Internal observable identifier.
            u_type: Type of uncertainty to retrieve. The default is the combined
                uncertainty.

        Returns:
            Requested experimental uncertainty.
        """
        return float(self._cpp_obj.get_exp_uncertainty(_cpp_observable_id(obs), _cpp_uncertainty_type(u_type)))

    def get_current_observables(self) -> List[BinnedObservableId]:
        """Return the currently registered observable/bin identifiers.

        Returns:
            List of ``BinnedObservableId`` objects currently selected in the
            C++ observable manager.
        """
        return [BinnedObservableId.from_cpp(x) for x in self._cpp_obj.get_current_observables()]

    def get_all_ops_deps(self, obs: Observables) -> Set[ParamId]:
        """Return all known parameter dependencies for an enum observable.

        Args:
            obs: Observable enum value.

        Returns:
            Set of parameter identifiers allowed as dependencies for this
            observable.
        """
        return {_param_from_cpp(x) for x in self._cpp_obj.get_all_ops_deps(_cpp_observable_enum(obs))}

    def get_all_ops_deps_id(self, obs: ObservableId) -> Set[ParamId]:
        """Return all known parameter dependencies for an id-based observable.

        Args:
            obs: Internal observable identifier.

        Returns:
            Set of parameter identifiers allowed as dependencies for this
            observable.
        """
        return {_param_from_cpp(x) for x in self._cpp_obj.get_all_ops_deps(_cpp_observable_id(obs))}

    def set_decay_config(self, decay: Decays, config: DecayConfig) -> "ObservableInterface":
        """Set the concrete configuration object used by one decay engine.

        Args:
            decay: Decay family to configure.
            config: Python decay-configuration wrapper. The concrete wrapper class
                must match the selected decay engine, for example ``BKllConfig`` for
                ``Decays.B__K_l_l`` or ``BKstarllConfig`` for
                ``Decays.B__Kstar_l_l``.

        Returns:
            ``self`` to support fluent chaining.
        """
        if not isinstance(config, DecayConfig):
            raise TypeError(f"config must inherit DecayConfig, got {type(config).__name__}.")
        self._cpp_obj.set_decay_config(_cpp_decay(decay), config.to_cpp())
        return self


    def set_param(self, pid: ParamId, value: float) -> None:
        """Set one model parameter in the C++ parameter store.

        Args:
            pid: Parameter identifier. ``pid.type`` must be defined and
                ``pid.code`` must contain exactly one LHA code entry.
            value: New parameter value.

        Raises:
            ValueError: If ``pid.type`` is missing or if the LHA code is not a
                single integer code.
            TypeError: If ``pid`` or ``pid.type`` has an invalid Python type.

        Notes:
            After mutating parameters, call :meth:`reload_params` and possibly
            :meth:`enable_obs` if the selected decays cache input parameters.
        """
        _require(pid, ParamId, "pid")
        if pid.type is None:
            raise ValueError("pid.type doit être défini pour set_param().")
        _require(pid.type, ParameterType, "pid.type")
        self._cpp_obj.set_param(str(pid.block), _single_lha_code(pid.code), float(value), pid.type.value)

    def get_param(self, pid: ParamId):
        """Read one model parameter from the C++ parameter store.

        Args:
            pid: Parameter identifier. ``pid.type`` must be defined and
                ``pid.code`` must contain exactly one LHA code entry.

        Returns:
            Parameter value returned by the C++ parameter provider.

        Raises:
            ValueError: If ``pid.type`` is missing or if the LHA code is not a
                single integer code.
            TypeError: If ``pid`` or ``pid.type`` has an invalid Python type.
        """
        _require(pid, ParamId, "pid")
        if pid.type is None:
            raise ValueError("pid.type doit être défini pour get_param().")
        _require(pid.type, ParameterType, "pid.type")
        return self._cpp_obj.get_param(str(pid.block), _single_lha_code(pid.code), pid.type.value)

    def reload_params(self) -> None:
        """Reload cached parameters for every registered decay.

        This forwards to the C++ manager and is typically used after calling
        :meth:`set_param`.
        """
        self._cpp_obj.reload_params()

    def enable_obs(self) -> None:
        """Force the C++ manager to enable/re-enable selected observables.

        This can rebuild or refresh the internal observable state after a larger
        configuration change.
        """
        self._cpp_obj.enable_obs()

    def set_bkstarll_threads(self, n_threads: int) -> None:
        """Configure the number of threads used by the ``B -> K* ll`` decay.

        Args:
            n_threads: Number of worker threads requested by the C++ decay
                implementation.
        """
        self._cpp_obj.set_bkstarll_threads(int(n_threads))

    def set_bkll_threads(self, n_threads: int) -> None:
        """Configure the number of threads used by the ``B -> K ll`` decay.

        Args:
            n_threads: Number of worker threads requested by the C++ decay
                implementation.
        """
        self._cpp_obj.set_bkll_threads(int(n_threads))

    def set_bsphi_threads(self, n_threads: int) -> None:
        """Configure the number of threads used by the ``Bs -> phi ll`` decay.

        Args:
            n_threads: Number of worker threads requested by the C++ decay
                implementation.
        """
        self._cpp_obj.set_bsphi_threads(int(n_threads))


__all__ = ["ObservableInterface"]

    
if __name__ == "__main__":
    from pyhyperiso.core.Core.HyperisoMaster import HyperisoMaster
    from pathlib import Path
    from pyhyperiso.core.Core.HyperisoConfig import HyperisoConfig, ExternalFlag
    from pyhyperiso.core.Common.GeneralEnum import Model
    from pyhyperiso.core.Core.ParamaterProvider import ParameterProvider
    print("🔧 Initializing PyHyperisoMaster with custom PyHyperisoConfig...")

    config = HyperisoConfig(
        flags={
            ExternalFlag.IS_LHA_SPECTRUM: True,
            ExternalFlag.HAS_WILSON_INPUT: False,
            ExternalFlag.HAS_TH_OBSERVABLE_INPUT: False,
            ExternalFlag.HYP_AS_SM_MARTY: True,
            # ExternalFlag.USE_MARTY: False
        },
        model=Model.SM,
        mty_model_name="MSSM_UFO",
        mty_model_path=Path("/my/custom/marty/path")
    )

    print("🔧 PyHyperisoConfig content:")
    print(config)

    hyp = HyperisoMaster()
    lha_file_path = "lha/si_input.flha" 

    print("\n🚀 Calling init with config...")
    hyp.init(lha_file=lha_file_path, config=config)
    

    from pyhyperiso.core.BusinessLogic.DecayConfig import BKstarllConfig
    
    kstar_conf = BKstarllConfig()
    
    interface = ObservableInterface()
    interface.add_observable(Observables.BR_B_XS_GAMMA, QCDOrder.NNLO, True)
    interface.add_observable(Observables.BR_BS_MUMU, QCDOrder.NNLO, True)
    interface.add_observable(Observables.BR_BD_MUMU, QCDOrder.NNLO, True)
    interface.add_observable(Observables.BR_B0__D_TAU_NU, QCDOrder.NNLO, True)
    interface.add_observable(Observables.BR_B__KSTAR_GAMMA, QCDOrder.NNLO, True)
    interface.set_decay_config(Decays.B__Kstar_l_l, kstar_conf)
    print(interface.compute_observable(Observables.BR_B_XS_GAMMA))  # Scalar(...)
    print(interface.compute_observable(Observables.BR_BS_MUMU))  # Scalar(...)
    print(interface.compute_all())
    
    test_values = []
    
    from pyhyperiso.core.Core.ParameterSetter import ParameterSetter, ParamId
    py_set = ParameterSetter()
    for i in range(50, 90, 10):
        py_set.mutate(ParamId(ParameterType.SM, "MASS", 24), i)
        print("aah : ", ParameterProvider().get_by_pid(ParamId(ParameterType.SM, "MASS", 24)))
        # interface.enable_obs()
        print("aah : ", ParameterProvider().get_by_pid(ParamId(ParameterType.SM, "MASS", 24)))
        test_values.append(interface.compute_observable(Observables.BR_B_XS_GAMMA)[0].value)
    
    import matplotlib
    matplotlib.use("TkAgg")
    import matplotlib.pyplot as plt
    
    plt.scatter(range(50, 90, 10), test_values)
    
    plt.show()
    # test_values = []
    # from pyhyperiso.core.Core.ParameterSetter import PyParameterSetter, ParamId, ParameterType
    # py_set = PyParameterSetter()
    # for i in range(1, 81):
    #     py_set.mutate(ParamId(ParameterType.WILSON, "B_SCALE", 1), i)
    #     test_values.append(interface.get_R(WilsonRequest(WGroup.B, WCoeff.C9, QCDOrder.LO, ContributionType.TOTAL)))
    
    # import matplotlib
    # matplotlib.use("TkAgg")
    # import matplotlib.pyplot as plt
    
    # plt.scatter(range(1, 81), test_values)
    
    # plt.show()