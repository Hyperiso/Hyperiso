"""Regression tests for installed-wheel statistical path resolution."""

from __future__ import annotations

import subprocess
import sys
import textwrap


def test_statistic_interface_uses_active_runtime_nuisance_paths(tmp_path):
    """A fresh process must resolve nuisance data from the active package assets.

    HyperIso owns process-global C++ runtime state.  Running this regression in a
    subprocess prevents earlier tests from leaving ``MemoryManager`` initialized
    with unrelated dummy blocks or path providers.
    """
    program = textwrap.dedent(
        r"""
        from pathlib import Path

        import pyhyperiso
        from pyhyperiso.Common import BinnedObservableId, Model, Observables, QCDOrder
        from pyhyperiso.Core import HyperisoConfig, HyperisoMaster
        from pyhyperiso.Observable import ObservableInterface
        from pyhyperiso.Statistic import ExperimentObs, StatisticConfig, StatisticInterface

        assets_root = Path(pyhyperiso.__file__).resolve().parent / "assets"
        lha_path = assets_root / "lha" / "si_input.flha"
        nuisance_path = assets_root / "default" / "nuisances.json"

        assert lha_path.is_file(), lha_path
        assert nuisance_path.is_file(), nuisance_path

        # Keep the master alive while interfaces use the active runtime providers.
        hyperiso = HyperisoMaster()
        hyperiso.init(
            lha_file=str(lha_path),
            config=HyperisoConfig(model=Model.SM),
        )

        observables = ObservableInterface()
        observables.add_observable(Observables.BR_BS_MUMU, QCDOrder.NNLO)

        statistic = StatisticInterface(
            StatisticConfig(MC_draws=4, MC_threads=1),
            observables,
        )

        assert statistic is not None

        exact = ExperimentObs(
            "DEFAULT",
            BinnedObservableId(Observables.BR_BS_MUMU),
        )
        statistic.select_experiment_observables([exact])
        assert statistic.has_experiment_observable_selection()

        selected = statistic.selected_experiment_observables()
        assert len(selected) == 1
        assert selected[0].experiment == "DEFAULT"
        assert selected[0].obs == exact.obs

        statistic.select_experiment_observables_all()
        assert not statistic.has_experiment_observable_selection()

        print("STATISTIC_RUNTIME_PATH_OK")
        print("EXACT_EXPERIMENT_SELECTION_OK")
        """
    )

    completed = subprocess.run(
        [sys.executable, "-c", program],
        cwd=tmp_path,
        check=False,
        capture_output=True,
        text=True,
    )
    output = completed.stdout + completed.stderr

    assert completed.returncode == 0, (
        "StatisticInterface failed in a fresh Python process.\n"
        f"stdout:\n{completed.stdout}\n"
        f"stderr:\n{completed.stderr}"
    )
    assert "STATISTIC_RUNTIME_PATH_OK" in completed.stdout
    assert "EXACT_EXPERIMENT_SELECTION_OK" in completed.stdout
    assert "/project/Assets/default/nuisances.json" not in output
