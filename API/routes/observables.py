from fastapi import APIRouter, Query

router = APIRouter()

mock_observables = {
    "BR_BS_MUMU": 3.4e-9,
    "BR_BD_MUMU": 1.2e-10,
}

@router.get("/calculate")
def calculate_observable(name: str):
    """Observable calculation."""
    if name not in mock_observables:
        return {"error": f"Observable {name} not found"}
    return {"name": name, "value": mock_observables[name]}

@router.get("/plot_variation")
def plot_observable_variation(
    param_name: str,
    min_value: float,
    max_value: float,
    steps: int = 10
):
    """Parameters variance."""
    values = [{"param": v, "value": v * 0.05} for v in range(int(min_value), int(max_value + 1))]
    return {"param_name": param_name, "values": values}

@router.get("/chi2")
def calculate_chi2():
    """Chi2 calculation."""
    chi2_value = 12.34
    return {"value": chi2_value}
