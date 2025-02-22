from fastapi import APIRouter, Query
from Hyperiso.Phyperiso import ObservableInterface, ObservableMapper, OrderMapper
router = APIRouter()

mock_observables = {
    "BR_BS_MUMU": 3.4e-9,
    "BR_BD_MUMU": 1.2e-10,
}

obs_interface = ObservableInterface()

@router.get("/compute_observable")
def calculate_observable(name: str):
    """Observable calculation."""
    if name not in ObservableMapper.get_model_str_list():
        return {"error": f"Observable {name} not found"}
    return {"name": name, "value": obs_interface.compute_observable(ObservableMapper.enum_elt(name))}

@router.post("/add_observable")
def add_observable(obs : str, order : str):
    obs_enum = ObservableMapper.enum_elt(obs)
    order_enum = OrderMapper.enum_elt(order)

    obs_interface.add_observable(obs_enum, order_enum)
    return {"message": f"Observable {obs} added in ObservableInterface with order {order}"}

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
