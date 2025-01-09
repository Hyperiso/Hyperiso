from fastapi import APIRouter, Query
from Python.Phyperiso import WilsonManager, Model, ParameterType
from API.utils.ParametersCache import MemoryManagerCache, ParametersCache

router = APIRouter()
mmCache = MemoryManagerCache("Test/InputFiles/testInput.flha", model=Model.SM)
import os

wilson_manager = WilsonManager("SM")
parameters = ParametersCache(ParameterType.SM)
mock_coefficients = {
    "C1": {"matching_value": 0.12, "running_value": 0.15},
    "C2": {"matching_value": 0.34, "running_value": 0.30},
}

@router.post("/set_matching_scale")
def set_matching_scale(scale: float):
    """Define matching scale for wilson coefficient."""
    return {"message": f"Matching scale set to {scale} GeV"}

@router.get("/coefficient")
def get_coefficient(name: str):
    """Retrieve coefficient value."""
    if name not in mock_coefficients:
        return {"error": f"Coefficient {name} not found"}
    return mock_coefficients[name]

@router.get("/plot_variation")
def plot_coefficient_variation(
    param_name: str,
    min_value: float,
    max_value: float,
    steps: int = 10
):
    """Obtain variation of coefficient given one parameters."""
    values = [{"param": v, "value": v * 0.1} for v in range(int(min_value), int(max_value + 1))]
    return {"param_name": param_name, "values": values}

@router.get("/calculate_coefficient")
def calculate_coefficient(group: str, name: str, order: str):
    """Calculate a specific coefficient."""
    coefficient = wilson_manager.get_matching_coefficient(group, name, order)
    return {"group": group, "name": name, "value": coefficient}

@router.get("/calculate_all_lhas")
def calculate_all_lhas(group: str, name: str, order: str):
    """Calculate a coefficient for all LHAs in the directory."""
    results = []
    for lha in os.listdir("DataBase/lha"):
        mmCache.switch_lha(lha, Model.SM)
        coefficient = wilson_manager.get_matching_coefficient(group, name, order)
        results.append({"lha": lha, "value": coefficient})
    return {"results": results}


@router.get("/plot_coefficients")
def plot_coefficients(group: str, name: str, order: str, param_name: str, min_value: float, max_value: float, steps: int = 10):
    """Plot coefficient variation."""
    values = []
    for step in range(steps):
        value = min_value + step * (max_value - min_value) / (steps - 1)
        parameters.set_block_value(param_name, 0, value)  # Example with pdgcode = 0
        coefficient = wilson_manager.get_matching_coefficient(group, name, order)
        values.append({"param": value, "coefficient": coefficient})
    return {"param_name": param_name, "values": values}