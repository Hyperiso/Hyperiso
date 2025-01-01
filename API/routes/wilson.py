from fastapi import APIRouter, Query

router = APIRouter()

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
