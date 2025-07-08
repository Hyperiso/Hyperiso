# Hyperiso/api/scan.py
from fastapi import APIRouter

router = APIRouter()

@router.get("/mock")
def get_mock_scan():
    # Exemple simple : param ↦ observable
    x_vals = [i for i in range(40, 120)]
    y_vals = [0.1 * x + 5 for x in x_vals]
    return {"x": x_vals, "y": y_vals}
