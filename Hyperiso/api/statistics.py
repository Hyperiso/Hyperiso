# Hyperiso/api/statistics.py
from fastapi import APIRouter

router = APIRouter()

@router.get("/chi2")
def get_chi2():
    chi2 = 12.3
    dof = 6
    return {"chi2": chi2, "dof": dof}
