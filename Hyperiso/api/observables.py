# Hyperiso/api/observables.py
from fastapi import APIRouter

router = APIRouter()

@router.get("")
def get_observables():
    return [
        {"name": "BR(Bs -> mu mu)", "value": 3.5e-9, "uncertainty": 0.2e-9},
        {"name": "Delta M_s", "value": 17.8, "uncertainty": 0.1},
        {"name": "RK", "value": 0.85, "uncertainty": 0.03},
    ]
