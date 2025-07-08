# Hyperiso/api/wilson.py
from fastapi import APIRouter
from pydantic import BaseModel
from pyhyperiso.core.Common.GeneralEnum import WGroup, WCoeff, QCDOrder, ContributionType
from pyhyperiso.core.Common.Configs import PyWilsonRequest
from pyhyperiso.phyperiso.pyhyperiso.wilson.wilson_interface import WilsonInterface

router = APIRouter()
interface = WilsonInterface()

class WilsonQuery(BaseModel):
    group: str
    coefficient: str
    order: str
    contribution: str
    mu_W: float

@router.post("/get_M")
def get_m(req: WilsonQuery):
    interface.set_matching_scale(req.mu_W)
    wilson_req = PyWilsonRequest(
        group=WGroup[req.group],
        coefficient=WCoeff[req.coefficient],
        order=QCDOrder[req.order],
        contribution=ContributionType[req.contribution]
    )
    value = interface.get_M(wilson_req)
    return {"value": float(value)}
