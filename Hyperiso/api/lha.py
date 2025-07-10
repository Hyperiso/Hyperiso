# Hyperiso/api/lha.py
from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from pyhyperiso.core.Core.HyperisoMaster import PyHyperisoMaster
from pyhyperiso.core.Core.Config import PyConfig, ExternalFlag
from pyhyperiso.core.Common.GeneralEnum import Model
from pathlib import Path

router = APIRouter()

hyperiso = PyHyperisoMaster()

DEFAULT_CONFIG = PyConfig(
    flags={
        ExternalFlag.IS_LHA_SPECTRUM: True,
        ExternalFlag.HAS_WILSON_INPUT: False,
        ExternalFlag.HAS_TH_OBSERVABLE_INPUT: False,
        ExternalFlag.USE_MARTY: False
    },
    model=Model.SM,
    mty_model_name="MSSM_UFO",
    mty_model_path=Path("/my/custom/marty/path")
)

class LhaInput(BaseModel):
    path: str

@router.post("/analyse")
def analyse_lha(req: LhaInput):
    try:
        hyperiso.init(lha_file=req.path, config=DEFAULT_CONFIG)
        return {
            "status": "LHA chargé avec succès",
            "model": hyperiso.model.name
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@router.post("/switch")
def switch_lha(req: LhaInput):
    try:
        hyperiso.switch_lha(lha_file=req.path, config=DEFAULT_CONFIG)
        return {
            "status": "LHA changé avec succès",
            "model": hyperiso.model.name
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@router.get("/status")
def get_current_model():
    try:
        model = hyperiso.model.name
        return {"status": "actif", "model": model}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
