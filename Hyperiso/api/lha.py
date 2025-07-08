# Hyperiso/api/lha.py
from fastapi import APIRouter
from pydantic import BaseModel
from pyhyperiso.core.Core.HyperisoMaster import PyHyperisoMaster
from pyhyperiso.core.Core.Config import PyConfig, ExternalFlag
from pyhyperiso.core.Common.GeneralEnum import Model
from pathlib import Path

router = APIRouter()
hyperiso = PyHyperisoMaster()

class LhaInput(BaseModel):
    path: str

@router.post("/analyse")
def analyse_lha(req: LhaInput):
    config = PyConfig(
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
    hyperiso.init(lha_file=req.path, config=config)
    return {"status": "LHA chargé", "model": hyperiso.model.name}
