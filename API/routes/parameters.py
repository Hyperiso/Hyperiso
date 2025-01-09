from fastapi import APIRouter, UploadFile, File, HTTPException
import shutil
from pathlib import Path
import sys, os
sys.path.append(os.path.join(os.getcwd(), "Python"))
from Python.Phyperiso import Parameters, ParameterType, Model, MemoryManager
from API.utils.ParametersCache import MemoryManagerCache, ParametersCache
from pydantic import BaseModel

# mmCache = MemoryManagerCache("Test/InputFiles/testInput.flha", model=Model.SM)

router = APIRouter()

LHA_DIR = Path("DataBase/lha/")
LHA_DIR.mkdir(parents=True, exist_ok=True)
model = Model.SM
param_type = ParameterType.SM

mem_cache = MemoryManagerCache("Test/InputFiles/testInput.flha", model)

param_cache = ParametersCache(param_type)

map_paramtype = {
    "SM" : ParameterType.SM,
    "SUSY" : ParameterType.SUSY,
    "THDM" : ParameterType.THDM,
    "CUSTOM" : ParameterType.CUSTOM,
    "FLAVOR" : ParameterType.FLAVOR,
    "FF" : ParameterType.FF,
    "WILSON" : ParameterType.WILSON
}

map_model = {
    "SM" : ParameterType.SM,
    "SUSY" : ParameterType.SUSY,
    "THDM" : ParameterType.THDM
}

class SetLHAModelRequest(BaseModel):
    lha_file: str
    model: str
    use_marty: bool = False
    is_spectrum: bool = False
    has_wilsons: bool = False
    has_obs: bool = False

@router.post("/upload")
async def upload_lha(file: UploadFile = File(...)):
    """SLHA file or folder upload."""
    file_path = LHA_DIR / file.filename
    with open(file_path, "wb") as buffer:
        shutil.copyfileobj(file.file, buffer)
    return {"message": f"File {file.filename} uploaded successfully", "path": str(file_path)}

@router.get("/value")
def get_parameter_value(block: str, code: int):
    """Retrieve parameters value."""
    print("BLLOOOOOOK", block, code)
    if not param_cache.exists(block, code):
        raise HTTPException(status_code=404, detail=f"Parameter {block}/{code} not found")
    value = param_cache(block, code)
    return {"block": block, "code": code, "value": value}

@router.get("/lha")
def get_lha():
    return {"lha" : mem_cache.get_lha()}

@router.get("/blocks_list")
def get_blocks(param_type : str):
    return {"blocks" : mem_cache.get_blocks_list(map_paramtype[param_type])}

@router.get("/block_info")
def get_block_info(block : str, param_type : str):
    if not block in mem_cache.get_blocks_list(map_paramtype[param_type]):
        raise HTTPException(status_code=404, detail = f"Block {block} not found")
    return {block : mem_cache.get_block_infos(block, map_paramtype[param_type])}

from fastapi import Request

@router.post("/set_lha_model")
def set_lha_model(request: SetLHAModelRequest):
    print("Requête reçue : ", request.dict())
    if request.model == "SM":
        mem_cache.switch_lha(
            request.lha_file,
            map_model[request.model],
            request.use_marty,
            request.is_spectrum,
            request.has_wilsons,
            request.has_obs
        )
    param_cache.switch_param(ParameterType.SM)
    return {"message": f"LHA file switched to {request.lha_file} with model {request.model}"}

@router.get("/model_enum")
def model_enum(model: str):
    if model == "SM":
        return {"model" : Model.SM}
    
@router.post("/set_parameter")
def set_parameter(block: str, code: int, value: float):
    """Set a parameter value in a specific block."""
    param_cache.set_block_value(block, code, value, force=True)
    return {"message": f"Parameter set: {block}[{code}] = {value}"}
