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

model_available = {"SM"}

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
    "SM" : Model.SM,
    "SUSY" : Model.SUSY,
    "THDM" : Model.THDM,
    "CUSTOM" : Model.CUSTOM
}

class SetLHAModelRequest(BaseModel):
    lha_file: str = None
    model: str = None
    use_marty: bool = None
    is_spectrum: bool = None
    has_wilsons: bool = None
    has_obs: bool = None

def check_model_lha(lha : str, model : str) -> bool:
    with open(lha, "r") as f:
        data = f.readlines()

    if model == "SUSY":
        mandatory = ["MODSEL", "SMINPUTS", "MINPAR", "EXTPAR"]
    elif model == "THDM":
        mandatory = ["MODSEL", "SMINPUTS", "MINPAR", "MASS"]
    else:
        return True
    is_okay = [False,False,False,False]
    for line in data:
        for i, elem in enumerate(mandatory):
            if elem in line:
                is_okay[i] = True
    
    if False in is_okay:
        return False
    else:
        return True
    

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
        return {"block" : block, "pdgcode" : code, "value" : ""}
        raise HTTPException(status_code=404, detail=f"Parameter {block}/{code} not found")
    value = param_cache(block, code)
    return {"block": block, "code": code, "value": value}

@router.get("/lha")
def get_lha():
    return {"lha" : mem_cache.get_lha()}

@router.get("/blocks_list")
def get_blocks(param_type : str):
    print("trying to get a block list")
    if param_type in map_model.keys():
        if param_type not in model_available:
            return {"blocks" : []}
    return {"blocks" : mem_cache.get_blocks_list(map_paramtype[param_type])}

@router.get("/block_info")
def get_block_info(block : str, param_type : str):
    if not block in mem_cache.get_blocks_list(map_paramtype[param_type]):
        raise HTTPException(status_code=404, detail = f"Block {block} not found")
    return {block : mem_cache.get_block_infos(block, map_paramtype[param_type])}

from fastapi import Request

@router.post("/set_memory_manager")
def set_lha_model(request: SetLHAModelRequest):
    print("Requête reçue : ", request.dict())
    infos = {"lha": request.lha_file, "model" : map_model[request.model], "use_marty" : request.use_marty,
             "is_spectrum" : request.is_spectrum, "has_wilson" : request.has_wilsons, "has_obs":  request.has_obs}
    if not check_model_lha(request.lha_file, request.model):
        print("not compatible : ", request.lha_file, request.model)
        return {"message", "LHA no switched, model and lha not compatible."}
    
    mem_cache.switch_info(infos)
    param_cache.switch_param(map_paramtype[request.model])
    model_available.add(request.model)
    return {"message": f"LHA file switched to {request.lha_file} with model {request.model}"}

@router.get("/model_enum")
def model_enum(model: str):
    print("wtf")
    if model == "SM":
        return {"model" : Model.SM}
    
@router.post("/set_parameter")
def set_parameter(block: str, code: int, value: float):
    print("trying to set a param")
    """Set a parameter value in a specific block."""
    param_cache.set_block_value(block, code, value, force=True)
    return {"message": f"Parameter set: {block}[{code}] = {value}"}
