from fastapi import APIRouter, UploadFile, File, HTTPException
import shutil
from pathlib import Path
import sys, os
sys.path.append(os.path.join(os.getcwd(), "Python"))
from Python.Phyperiso import Parameters, ParameterType, Model, MemoryManager
from API.utils.ParametersCache import MemoryManagerCache, ParametersCache

mmCache = MemoryManagerCache("Test/InputFiles/testInput.flha", model=Model.SM)

router = APIRouter()

LHA_DIR = Path("DataBase/lha/")
LHA_DIR.mkdir(parents=True, exist_ok=True)

mem_cache = MemoryManagerCache("Test/InputFiles/testInput.flha", Model.SM)

param_cache = ParametersCache(ParameterType.SM)


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
    if not param_cache.exists(block, code):
        raise HTTPException(status_code=404, detail=f"Parameter {block}/{code} not found")
    value = param_cache(block, code)
    return {"block": block, "code": code, "value": value}

@router.get("/lha")
def get_lha():
    return {"lha" : mem_cache.get_lha()}

@router.get("/blocks_list")
def get_blocks():
    return {"blocks" : mem_cache.get_blocks_list()}

@router.get("/block_info")
def get_block_info(block : str):
    if not block in mem_cache.get_blocks_list():
        raise HTTPException(status_code=404, detail = f"Block {block} not found")
    return {block : mem_cache.get_block_infos(block)}


