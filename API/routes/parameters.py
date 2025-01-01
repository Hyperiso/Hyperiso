from fastapi import APIRouter, UploadFile, File, HTTPException
import shutil
from pathlib import Path
import sys, os
sys.path.append(os.path.join(os.getcwd(), "Python"))
from Python.Phyperiso import Parameters, ParameterType, Model, MemoryManager

router = APIRouter()

LHA_DIR = Path("DataBase/lha/")
LHA_DIR.mkdir(parents=True, exist_ok=True)

mem = MemoryManager()
mem.init("Test/InputFiles/testInput.flha", Model.SM)
parameters = Parameters(ParameterType.SM)

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
    if not parameters.exists(block, code):
        raise HTTPException(status_code=404, detail=f"Parameter {block}/{code} not found")
    value = parameters(block, code)
    return {"block": block, "code": code, "value": value}
