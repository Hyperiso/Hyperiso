# Hyperiso/api/parameters.py
from fastapi import APIRouter
from pyhyperiso.core.Core.ParamaterProvider import PyParameterProvider
# from pyhyperiso.core.Core.ParameterProvider import PyParameterProvider

router = APIRouter()

@router.get("/all")
def get_all_params():
    provider = PyParameterProvider()
    data = []
    for block in ["MASS", "SMINPUTS"]:
        for code in range(1, 50):
            try:
                val = provider.get_by_block(block, [code])
                data.append({"block": block, "code": code, "value": val})
            except Exception:
                continue
    return data