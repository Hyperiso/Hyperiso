from fastapi import APIRouter, Query
from Python.Phyperiso import WilsonManager, Model, ParameterType
from Python.Phyperiso import BCoefficientGroup, BScalarCoefficientGroup, BPrimeCoefficientGroup
from API.utils.ParametersCache import MemoryManagerCache, ParametersCache
from pydantic import BaseModel

router = APIRouter()
mmCache = MemoryManagerCache("Test/InputFiles/testInput.flha", model=Model.SM)
import os

wilson_managers = {"SM" : WilsonManager("SM"), "SUSY" : WilsonManager("SUSY"),
                   "THDM" : WilsonManager("THDM"), "CUSTOM" : WilsonManager("CUSTOM")}
wilson_manager = WilsonManager("SM")
parameters = ParametersCache(ParameterType.SM)


map_group = {"BCoefficientGroupSM" : BCoefficientGroup(),
             "BScalarCoefficientGroupSM" : BScalarCoefficientGroup(),
             "BPrimeCoefficientGroupSM" : BPrimeCoefficientGroup()}

class GroupRequest(BaseModel):
    group : str = "BCoefficientGroup"
    model : str = "SM"

class WilsonRequest(BaseModel):
    model : str = "SM"
    group: str = "BCoefficientGroup"
    scale : float = 81
    qcd_order : str = "LO"
    name : str = "C1"

class PlotWilsonRequest(BaseModel):
    model : str = "SM"
    group : str = "BCoefficientGroup"
    name : str = "C1"
    order : str = "LO"
    matching_scale : float = 81
    running_scale : float = 42
    param_block : str = "SMINPUTS"
    param_pdgcode : int = 6
    min_value : float = 150
    max_value : float = 200
    steps : int = 10

@router.post("/register_group")
def register_group(request : GroupRequest):
    wilson_manager.register_coefficient_group(request.group, map_group[request.group + request.model])
    return {"message" : f"group {request.group} initialized for model {request.model}"}

@router.post("/set_matching_scale")
def set_matching_scale(request: WilsonRequest):
    """Define matching scale for wilson coefficient."""
    wilson_managers[request.model].set_q_match(request.group, request.scale)
    wilson_managers[request.model].set_matching_coefficient(request.group, request.qcd_order)
    return {"message": f"Matching scale set to {request.scale} GeV"}

@router.post("/set_group_scale")
def set_matching_scale(request: WilsonRequest):
    """Define matching scale for wilson coefficient."""
    wilson_managers[request.model].set_group_scale(request.group, request.scale)
    wilson_managers[request.model].set_run_coefficient(request.group, request.qcd_order)
    return {"message": f"Running scale set to {request.scale} GeV"}

@router.get("/get_coefficient")
def get_coefficient(request : WilsonRequest):
    """Retrieve coefficient value."""
    print("I was here")
    wilson_managers[request.model].set_q_match(request.group, request.scale)
    wilson_managers[request.model].set_matching_coefficient(request.group, request.qcd_order)
    value = wilson_managers[request.model].get_matching_coefficient(request.group, request.name, request.qcd_order)
    return {"coeff_real" : value.real, "coeff_img" : value.imag}


@router.get("/get_run_coefficient")
def get_coefficient(model : str, group : str, name: str, order : str):
    """Retrieve coefficient value."""
    print("I was here")
    value = wilson_managers[model].get_run_coefficient(group, name, order)
    return {"coeff_real" : value.real, "coeff_img" : value.imag}


@router.get("/plot_variation")
def plot_coefficient_variation(
    param_name: str,
    min_value: float,
    max_value: float,
    steps: int = 10
):
    """Obtain variation of coefficient given one parameters."""
    values = [{"param": v, "value": v * 0.1} for v in range(int(min_value), int(max_value + 1))]
    return {"param_name": param_name, "values": values}

@router.get("/calculate_coefficient")
def calculate_coefficient(group: str, name: str, order: str):
    """Calculate a specific coefficient."""
    coefficient = wilson_manager.get_matching_coefficient(group, name, order)
    return {"group": group, "name": name, "value": coefficient}

@router.get("/calculate_all_lhas")
def calculate_all_lhas(group: str, name: str, order: str, scale : float):
    """Calculate a coefficient for all LHAs in the directory."""
    results = []
    for lha in os.listdir("DataBase/lha"):
        infos = {"lha": os.path.join("DataBase/lha", lha), "model" : Model.SM, "use_marty" : False,
             "is_spectrum" : False, "has_wilson" : False, "has_obs":  False}
        mmCache.switch_info(infos)
        wilson_managers["SM"].set_q_match(group, scale)
        wilson_managers["SM"].set_matching_coefficient(group, order)
        coefficient = wilson_manager.get_matching_coefficient(group, name, order)
        results.append({"lha": lha, "value": coefficient.real})
    return {"results": results}


@router.get("/plot_coefficients")
def plot_coefficients(request : PlotWilsonRequest):
    """Plot coefficient variation."""
    values = []
    for step in range(request.steps):
        value = request.min_value + step * (request.max_value - request.min_value) / (request.steps - 1)
        wilson_managers[request.model].set_params(request.group, request.param_block, request.param_pdgcode, value)
        # parameters.set_block_value(request.param_block, request.param_pdgcode, value)
        wilson_managers[request.model].set_q_match(request.group, request.matching_scale)
        wilson_managers[request.model].set_matching_coefficient(request.group, request.order)
        coefficient = wilson_managers[request.model].get_matching_coefficient(request.group, request.name, request.order)
        values.append({"param": value, "coefficient": coefficient.real, "coefficient_imag" : coefficient.imag})
    print(values)
    return {"param_name": request.param_block+"_"+str(request.param_pdgcode), "values": values}

@router.get("/plot_run_coefficients")
def plot_coefficients(request : PlotWilsonRequest):
    """Plot coefficient variation."""
    values = []
    for step in range(request.steps):
        value = request.min_value + step * (request.max_value - request.min_value) / (request.steps - 1)
        wilson_managers[request.model].set_params(request.group, request.param_block, request.param_pdgcode, value)
        # parameters.set_block_value(request.param_block, request.param_pdgcode, value)
        wilson_managers[request.model].set_q_match(request.group, request.matching_scale)
        wilson_managers[request.model].set_matching_coefficient(request.group, request.order)
        wilson_managers[request.model].set_group_scale(request.group, request.running_scale)
        wilson_managers[request.model].set_run_coefficient(request.group, request.order)
        coefficient = wilson_managers[request.model].get_run_coefficient(request.group, request.name, request.order)
        values.append({"param": value, "coefficient": coefficient.real, "coefficient_imag" : coefficient.imag})
    print(values)
    return {"param_name": request.param_block+"_"+str(request.param_pdgcode), "values": values}

@router.get("/plot_run_coefficients_Q")
def plot_coefficients(request : PlotWilsonRequest):
    """Plot coefficient variation."""
    values = []
    wilson_managers[request.model].set_q_match(request.group, request.matching_scale)
    wilson_managers[request.model].set_matching_coefficient(request.group, request.order)
    wilson_managers[request.model].set_group_scale(request.group, request.running_scale)
    wilson_managers[request.model].set_run_coefficient(request.group, request.order)
    for step in range(request.steps):
        value = request.min_value + step * (request.max_value - request.min_value) / (request.steps - 1)
        wilson_managers[request.model].set_group_scale(request.group, value)
        # wilson_managers[request.model].set_run_coefficient(request.group, request.order)
        coefficient = wilson_managers[request.model].get_run_coefficient(request.group, request.name, request.order)
        values.append({"param": value, "coefficient": coefficient.real, "coefficient_imag" : coefficient.imag})
    return {"param_name": request.param_block+"_"+str(request.param_pdgcode), "values": values}