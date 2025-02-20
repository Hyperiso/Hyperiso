from fastapi import HTTPException
from pydantic import BaseModel
import os
import json
from fastapi import APIRouter, Query

SLHA_DIRECTORY = "DataBase/lha/"
router = APIRouter()

if not os.path.exists(SLHA_DIRECTORY):
    os.makedirs(SLHA_DIRECTORY)

class ParameterConfig(BaseModel):
    block: str
    number: int
    start: float
    stop: float
    step: float

class SLHAGenerationRequest(BaseModel):
    template: str
    parameters: list[ParameterConfig]

def generate_slha_files(template: str, parameters: list[ParameterConfig]):
    generated_files = []

    def replace_parameter(template_lines, block, number, value):
        """
        Remplace un paramètre dans un fichier SLHA.

        Args:
            template_lines (list): Contenu du fichier SLHA sous forme de liste de lignes.
            block (str): Nom du bloc où effectuer la recherche (ex : "MODSEL").
            number (int): Numéro du paramètre à modifier.
            value (float): Nouvelle valeur à assigner.

        Returns:
            list: Les lignes mises à jour.
        """
        updated_lines = []
        inside_block = False

        for line in template_lines:
            stripped_line = line.strip()

            if stripped_line.upper().startswith(f"BLOCK {block.upper()}"):
                inside_block = True
                updated_lines.append(line)
                continue

            if inside_block:
                if stripped_line == "" or stripped_line.startswith("#"):
                    updated_lines.append(line)
                    continue

                parts = stripped_line.split()
                if len(parts) >= 2 and parts[0].isdigit() and int(parts[0]) == number:
                    parts[1] = f"{value:.8e}"
                    updated_line = "    " + "    ".join(parts)
                    updated_lines.append(updated_line)
                    inside_block = False
                    continue

            updated_lines.append(line)

        return updated_lines


    template_lines = template.splitlines()

    def recursive_generate(current_params, index):
        if index == len(parameters):
            filename = os.path.join(SLHA_DIRECTORY, f"SLHA_{len(generated_files) + 1}.lha")
            with open(filename, "w") as file:
                file.write("\n".join(current_params))
            generated_files.append(filename)
            return

        param = parameters[index]
        value = param.start
        while value <= param.stop:
            updated_lines = replace_parameter(current_params, param.block, param.number, value)
            recursive_generate(updated_lines, index + 1)
            value += param.step

    recursive_generate(template_lines, 0)
    return generated_files

@router.post("/generate_slha")
async def generate_slha(request: SLHAGenerationRequest):
    if not request.template:
        raise HTTPException(status_code=400, detail="Template content is missing.")
    if not request.parameters:
        raise HTTPException(status_code=400, detail="No parameters provided.")

    try:
        generated_files = generate_slha_files(request.template, request.parameters)
        return {"message": "SLHA files generated successfully!", "files": generated_files}
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"An error occurred: {str(e)}")

@router.get("/list_files")
async def list_slha_files():
    try:
        files = os.listdir(SLHA_DIRECTORY)
        return {"files": files, "count": len(files)}
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Unable to list files: {str(e)}")

@router.delete("/clean_directory")
async def clean_slha_directory():
    try:
        files = os.listdir(SLHA_DIRECTORY)
        for file in files:
            os.remove(os.path.join(SLHA_DIRECTORY, file))
        return {"message": "Directory cleaned successfully."}
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Unable to clean directory: {str(e)}")