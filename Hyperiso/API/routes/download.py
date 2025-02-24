from fastapi import APIRouter, HTTPException
import os
import requests
from fastapi.responses import StreamingResponse
router = APIRouter()

@router.get("/download")
def download_file():
    """
    Downlaod and send .zip from GitHub.
    """
    github_url = "https://github.com/Hyperiso/Hyperiso/archive/refs/heads/main.zip"
    
    response = requests.get(github_url, stream=True)
    if response.status_code != 200:
        raise HTTPException(status_code=response.status_code, detail="File not found on GitHub")
    
    headers = {"Content-Disposition": "attachment; filename=hyperiso_latest.zip"}
    return StreamingResponse(response.raw, media_type="application/zip", headers=headers)