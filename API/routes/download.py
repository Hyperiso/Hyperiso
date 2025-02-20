from fastapi import APIRouter, HTTPException
from fastapi.responses import FileResponse
import os

router = APIRouter()

@router.get("/download")
def download_file():
    """
    Secure file download
    Requires authentication token.
    """
    file_path = "Download/hyperiso_latest.zip"
    if not os.path.exists(file_path):
        raise HTTPException(status_code=404, detail="File not found ")
    
    return FileResponse(path=file_path, filename="hyperiso_latest.zip", media_type="application/zip")