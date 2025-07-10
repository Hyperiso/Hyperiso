# Hyperiso/api/main.py
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from Hyperiso.api import wilson, lha, parameters, scan, observables, statistics

app = FastAPI()

# Autoriser le front (localhost:3000)
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

app.include_router(lha.router, prefix="/api/lha")
app.include_router(wilson.router, prefix="/api/wilson")
app.include_router(parameters.router, prefix="/api/parameters")
app.include_router(scan.router, prefix="/api/scan")
app.include_router(observables.router, prefix="/api/observables")
app.include_router(statistics.router, prefix="/api/statistics")