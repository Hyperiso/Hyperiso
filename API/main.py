from fastapi import FastAPI
from API.routes import download, auth, generation, parameters, wilson, observables

app = FastAPI(title="HyperISO API")

app.include_router(download.router, prefix="/download", tags= ["Download"])
app.include_router(auth.router, prefix="/auth", tags=["Auth"])
app.include_router(generation.router, prefix="/generation", tags=["Generation"])
app.include_router(parameters.router, prefix="/parameters", tags=["Parameters"])
app.include_router(wilson.router, prefix="/wilson", tags=["Wilson Coefficients"])
app.include_router(observables.router, prefix="/observables", tags=["Observables"])