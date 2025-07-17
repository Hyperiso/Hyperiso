// GHyperiso/lib/api.ts
import axios from "axios";

const api = axios.create({
  baseURL: "http://localhost:8000/api",
});

export default api;

// Exemple d’usage :
// const res = await api.get("/observables")
