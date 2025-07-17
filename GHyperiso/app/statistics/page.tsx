import { useEffect, useState } from "react";
import { Card, CardContent } from "@/components/ui/card";
import axios from "axios";

export default function StatisticsPage() {
  const [chi2, setChi2] = useState<number | null>(null);
  const [dof, setDof] = useState<number | null>(null);

  useEffect(() => {
    axios.get("http://localhost:8000/api/statistics/chi2").then((res) => {
      setChi2(res.data.chi2);
      setDof(res.data.dof);
    });
  }, []);

  return (
    <div className="p-6">
      <h1 className="text-2xl font-bold mb-4">Statistiques (χ²)</h1>
      <Card>
        <CardContent className="p-4 text-lg">
          {chi2 !== null && dof !== null ? (
            <>
              <p>χ² = {chi2.toFixed(3)}</p>
              <p>degrés de liberté (dof) = {dof}</p>
              <p>χ²/dof = {(chi2 / dof).toFixed(3)}</p>
            </>
          ) : (
            <p>Chargement...</p>
          )}
        </CardContent>
      </Card>
    </div>
  );
}
