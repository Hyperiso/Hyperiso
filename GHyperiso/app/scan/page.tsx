import { useEffect, useState } from "react";
import Plot from "react-plotly.js";
import { Card, CardContent } from "@/components/ui/card";
import axios from "axios";

export default function ScanPage() {
  const [scanData, setScanData] = useState<{ x: number[]; y: number[] } | null>(null);

  useEffect(() => {
    axios.get("http://localhost:8000/api/scan/mock").then((res) => {
      setScanData(res.data);
    });
  }, []);

  return (
    <div className="p-6">
      <h1 className="text-2xl font-bold mb-4">Exploration de scan</h1>
      <Card>
        <CardContent className="p-4">
          {scanData ? (
            <Plot
              data={[
                {
                  x: scanData.x,
                  y: scanData.y,
                  type: "scatter",
                  mode: "lines+markers",
                },
              ]}
              layout={{ title: "Résultat de scan (mock)", xaxis: { title: "param" }, yaxis: { title: "observable" } }}
              style={{ width: "100%", height: "500px" }}
            />
          ) : (
            <p>Chargement des données...</p>
          )}
        </CardContent>
      </Card>
    </div>
  );
}
