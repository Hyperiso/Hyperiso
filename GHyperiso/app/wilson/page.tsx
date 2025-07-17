"use client";
import { useState } from "react";
import { Button } from "@/components/ui/button";
import { Card, CardContent } from "@/components/ui/card";
import { Input } from "@/components/ui/input";
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from "@/components/ui/select";
import Plot from "react-plotly.js";
import axios from "axios";

export default function WilsonPage() {
  const [group, setGroup] = useState("B");
  const [coeff, setCoeff] = useState("C7");
  const [order, setOrder] = useState("LO");
  const [contribution, setContribution] = useState("TOTAL");
  const [muWScan, setMuWScan] = useState<number[]>([]);
  const [values, setValues] = useState<number[]>([]);
  const [loading, setLoading] = useState(false);

  async function handleScan() {
    setLoading(true);
    try {
      const scanPoints = Array.from({ length: 40 }, (_, i) => 40 + i * 2);
      const results: number[] = [];
      for (const muW of scanPoints) {
        const res = await axios.post("http://localhost:8000/api/wilson/get_M", {
          group,
          coefficient: coeff,
          order,
          contribution,
          mu_W: muW
        });
        results.push(res.data.value); // Suppose Scalar has a .value
      }
      setMuWScan(scanPoints);
      setValues(results);
    } catch (err) {
      console.error(err);
    } finally {
      setLoading(false);
    }
  }

  return (
    <div className="p-6">
      <h1 className="text-2xl font-bold mb-6">Visualisation des coefficients de Wilson</h1>
      <Card className="mb-4">
        <CardContent className="flex flex-wrap gap-4 p-4">
          <Select onValueChange={setGroup} defaultValue="B">
            <SelectTrigger className="w-[150px]">
              <SelectValue placeholder="Groupe" />
            </SelectTrigger>
            <SelectContent>
              <SelectItem value="B">B</SelectItem>
              <SelectItem value="BPrime">B'</SelectItem>
              <SelectItem value="K">K</SelectItem>
            </SelectContent>
          </Select>

          <Input placeholder="Coefficient" value={coeff} onChange={(e) => setCoeff(e.target.value)} className="w-[100px]" />

          <Select onValueChange={setOrder} defaultValue="LO">
            <SelectTrigger className="w-[100px]">
              <SelectValue placeholder="Ordre" />
            </SelectTrigger>
            <SelectContent>
              <SelectItem value="LO">LO</SelectItem>
              <SelectItem value="NLO">NLO</SelectItem>
              <SelectItem value="NNLO">NNLO</SelectItem>
            </SelectContent>
          </Select>

          <Select onValueChange={setContribution} defaultValue="TOTAL">
            <SelectTrigger className="w-[150px]">
              <SelectValue placeholder="Contribution" />
            </SelectTrigger>
            <SelectContent>
              <SelectItem value="TOTAL">Total</SelectItem>
              <SelectItem value="EW">EW</SelectItem>
              <SelectItem value="QCD">QCD</SelectItem>
            </SelectContent>
          </Select>

          <Button onClick={handleScan} disabled={loading}>
            {loading ? "Chargement..." : "Lancer le scan mu_W"}
          </Button>
        </CardContent>
      </Card>

      {muWScan.length > 0 && (
        <Plot
          data={[
            {
              x: muWScan,
              y: values,
              type: "scatter",
              mode: "lines+markers",
              marker: { size: 6 },
            },
          ]}
          layout={{ title: `C${coeff} vs mu_W`, xaxis: { title: "mu_W" }, yaxis: { title: `C${coeff}` } }}
          style={{ width: "100%", height: "500px" }}
        />
      )}
    </div>
  );
}
