"use client";
import { useState } from "react";
import { Button } from "@/components/ui/button";
import { Input } from "@/components/ui/input";
import { Card, CardContent } from "@/components/ui/card";
import axios from "axios";

export default function LhaAnalysePage() {
  const [filePath, setFilePath] = useState("lha/camilia.flha");
  const [output, setOutput] = useState<string>("");

  const handleAnalyse = async () => {
    try {
      const res = await axios.post("http://localhost:8000/api/lha/analyse", {
        path: filePath,
      });
      setOutput(JSON.stringify(res.data, null, 2));
    } catch (err) {
      console.error(err);
      setOutput("Erreur lors de l'analyse du fichier LHA");
    }
  };

  return (
    <div className="p-6">
      <h1 className="text-2xl font-bold mb-4">Analyse de fichier LHA</h1>
      <Card className="mb-4">
        <CardContent className="p-4 space-y-4">
          <Input
            value={filePath}
            onChange={(e) => setFilePath(e.target.value)}
            placeholder="Chemin vers fichier .flha"
          />
          <Button onClick={handleAnalyse}>Analyser</Button>
        </CardContent>
      </Card>
      {output && (
        <pre className="bg-black text-white p-4 rounded overflow-auto text-sm">
          {output}
        </pre>
      )}
    </div>
  );
}