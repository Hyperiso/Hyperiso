import { useEffect, useState } from "react";
import { Card, CardContent } from "@/components/ui/card";
import { Table, TableBody, TableCell, TableHead, TableHeader, TableRow } from "@/components/ui/table";
import axios from "axios";

export default function ObservablesPage() {
  const [obs, setObs] = useState<any[]>([]);

  useEffect(() => {
    axios.get("http://localhost:8000/api/observables").then((res) => {
      setObs(res.data);
    });
  }, []);

  return (
    <div className="p-6">
      <h1 className="text-2xl font-bold mb-4">Observables</h1>
      <Card>
        <CardContent className="p-4">
          <Table>
            <TableHeader>
              <TableRow>
                <TableHead>Nom</TableHead>
                <TableHead>Valeur</TableHead>
                <TableHead>Incertitude</TableHead>
              </TableRow>
            </TableHeader>
            <TableBody>
              {obs.map((o, i) => (
                <TableRow key={i}>
                  <TableCell>{o.name}</TableCell>
                  <TableCell>{o.value}</TableCell>
                  <TableCell>{o.uncertainty}</TableCell>
                </TableRow>
              ))}
            </TableBody>
          </Table>
        </CardContent>
      </Card>
    </div>
  );
}