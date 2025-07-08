"use client";
import { useEffect, useState } from "react";
import { Card, CardContent } from "@/components/ui/card";
import { Table, TableBody, TableCell, TableHead, TableHeader, TableRow } from "@/components/ui/table";
import axios from "axios";

export default function ParametersPage() {
  const [params, setParams] = useState<any[]>([]);

  useEffect(() => {
    axios.get("http://localhost:8000/api/parameters/all").then((res) => {
      setParams(res.data);
    });
  }, []);

  return (
    <div className="p-6">
      <h1 className="text-2xl font-bold mb-4">Paramètres LHA chargés</h1>
      <Card>
        <CardContent className="p-4 overflow-auto">
          <Table>
            <TableHeader>
              <TableRow>
                <TableHead>Block</TableHead>
                <TableHead>Code</TableHead>
                <TableHead>Valeur</TableHead>
              </TableRow>
            </TableHeader>
            <TableBody>
              {params.map((p, i) => (
                <TableRow key={i}>
                  <TableCell>{p.block}</TableCell>
                  <TableCell>{p.code}</TableCell>
                  <TableCell>{p.value}</TableCell>
                </TableRow>
              ))}
            </TableBody>
          </Table>
        </CardContent>
      </Card>
    </div>
  );
}
