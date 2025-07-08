// GHyperiso/components/ParamTable.tsx
import { Table, TableBody, TableCell, TableHead, TableHeader, TableRow } from "@/components/ui/table";

interface ParamEntry {
  block: string;
  code: number;
  value: number;
}

export default function ParamTable({ data }: { data: ParamEntry[] }) {
  return (
    <Table>
      <TableHeader>
        <TableRow>
          <TableHead>Bloc</TableHead>
          <TableHead>Code</TableHead>
          <TableHead>Valeur</TableHead>
        </TableRow>
      </TableHeader>
      <TableBody>
        {data.map((row, i) => (
          <TableRow key={i}>
            <TableCell>{row.block}</TableCell>
            <TableCell>{row.code}</TableCell>
            <TableCell>{row.value}</TableCell>
          </TableRow>
        ))}
      </TableBody>
    </Table>
  );
}
