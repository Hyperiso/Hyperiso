// GHyperiso/components/WilsonForm.tsx
import { Input } from "@/components/ui/input";
import { Button } from "@/components/ui/button";
import { Select, SelectTrigger, SelectContent, SelectItem, SelectValue } from "@/components/ui/select";

interface WilsonFormProps {
  group: string;
  coeff: string;
  order: string;
  contribution: string;
  loading: boolean;
  setGroup: (v: string) => void;
  setCoeff: (v: string) => void;
  setOrder: (v: string) => void;
  setContribution: (v: string) => void;
  onRun: () => void;
}

export default function WilsonForm({
  group,
  coeff,
  order,
  contribution,
  loading,
  setGroup,
  setCoeff,
  setOrder,
  setContribution,
  onRun,
}: WilsonFormProps) {
  return (
    <div className="flex flex-wrap gap-4">
      <Select onValueChange={setGroup} defaultValue={group}>
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

      <Select onValueChange={setOrder} defaultValue={order}>
        <SelectTrigger className="w-[100px]">
          <SelectValue placeholder="Ordre" />
        </SelectTrigger>
        <SelectContent>
          <SelectItem value="LO">LO</SelectItem>
          <SelectItem value="NLO">NLO</SelectItem>
          <SelectItem value="NNLO">NNLO</SelectItem>
        </SelectContent>
      </Select>

      <Select onValueChange={setContribution} defaultValue={contribution}>
        <SelectTrigger className="w-[150px]">
          <SelectValue placeholder="Contribution" />
        </SelectTrigger>
        <SelectContent>
          <SelectItem value="TOTAL">Total</SelectItem>
          <SelectItem value="EW">EW</SelectItem>
          <SelectItem value="QCD">QCD</SelectItem>
        </SelectContent>
      </Select>

      <Button onClick={onRun} disabled={loading}>
        {loading ? "Chargement..." : "Lancer le scan mu_W"}
      </Button>
    </div>
  );
}
