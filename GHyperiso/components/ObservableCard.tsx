// GHyperiso/components/ObservableCard.tsx
interface ObservableCardProps {
  name: string;
  value: number;
  uncertainty: number;
}

export default function ObservableCard({ name, value, uncertainty }: ObservableCardProps) {
  return (
    <div className="border rounded-xl shadow-md p-4 bg-white">
      <h2 className="font-semibold text-lg mb-2">{name}</h2>
      <p className="text-sm">
        Valeur : <strong>{value}</strong>
      </p>
      <p className="text-sm">
        Incertitude : ±<strong>{uncertainty}</strong>
      </p>
    </div>
  );
}