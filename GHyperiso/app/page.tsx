import Link from "next/link";
import { Card, CardContent } from "@/components/ui/card";

const links = [
  { href: "/lha-analyse", label: "📂 Analyse d’un fichier LHA" },
  { href: "/parameters", label: "📊 Paramètres chargés" },
  { href: "/wilson", label: "🧮 Coefficients de Wilson" },
  { href: "/observables", label: "📈 Observables et erreurs" },
  { href: "/scan", label: "🔁 Exploration de scans" },
  { href: "/statistics", label: "📉 Statistiques (χ²)" },
];

export default function HomePage() {
  return (
    <div className="grid gap-6 md:grid-cols-2">
      {links.map((link) => (
        <Card key={link.href}>
          <CardContent className="p-6">
            <Link
              href={link.href}
              className="text-lg font-semibold hover:underline text-blue-600"
            >
              {link.label}
            </Link>
          </CardContent>
        </Card>
      ))}
    </div>
  );
}
