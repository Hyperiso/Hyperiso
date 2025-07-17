import "./globals.css";
import Link from "next/link";
import { ReactNode } from "react";

export const metadata = {
  title: "GHyperiso",
  description: "Interface graphique pour Hyperiso",
};

export default function RootLayout({ children }: { children: ReactNode }) {
  return (
    <html lang="fr">
      <body className="min-h-screen bg-gray-100 text-gray-900">
        <header className="bg-white shadow p-4">
          <nav className="flex gap-4 text-sm font-medium">
            <Link href="/">Accueil</Link>
            <Link href="/lha-analyse">LHA Analyse</Link>
            <Link href="/parameters">Paramètres</Link>
            <Link href="/wilson">Wilson</Link>
            <Link href="/observables">Observables</Link>
            <Link href="/scan">Scan</Link>
            <Link href="/statistics">Statistiques</Link>
          </nav>
        </header>
        <main className="max-w-7xl mx-auto p-6">{children}</main>
      </body>
    </html>
  );
}