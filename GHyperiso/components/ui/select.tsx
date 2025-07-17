// GHyperiso/components/ui/select.tsx
"use client";

import * as React from "react";
import { cn } from "@/lib/utils";

export function Select({ children, onValueChange, defaultValue }: {
  children: React.ReactNode;
  onValueChange?: (val: string) => void;
  defaultValue?: string;
}) {
  return (
    <select
      className="border rounded-md px-3 py-2 text-sm"
      onChange={(e) => onValueChange?.(e.target.value)}
      defaultValue={defaultValue}
    >
      {children}
    </select>
  );
}

export function SelectTrigger({ children, className }: React.HTMLAttributes<HTMLDivElement>) {
  return <div className={cn("", className)}>{children}</div>;
}

export function SelectValue({ placeholder }: { placeholder?: string }) {
  return <option disabled>{placeholder}</option>;
}

export function SelectContent({ children }: { children: React.ReactNode }) {
  return <>{children}</>;
}

export function SelectItem({ value, children }: { value: string; children: React.ReactNode }) {
  return <option value={value}>{children}</option>;
}