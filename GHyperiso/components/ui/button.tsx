// GHyperiso/components/ui/button.tsx
import * as React from "react";
import { cn } from "@/lib/utils";

export interface ButtonProps extends React.ButtonHTMLAttributes<HTMLButtonElement> {}

const Button = React.forwardRef<HTMLButtonElement, ButtonProps>(({ className, ...props }, ref) => {
  return (
    <button
      ref={ref}
      className={cn("bg-blue-600 text-white px-4 py-2 rounded-md hover:bg-blue-700 transition", className)}
      {...props}
    />
  );
});
Button.displayName = "Button";

export { Button };