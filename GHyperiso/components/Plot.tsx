// GHyperiso/components/Plot.tsx
import Plot from "react-plotly.js";

interface PlotProps {
  x: number[];
  y: number[];
  title?: string;
  xLabel?: string;
  yLabel?: string;
}

export default function SimplePlot({ x, y, title, xLabel, yLabel }: PlotProps) {
  return (
    <Plot
      data={[
        {
          x,
          y,
          type: "scatter",
          mode: "lines+markers",
          marker: { size: 6 },
        },
      ]}
      layout={{
        title,
        xaxis: { title: xLabel || "x" },
        yaxis: { title: yLabel || "y" },
      }}
      style={{ width: "100%", height: "500px" }}
    />
  );
}