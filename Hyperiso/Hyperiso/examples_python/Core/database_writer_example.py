"""Export the initialized Core database to JSON, YAML and FLHA files."""

import argparse
from pathlib import Path

import pyhyperiso
from pyhyperiso.Common import Model, ParameterType
from pyhyperiso.Core import DatabaseWriter, HyperisoConfig, HyperisoMaster
from pyhyperiso.core.Common.ParamId import ParamId


def main() -> None:
    """Initialize HyperIso and demonstrate full and selective database exports."""
    default_lha = Path(pyhyperiso.__file__).resolve().parent / "assets" / "lha" / "si_input.flha"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lha", type=Path, default=default_lha, help="Input FLHA/LHA file")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("database_exports_python"),
        help="Directory receiving the exported database files",
    )
    args = parser.parse_args()

    hyperiso = HyperisoMaster()
    hyperiso.init(
        lha_file=str(args.lha),
        config=HyperisoConfig(model=Model.SM),
    )

    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    writer = DatabaseWriter()

    # Export the complete initialized Core database.
    writer.write(output_dir / "database.json")
    writer.write(output_dir / "database.yaml")
    writer.write(output_dir / "database.flha")

    # Export a block subset.
    writer.write_blocks(
        output_dir / "sm_inputs.yaml",
        ["SMINPUTS", "MASS"],
    )

    # Export selected block/id entries.
    writer.write_parameters(
        output_dir / "selected_parameters.json",
        [
            ParamId(type=ParameterType.SM, block="MASS", code=25),
            ParamId(type=ParameterType.SM, block="SMINPUTS", code=3),
        ],
    )

    print(f"Database exports written to {output_dir.resolve()}")


if __name__ == "__main__":
    main()
