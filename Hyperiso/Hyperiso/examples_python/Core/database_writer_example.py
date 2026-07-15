"""Export the initialized Core database to JSON, YAML and FLHA files."""

from pathlib import Path

from pyhyperiso.Common import Model, ParameterType
from pyhyperiso.Core import DatabaseWriter, HyperisoConfig, HyperisoMaster
from pyhyperiso.core.Common.ParamId import ParamId


def main() -> None:
    """Initialize HyperIso and demonstrate full and selective database exports."""
    hyperiso = HyperisoMaster()
    hyperiso.init(
        lha_file="lha/si_input.flha",
        config=HyperisoConfig(model=Model.SM),
    )

    output_dir = Path("database_exports_python")
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
