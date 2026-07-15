#include <filesystem>
#include <iostream>
#include <vector>

#include "Config.h"
#include "FileWriter.h"
#include "HyperisoMaster.h"
#include "Include.h"
#include "ParamID.h"

int main(int argc, char** argv) {
    HyperisoConfig config;
    config.model = Model::SM;
    config.flags[ExternalFlag::IS_LHA_SPECTRUM] = false;
    config.flags[ExternalFlag::HAS_WILSON_INPUT] = false;
    config.flags[ExternalFlag::HAS_TH_OBSERVABLE_INPUT] = false;
    config.flags[ExternalFlag::HYP_AS_SM_MARTY] = true;

    const std::filesystem::path lha_path =
        argc > 1 ? std::filesystem::path(argv[1])
                 : std::filesystem::path("Assets/lha/si_input.flha");
    const std::filesystem::path output_dir =
        argc > 2 ? std::filesystem::path(argv[2])
                 : std::filesystem::path("database_exports_cpp");

    HyperisoMaster hyperiso;
    hyperiso.init(lha_path.string(), config);

    std::filesystem::create_directories(output_dir);

    FileWriter writer;

    // Export the complete initialized Core database in the three supported
    // serialization families.
    writer.write((output_dir / "database.json").string());
    writer.write((output_dir / "database.yaml").string());
    writer.write((output_dir / "database.flha").string());

    // Export a block subset.
    writer.write_blocks(
        (output_dir / "sm_inputs.yaml").string(),
        {BlockName("SMINPUTS"), BlockName("MASS")}
    );

    // Export selected block/id entries.
    writer.write_parameters(
        (output_dir / "selected_parameters.json").string(),
        {
            ParamId(ParameterType::SM, BlockName("MASS"), LhaID(25)),
            ParamId(ParameterType::SM, BlockName("SMINPUTS"), LhaID(3)),
        }
    );

    std::cout << "Database exports written to " << output_dir << '\n';
    return 0;
}
