#include <filesystem>
#include <iostream>
#include <vector>

#include "Config.h"
#include "FileWriter.h"
#include "HyperisoMaster.h"
#include "Include.h"
#include "ParamID.h"

int main() {
    HyperisoConfig config;
    config.model = Model::SM;
    config.flags[ExternalFlag::IS_LHA_SPECTRUM] = false;
    config.flags[ExternalFlag::HAS_WILSON_INPUT] = false;
    config.flags[ExternalFlag::HAS_TH_OBSERVABLE_INPUT] = false;
    config.flags[ExternalFlag::HYP_AS_SM_MARTY] = true;

    HyperisoMaster hyperiso;
    hyperiso.init("lha/si_input.flha", config);

    const std::filesystem::path output_dir = "database_exports_cpp";
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
