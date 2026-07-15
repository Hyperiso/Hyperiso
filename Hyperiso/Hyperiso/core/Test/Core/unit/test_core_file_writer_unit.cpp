#include <cassert>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <memory>
#include <stdexcept>
#include <string>

#include "Block.h"
#include "BlockAccessor.h"
#include "FileWriter.h"
#include "ParamBlockLoader.h"
#include "Parameter.h"

namespace {

std::string read_text(const std::filesystem::path& path) {
    std::ifstream input(path);
    return std::string(
        std::istreambuf_iterator<char>(input),
        std::istreambuf_iterator<char>()
    );
}

std::shared_ptr<BlockAccessor> make_accessor() {
    auto accessor = std::make_shared<BlockAccessor>();

    auto mass = std::make_shared<Block>();
    mass->blockname = BlockName("MASS");
    mass->store(
        LhaID(25),
        std::make_shared<Parameter>(
            ParamId(BlockName("MASS"), LhaID(25)),
            scalar_t(125.25, 0.125),
            0.1,
            0.2
        )
    );
    accessor->emplace(mass->blockname, mass);

    auto inputs = std::make_shared<Block>();
    inputs->blockname = BlockName("SMINPUTS");
    inputs->store(
        LhaID(3),
        std::make_shared<Parameter>(
            ParamId(BlockName("SMINPUTS"), LhaID(3)),
            scalar_t(0.1181),
            0.0,
            0.0
        )
    );
    accessor->emplace(inputs->blockname, inputs);

    auto custom = std::make_shared<Block>();
    custom->blockname = BlockName("CUSTOMBLOCK");
    custom->store(
        LhaID(1, 2),
        std::make_shared<Parameter>(
            ParamId(BlockName("CUSTOMBLOCK"), LhaID(1, 2)),
            scalar_t(4.5),
            0.0,
            0.0
        )
    );
    accessor->emplace(custom->blockname, custom);

    return accessor;
}

} // namespace

int main() {
    const auto output_dir =
        std::filesystem::temp_directory_path() / "hyperiso_core_file_writer_unit";
    std::filesystem::remove_all(output_dir);

    const auto accessor = make_accessor();
    FileWriter writer;

    const auto json_path = output_dir / "database.json";
    const auto yaml_path = output_dir / "database.yaml";
    const auto lha_path = output_dir / "database.slha";

    writer.write(json_path.string(), accessor);
    writer.write(yaml_path.string(), accessor);
    writer.write(lha_path.string(), accessor);

    assert(std::filesystem::is_regular_file(json_path));
    assert(std::filesystem::is_regular_file(yaml_path));
    assert(std::filesystem::is_regular_file(lha_path));

    const std::string json = read_text(json_path);
    const std::string yaml = read_text(yaml_path);
    const std::string lha = read_text(lha_path);

    assert(json.find("\"MASS\"") != std::string::npos);
    assert(json.find("\"imaginary_value\"") != std::string::npos);
    assert(yaml.find("imaginary_value") != std::string::npos);
    assert(lha.find("BLOCK MASS") != std::string::npos);
    assert(lha.find("BLOCK SMINPUTS") != std::string::npos);
    assert(lha.find("Block CUSTOMBLOCK") != std::string::npos);
    assert(lha.find("Block IMMASS") != std::string::npos);
    assert(lha.find("1.25000000e-01") != std::string::npos);

    // JSON/YAML exports preserve complex values through the database loader.
    auto loaded = std::make_shared<BlockAccessor>();
    ParamBlockLoader loader;
    loader.load(loaded, json_path);
    const scalar_t loaded_mass = loaded->getValue(BlockName("MASS"), LhaID(25));
    assert(std::abs(loaded_mass.real() - 125.25) < 1e-12);
    assert(std::abs(loaded_mass.imag() - 0.125) < 1e-12);

    const auto blocks_path = output_dir / "mass_only.json";
    writer.write_blocks(blocks_path.string(), {BlockName("MASS")}, accessor);
    const std::string blocks_json = read_text(blocks_path);
    assert(blocks_json.find("\"MASS\"") != std::string::npos);
    assert(blocks_json.find("\"SMINPUTS\"") == std::string::npos);

    const auto params_path = output_dir / "parameter_only.yaml";
    writer.write_parameters(
        params_path.string(),
        {ParamId(BlockName("SMINPUTS"), LhaID(3))},
        accessor
    );
    const std::string params_yaml = read_text(params_path);
    assert(params_yaml.find("SMINPUTS") != std::string::npos);
    assert(params_yaml.find("MASS") == std::string::npos);

    bool unsupported_extension_threw = false;
    try {
        writer.write((output_dir / "database.txt").string(), accessor);
    } catch (const std::invalid_argument&) {
        unsupported_extension_threw = true;
    }
    assert(unsupported_extension_threw);

    std::filesystem::remove_all(output_dir);
    return 0;
}
