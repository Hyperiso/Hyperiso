#include "FileWriter.h"

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <stdexcept>
#include <utility>

#include "JsonParser.h"
#include "LhaParser.h"
#include "MemoryManager.h"
#include "ParamBlockWriter.h"
#include "Parameters.h"
#include "YamlParser.h"

namespace {

std::string lower_extension(const std::string& destination) {
    std::string extension = std::filesystem::path(destination).extension().string();
    std::transform(extension.begin(), extension.end(), extension.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return extension;
}

void ensure_parent_directory(const std::string& destination) {
    const auto parent = std::filesystem::path(destination).parent_path();
    if (!parent.empty()) {
        std::error_code error;
        std::filesystem::create_directories(parent, error);
        if (error) {
            throw std::runtime_error(
                "FileWriter: cannot create output directory '" + parent.string() +
                "': " + error.message()
            );
        }
    }
}

std::shared_ptr<BlockAccessor> source_for_parameter(
    const ParamId& parameter_id,
    const std::shared_ptr<BlockAccessor>& fallback
) {
    if (!parameter_id.type.has_value()) {
        return fallback;
    }

    auto parameters = Parameters::GetInstance(parameter_id.type.value());
    if (!parameters || !parameters->get_block_accessor()) {
        throw std::invalid_argument(
            "FileWriter: no parameter repository is available for the requested ParamId"
        );
    }
    return parameters->get_block_accessor();
}

} // namespace

std::shared_ptr<BlockAccessor> FileWriter::resolve_source(std::shared_ptr<BlockAccessor> src) {
    if (src) {
        return src;
    }

    auto* memory_manager = MemoryManager::GetInstance();
    if (!memory_manager || !memory_manager->is_ready()) {
        throw std::logic_error(
            "FileWriter: HyperIso must be initialized before exporting the current database"
        );
    }
    return memory_manager->extract_block_accessor();
}

void FileWriter::write_accessor(
    const std::string& dest,
    const std::shared_ptr<BlockAccessor>& src
) {
    if (!src) {
        throw std::invalid_argument("FileWriter: source BlockAccessor is null");
    }

    ensure_parent_directory(dest);

    std::error_code remove_error;
    std::filesystem::remove(dest, remove_error);
    if (remove_error) {
        throw std::runtime_error(
            "FileWriter: cannot replace output file '" + dest + "': " +
            remove_error.message()
        );
    }

    auto node = std::make_shared<DBNode>();
    ParamBlockWriter().write(node, src);

    const std::string extension = lower_extension(dest);
    if (extension == ".json") {
        JSONParser().writeToFile(dest, node);
    } else if (extension == ".yaml" || extension == ".yml") {
        YAMLParser().writeToFile(dest, node);
    } else if (extension == ".lha" || extension == ".slha" || extension == ".flha") {
        LhaParser().writeToFile(dest, node);
    } else {
        throw std::invalid_argument(
            "FileWriter: unsupported output extension '" + extension +
            "'; expected .json, .yaml, .yml, .lha, .slha or .flha"
        );
    }

    std::error_code error;
    const bool exists = std::filesystem::is_regular_file(dest, error);
    if (error || !exists) {
        throw std::runtime_error("FileWriter: output file was not created: " + dest);
    }
}

void FileWriter::write(const std::string& dest, std::shared_ptr<BlockAccessor> src) {
    write_accessor(dest, resolve_source(std::move(src)));
}

void FileWriter::write_blocks(
    const std::string& dest,
    const std::vector<BlockName>& block_names,
    std::shared_ptr<BlockAccessor> src
) {
    if (block_names.empty()) {
        throw std::invalid_argument("FileWriter: block_names must not be empty");
    }

    const auto source = resolve_source(std::move(src));
    auto selected = std::make_shared<BlockAccessor>();

    for (const auto& block_name : block_names) {
        if (!source->contains(block_name)) {
            throw std::invalid_argument(
                "FileWriter: requested block '" + block_name.canonical() + "' does not exist"
            );
        }
        selected->emplace(block_name, source->at(block_name)->deep_clone_plain());
    }

    write_accessor(dest, selected);
}

void FileWriter::write_parameters(
    const std::string& dest,
    const std::vector<ParamId>& parameter_ids,
    std::shared_ptr<BlockAccessor> src
) {
    if (parameter_ids.empty()) {
        throw std::invalid_argument("FileWriter: parameter_ids must not be empty");
    }

    const auto fallback = resolve_source(std::move(src));
    auto selected = std::make_shared<BlockAccessor>();

    for (const auto& parameter_id : parameter_ids) {
        const auto source = source_for_parameter(parameter_id, fallback);
        if (!source->has_param(parameter_id.block, parameter_id.code)) {
            throw std::invalid_argument(
                "FileWriter: requested parameter '" + parameter_id.block.canonical() +
                "/" + parameter_id.code.to_string() + "' does not exist"
            );
        }

        if (!selected->contains(parameter_id.block)) {
            auto block = std::make_shared<Block>();
            block->blockname = parameter_id.block;

            const auto source_block = source->at(parameter_id.block);
            if (source_block->has_scale()) {
                block->set_scale(source_block->get_scale());
            }
            selected->emplace(parameter_id.block, block);
        }

        auto destination_block = selected->at(parameter_id.block);
        if (!destination_block->contains(parameter_id.code)) {
            const auto source_parameter = source->getParameter(
                parameter_id.block,
                parameter_id.code
            );
            destination_block->store(
                parameter_id.code,
                std::make_shared<Parameter>(*source_parameter)
            );
        }
    }

    write_accessor(dest, selected);
}
