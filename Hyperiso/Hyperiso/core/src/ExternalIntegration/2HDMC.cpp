#include "2HDMC.h"

#include <filesystem>
#include <stdexcept>

#include "THDM.h"

namespace fs = std::filesystem;

void TwoHDMCalculator::calculateSpectrum(const std::string& inputFilePath, const std::string& outputFilePath)
{
    const fs::path input(inputFilePath);
    const fs::path output(outputFilePath);

    if (!fs::exists(input)) {
        throw std::runtime_error("2HDMC input LHA file does not exist: " + input.string());
    }

    std::error_code ec;
    if (output.has_parent_path()) {
        fs::create_directories(output.parent_path(), ec);
        if (ec) {
            throw std::runtime_error("Cannot create 2HDMC output directory '" + output.parent_path().string() + "': " + ec.message());
        }
    }

    THDM model;
    model.read_LesHouches(input.string().c_str());
    model.write_LesHouches(output.string().c_str(), true, true, true);
}
