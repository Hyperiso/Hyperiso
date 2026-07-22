#include <cassert>
#include <filesystem>
#include <fstream>
#include <iostream>

#include "CinematicExtractor.h"

namespace fs = std::filesystem;

int main() {
    std::cout << "== CinematicExtractor UNIT ==\n";
    const fs::path source = fs::temp_directory_path() / "hyperiso_cinematic_extractor_unit.cpp";
    {
        std::ofstream out(source);
        out << "auto wil = model.computeWilsonCoefficients(\n"
            << "  mty::Order::TreeLevel,\n"
            << "  {Incoming(\"b\"), Outgoing(\"s\"), Outgoing(\"mu\"), "
            << "Outgoing(AntiPart(\"mu\"))}, opts);\n";
    }

    CinematicExtractor extractor;
    const auto process = extractor.extract_process(source.string());
    assert(process.incoming_count() == 1);
    assert(process.outgoing_count() == 3);
    assert(process.incoming.at(0) == "b");
    assert(process.outgoing.at(0) == "s");
    assert(process.outgoing.at(1) == "mu");
    assert(process.outgoing.at(2) == "mu");

    const fs::path repeated_source =
        fs::temp_directory_path() / "hyperiso_cinematic_extractor_repeated_unit.cpp";
    {
        std::ofstream out(repeated_source);
        out << "std::vector<Insertion> insertions() {\n"
            << "  return {Incoming(\"b\"), Outgoing(\"s\"), Outgoing(\"mu\"), "
            << "Outgoing(AntiPart(\"mu\"))};\n}\n"
            << "auto first = model.computeWilsonCoefficients(order, insertions(), opts);\n"
            << "auto second = model.computeWilsonCoefficients(order, insertions(), opts);\n";
    }
    const auto repeated = extractor.extract_process(repeated_source.string());
    assert(repeated.incoming_count() == 1);
    assert(repeated.outgoing_count() == 3);
    assert(repeated.incoming.at(0) == "b");
    assert(repeated.outgoing.at(0) == "s");

    std::cout << "CinematicExtractor wrapped-antiparticle and duplicate-list tests OK\n";
    return 0;
}
