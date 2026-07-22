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

    std::cout << "CinematicExtractor wrapped-antiparticle test OK\n";
    return 0;
}
