#include "TemplateManager.h"

#include <cassert>
#include <filesystem>
#include <fstream>
#include <iostream>

namespace fs = std::filesystem;

class TemplateManagerProbe final : public TemplateManagerBase {
public:
    explicit TemplateManagerProbe(const std::string& templates_dir)
        : TemplateManagerBase(templates_dir) {}

    using TemplateManagerBase::already_generated;

private:
    void generateTemplateImpl(const std::string&, const std::string&) override {}
};

static void write_file(const fs::path& path, const std::string& content) {
    fs::create_directories(path.parent_path());
    std::ofstream output(path);
    output << content;
}

int main() {
    std::cout << "== TemplateManager UNIT ==\n";

    const fs::path root = fs::temp_directory_path() / "template_manager_unit";
    fs::remove_all(root);
    fs::create_directories(root);

    TemplateManagerProbe probe(root.string());

    const fs::path numeric = root / "numeric.cpp";
    write_file(numeric, "//42\nint main() {}\n");
    assert(probe.already_generated(numeric.string()));

    const fs::path analytical = root / "analytical.cpp";
    write_file(analytical, "#include <iostream>\n//42\nint main() {}\n");
    assert(probe.already_generated(analytical.string()));

    const fs::path missing_marker = root / "plain.cpp";
    write_file(missing_marker, "#include <iostream>\nint main() {}\n");
    assert(!probe.already_generated(missing_marker.string()));
    assert(!probe.already_generated((root / "missing.cpp").string()));

    fs::remove_all(root);
    std::cout << "UNIT OK\n";
    return 0;
}
