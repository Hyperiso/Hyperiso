#include "FileWriter.h"

#include <unordered_set>

namespace {
const std::unordered_set<std::string>& wilsons_without_mudim() {
    static const std::unordered_set<std::string> values = {
        "C10"
    };
    return values;
}
}

FileWriter::FileWriter(const std::string& wilson, const std::string& model, bool bsm_split_generation) :
    wilson(wilson), model(model), bsm_split_generation(bsm_split_generation) {}

bool FileWriter::should_set_mudim() const {
    return wilsons_without_mudim().find(this->wilson) == wilsons_without_mudim().end();
}

void FileWriter::add_output_writer(std::ofstream& outputFile) {
    outputFile << "\tstd::string path = \"" << FileNameManager::getInstance(this->wilson, this->model)->getCsvWilsonFileName() <<"\";\n";

    if (should_set_mudim()) {
        outputFile << "\tsetMu(Q_match);\n";
    } else {
        outputFile << "\t// setMu(Q_match) intentionally skipped for " << this->wilson
                   << " because LoopTools mudim must stay at its default value.\n";
    }

    outputFile << "\twriteWilsonCoefficients(\"" + wilson + "\", " + wilson + "(param), Q_match, path);\n";

    if (bsm_split_generation) {
        outputFile << "\twriteWilsonCoefficients(\"" + wilson + "_SM_COMPONENT\", " + wilson + "_SM(param), Q_match, path);\n";
        outputFile << "\twriteWilsonCoefficients(\"" + wilson + "_TOTAL_COMPONENT\", " + wilson + "_TOT(param), Q_match, path);\n";
    }

}

void FileWriter::add_argpars(std::ofstream& outputFile) {
    
    outputFile << "\tdouble Q_match = 80.379;\n";
    outputFile << "\tfor (int i = 1; i < argc; i++) {\n";
    outputFile << "\t\tif (std::string(argv[i]) == \"--Q_match\" || std::string(argv[i]) == \"-Q\") {\n";
    outputFile << "\t\t\tQ_match = std::stod(argv[i + 1]);\n";
    outputFile << "\t\t\ti++;\n";
    outputFile << "\t\t} else if (std::string(argv[i]) == \"--help\" || std::string(argv[i]) == \"-h\") {\n";
    outputFile << "\t\t\tstd::cout << \"Options availables :\" << std::endl;\n";
    outputFile << "\t\t\tstd::cout << \"--Q_match/-Q : Value of Q_match (default 80.379)\" << std::endl;\n";
    outputFile << "\t\t\tstd::cout << \"--help/-h : Affiche ce message.\" << std::endl;\n";
    outputFile << "\t\t\treturn 0;\n";
    outputFile << "\t\t}\n";
    outputFile << "\t}\n\n";
    
}

void FileWriter::add_input_reader(std::ofstream& outputFile) {
    outputFile << "\tparam_t param;\n";

    outputFile << "\tstd::string ParamFilePath = \"" << FileNameManager::getInstance(this->wilson, this->model)->getParamFileName() << "\";\n";
    outputFile << "\tstd::ifstream ParamFile(ParamFilePath);" << "\n";
    outputFile << "\treadParams(ParamFile, param.realParams, param.complexParams);" << "\n";


}