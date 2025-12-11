#include "FileWriter.h"

FileWriter::FileWriter(const std::string& wilson, const std::string& model) :
    wilson(wilson), model(model) {}

void FileWriter::add_output_writer(std::ofstream& outputFile) {
    outputFile << "\tstd::string path = \"" << FileNameManager::getInstance(this->wilson, this->model)->getCsvWilsonFileName() <<"\";\n";
    outputFile << "\tsetMu(Q_match);\n";
    outputFile << "\twriteWilsonCoefficients(\"" + wilson + "\", " + wilson + "(param), Q_match, path);\n";

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