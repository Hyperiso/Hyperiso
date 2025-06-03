//42
#include "c6_thdm.h"

#include <fstream>
#include "csv_helper.h"
using namespace c6_thdm;

int main(int argc, char** argv) {

	param_t param;
	std::string ParamFilePath = "/home/theo/hyperiso/Assets/MartyTemp/libs/C6_THDM/bin/paramlist.csv";
	std::ifstream ParamFile(ParamFilePath);
	readParams(ParamFile, param.realParams, param.complexParams);
	double Q_match = 80.379;
	for (int i = 1; i < argc; i++) {
		if (std::string(argv[i]) == "--Q_match" || std::string(argv[i]) == "-Q") {
			Q_match = std::stod(argv[i + 1]);
			i++;
		} else if (std::string(argv[i]) == "--help" || std::string(argv[i]) == "-h") {
			std::cout << "Options availables :" << std::endl;
			std::cout << "--Q_match/-Q : Value of Q_match (default 80.379)" << std::endl;
			std::cout << "--help/-h : Affiche ce message." << std::endl;
			return 0;
		}
	}

	std::string path = "/home/theo/hyperiso/Assets/MartyTemp/THDM_wilson.csv";
	setMu(Q_match);
	writeWilsonCoefficients("C6", C6(param), Q_match, path);
    return 0;
}
