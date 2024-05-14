#include "2HDMC.h"
#include "THDM.h"
#include "Constraints.h"

void TwoHDMCalculator::calculateSpectrum(const std::string& inputFilePath, const std::string& outputFilePath) {

    THDM model;

    bool pset = model.read_LesHouches(inputFilePath);

    model.write_LesHouches(outputFilePath, true, true, true);

}