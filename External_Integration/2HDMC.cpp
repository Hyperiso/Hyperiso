#ifdef BUILD_WITH_2HDMC
#include "2HDMC.h"
#include "THDM.h"
#include "Constraints.h"
#include "MemoryManager.h"

void TwoHDMCalculator::calculateSpectrum(const std::string& inputFilePath, const std::string& outputFilePath) {

    THDM model;
    std::string root_file = MemoryManager::findNearestHyperisoDirectory();

    bool pset = model.read_LesHouches((inputFilePath).c_str());

    model.write_LesHouches((outputFilePath).c_str(), true, true, true);

}

#endif