#include "HyperisoInterface.h"

void ProcessFile(HyperisoInterface& instance, const std::string& lhaFile) {
    MemoryManager* mm = instance.GetMemoryManager();
    WilsonInterface* wi = instance.GetWilsonInterface();

    mm->init(lhaFile, {0});
    wi->AddWilsonGroup(WilsonGroups::BCoefficients);
    std::cout << "Processing file: " << lhaFile << std::endl;
    wi->setQMatch(WilsonGroups::BCoefficients, 81.);
    wi->setMatchingCoefficient(WilsonGroups::BCoefficients, CoefficientOrder::LO);

    std::cout << "C7 in SM : " << wi->getMatchingCoefficient(WilsonGroups::BCoefficients, WilsonCoefficientList::C7, CoefficientOrder::LO) << std::endl;
}

int main() {
    std::vector<std::string> lhaFiles = {"file1.slha", "file2.slha", "file3.slha"};

    HyperisoInterface& hyperiso = HyperisoInterface::GetInstance();

    hyperiso.RunInParallel(
        lhaFiles,
        ProcessFile, 
        "SM",
        {0},
        4
    );

}