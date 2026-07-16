#include "HyperisoInterface.h"

void ProcessFile(HyperisoInterface& instance, const std::string& lhaFile) {
    MemoryManager* mm = instance.GetMemoryManager();
    WilsonInterface* wi = instance.GetWilsonInterface();

    mm->init(lhaFile, {0});
    wi->addWilsonGroup(WGroup::B);
    std::cout << "Processing file: " << lhaFile << std::endl;
    wi->set_matching_scale(WGroup::B, 81.);
    wi->init_group(WGroup::B, CoefficientOrder::LO);

    std::cout << "C7 in SM : " << wi->getMatchingCoefficient(WGroup::B, WilsonCoefficientList::C7, CoefficientOrder::LO) << std::endl;
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