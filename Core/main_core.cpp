#include "Parameters.h"

int main() {
    MemoryManager::GetInstance()->init("Test/InputFiles/testinput_thdm.lha", Model::CUSTOM);
    Parameters::GetInstance(ParameterType::CUSTOM);
    return 0;
}