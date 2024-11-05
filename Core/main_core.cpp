#include "Parameters.h"

int main() {
    MemoryManager::GetInstance("Test/InputFiles/testinput_thdm.lha", {0,4})->init();
    Parameters::GetInstance(4);
    return 0;
}