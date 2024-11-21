#include "Parameters.h"

int main() {
    MemoryManager::GetInstance()->init("Test/InputFiles/testinput_thdm.lha", {0,4});
    Parameters::GetInstance(4);
    return 0;
}