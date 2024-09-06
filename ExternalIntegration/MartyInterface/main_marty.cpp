#include "MartyInterface.h"

int main() {

    MartyInterface MartyInterface;
    // MartyInterface.generate("C7", "SM");
    // MartyInterface.compile_run("C7", "SM");
    MartyInterface.generate_numlib("C7", "SM");
    MartyInterface.compile_run_libs("C7", "SM");
    return 0;
}