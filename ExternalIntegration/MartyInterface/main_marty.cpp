#include "MartyInterface.h"

int main() {

    MartyInterface MartyInterface;
    MartyInterface.generate("C7", "SM");
    MartyInterface.compile_run("C7", "SM");
    return 0;
}