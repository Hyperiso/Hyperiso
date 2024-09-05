#include "MartyInterface.h"

int main() {

    MartyInterface MartyInterface;
    MartyInterface.generate("C7", "SM");
    MartyInterface.compile("C7", "SM");
    return 0;
}