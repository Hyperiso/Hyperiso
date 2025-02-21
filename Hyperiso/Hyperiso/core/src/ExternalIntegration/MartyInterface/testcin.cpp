#include "CinematicExtractor.h"
#include <iostream>

int main() {
    CinematicExtractor cinext;
    std::cout << "Number of incoming particle : " << cinext.extract("build/generated_C7.cpp").first << std::endl;
    std::cout << "Number of outgoing particle : " << cinext.extract("build/generated_C7.cpp").second << std::endl;
    return 0;
}