#include <iostream>
#include "testQCDParameters.cpp"
#include "testWilson.cpp"

bool testQCDParameters(); 
bool testWilson(); 

int main() {
    bool allTestsPassed = true;

    if (!testQCDParameters()) {
        std::cerr << "QCDParameters test failed." << std::endl;
        allTestsPassed = false;
    }

    if (!testWilson()) {
        std::cerr << "Wilson test failed." << std::endl;
        allTestsPassed = false;
    }


    if (allTestsPassed) {
        std::cout << "All tests passed." << std::endl;
        return 0; // Succès
    } else {
        return 1; // Échec
    }
}