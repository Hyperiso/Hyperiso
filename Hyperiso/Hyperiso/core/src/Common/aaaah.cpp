#include <iostream>
#include "General.h"

void printBlock(const BlockName& block) {
    std::cout << "Block: " << block << "\n";
}


int main() {
    BlockName b1{"stone"};
    b1.addAlias("rock").addAlias("pierre");

    BlockName b2{"rock"};

    if (b1 == b2) {
        std::cout << "b1 and b2 are equal!\n";
    }

    if (b1 == "pierre") {
        std::cout << "b1 can also be called 'pierre'\n";
    }

    printBlock("cobblestone");
    printBlock({"grass", "herbe"});
    std::cout << b1 << std::endl;
    return 0;
}