#include "lha_reader.h"
#include "config.hpp"
#include <iostream>
#include <cassert>

int main() {
    Logger* logger = Logger::getInstance();
    logger->setLevel(Logger::LogLevel::DEBUG);

    std::string root_path = project_root.data();
    LhaParser reader(root_path + "Test/InputFiles/testInput.flha");
    reader.addBlockType("testadd", 2, 1);
    reader.readAll();

    // Check that TESTSKIP has indeed been skipped
    assert(reader.getBlockCount() == 13);
    
    // Check complex ID and global scale parsing
    auto wilsonC1 = static_cast<LhaElement<double>*>(reader.getBlock("FWCOEF")->get({3040405L, 6161, 00, 2}));
    LOG_INFO(wilsonC1);
    assert(wilsonC1->getValue() == 0.);
    assert(wilsonC1->getScale() == 1.60846e+02);

    // Check custom block addition and case-insentitive search
    auto t_12_L = static_cast<LhaElement<double>*>(reader.getBlock("Testadd")->get(1));
    assert(t_12_L->getValue() == 2.4579); 

    // Check minus sign detection
    auto delta_CP = static_cast<LhaElement<double>*>(reader.getBlock("upmnsin")->get(4));
    assert(delta_CP->getValue() == -2.84); 

    return 0;
}