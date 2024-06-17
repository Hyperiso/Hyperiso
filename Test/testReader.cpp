#include "lha_reader.h"

#include <iostream>
#include <cassert>

int main() {
    Logger* logger = Logger::getInstance();
    logger->setLogFile("app.log");
    logger->setLevel(Logger::LogLevel::DEBUG);

    LOG_WARN("Warning message with multiple", "arguments", 123, 45.6);
    LhaReader reader("../testInput.flha");
    reader.addBlockType("testadd", 2, 1);
    reader.readAll();

    // Check that TESTSKIP has indeed been skipped
    assert(reader.getBlockCount() == 13);
    
    // Check complex ID and global scale parsing
    auto wilsonC1 = static_cast<LhaElement<double>*>(reader.getBlock("FWCOEF")->get("03040405|6161|00|2"));
    assert(wilsonC1->getValue() == 0.);
    assert(wilsonC1->getScale() == 1.60846e+02);

    // Check custom block addition and case-insentitive search
    auto t_12_L = static_cast<LhaElement<double>*>(reader.getBlock("Testadd")->get("1"));
    assert(t_12_L->getValue() == 2.4579); 

    // Check minus sign detection
    auto delta_CP = static_cast<LhaElement<double>*>(reader.getBlock("upmnsin")->get("4"));
    assert(delta_CP->getValue() == -2.84); 

    return 0;
}