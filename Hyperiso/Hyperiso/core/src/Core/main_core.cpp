#include "Parameters.h"

int main() {
    Logger::getInstance()->setLevel(Logger::LogLevel::INFO);
    MemoryManager::GetInstance()->init("lha/testInput.flha", Model::SM);
    return 0;
}