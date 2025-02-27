#include "Parameters.h"
#include "MemoryManager.h"

int main() {
    Logger::getInstance()->setLevel(Logger::LogLevel::INFO);

    Config config;

    config.model = Model::SM;

    MemoryManager::GetInstance()->init("lha/testInput.flha", config);

    Parameters::GetInstance();
    return 0;
}