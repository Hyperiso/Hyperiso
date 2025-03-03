#include "Parameters.h"
#include "MemoryManager.h"

int main() {
    Logger::getInstance()->setLevel(Logger::LogLevel::INFO);

    Config config;

    config.model = Model::SM;

    MemoryManager::GetInstance()->init("lha/testInput.flha", config);

    Parameters::GetInstance();

    DBMemento memento;
    for (size_t i = 0; i < memento.stack_size(); i++){
        memento.print_snapshot_content(i);
    }

    return 0;
}