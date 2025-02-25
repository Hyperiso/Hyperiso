#include "Parameters.h"

int main() {
    MemoryManager::GetInstance()->init("lha/testInput.flha", Model::SM);
    return 0;
}