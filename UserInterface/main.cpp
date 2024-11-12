#include <iostream>
#include <string>
#include "WilsonHandler.h"

void print_usage() {
    std::cout << "Usage: ./main <command> [options]\n"
              << "\nCommands:\n"
              << "  wilson      : Manage Wilson coefficients\n"
              << "  other_cmd   : Other commands for the program\n"
              << "\nUse './main <command> --help' for more information on a specific command.\n";
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        print_usage();
        return 1;
    }

    std::string command = argv[1];

    if (command == "wilson") {
        return handleWilsonOptions(argc - 1, argv + 1);
    } else {
        std::cerr << "Unknown command: " << command << std::endl;
        print_usage();
        return 1;
    }

    return 0;
}
