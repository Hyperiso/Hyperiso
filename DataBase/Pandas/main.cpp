#include "CSVReader.h"
#include "DataFrame.h"
#include <iostream>

int main() {
    CSVReader pd;

    CSVOptions options;
    options.hasIndex = true;
    options.indexType = typeid(std::string);
    options["truc1"] = typeid(int);
    options["truc2"] = typeid(double);
    options["truc3"] = typeid(double);


    DataFrame df = pd.read_csv("data.csv", options);

    df.head();

    df.describe();
    std::cout << df.getColumnNames() << std::endl;;

    df.to_csv("output.csv");


    df.print();
    
    for (auto& i : df.getColumnNames()) {
        std::cout << i << "\n";
    }
    std::cout << df.getIndex<std::string>().iat(0) << std::endl;;
    std::cout << "\n";
    
    try {
        int value = df.iat<int>(0, "truc1");
        std::cout << "Value at (0, truc1): " << value << std::endl;
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
    }

    


    return 0;
}