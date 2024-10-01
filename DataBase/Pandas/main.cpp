#include "CSVReader.h"
#include "DataFrame.h"
#include <iostream>

int main() {
    CSVReader pd;

    CSVOptions options;
    options.hasIndex = false;
    options.indexType = typeid(std::string);
    options["truc1"] = typeid(int);
    options["truc2"] = typeid(double);
    options["truc3"] = typeid(double);


    DataFrame df = pd.read_csv("data.csv", options);
    df.print();

    df.head();

    df.describe();

    df.to_csv("output.csv");
    std::cout << df.columns << std::endl;;
    std::cout << df.index << std::endl;
    

    std::cout << df.shape << std::endl;


    return 0;
}