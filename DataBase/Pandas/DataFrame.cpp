#include "DataFrame.h"
#include <iostream>

void DataFrame::print() const {
    // Afficher l'index si disponible
    if (!index.empty()) {
        std::cout << "Index\t";
    }

    // Afficher les noms des colonnes
    for (const auto& colName : columns) {
        std::cout << colName << "\t";
    }
    std::cout << std::endl;

    // Afficher les données
    for (size_t i = 0; i < nRows; ++i) {
        // Afficher l'index si disponible
        if (!index.empty()) {
            std::cout << index[i] << "\t";
        }

        // Afficher les valeurs de chaque colonne
        for (const auto& colName : columns) {
            try {
                // std::cout << "tout va bien22" << std::endl;
                // std::cout << colName << std::endl;
                // std::cout << "name : "<<  this->csvOptions.columnTypes.at(colName).name() << std::endl;
                // const std::type_index& colType = csvOptions.columnTypes.at(colName);
                // std::cout << "tout va bien" << std::endl;
                // for (auto& key : csvOptions.columnTypes) {
                //     std::cout << key.first << " bite ";
                // }
                if (this->csvOptions.columnTypes.at(colName) == typeid(int)) {
                    // std::cout << "tout va bienssdf" << std::endl;
                    std::cout << iat<int>(i, colName) << "\t";
                } else if (this->csvOptions.columnTypes.at(colName) == typeid(double)) {
                    std::cout << iat<double>(i, colName) << "\t";
                } else if (this->csvOptions.columnTypes.at(colName) == typeid(std::string)) {
                    std::cout << iat<std::string>(i, colName) << "\t";
                } else {
                    std::cout << "bug" << std::endl;
                    std::cout << "NaN\t";  // Gestion des types inconnus
                }
            } catch (const std::exception& e) {
                std::cout << e.what() << " " << colName << std::endl;
                std::cout << "ERROR" << std::endl;
                std::cout << "NaN\t";  // Valeur manquante ou mauvais type
            }
        }
        std::cout << std::endl;
    }
}


const std::vector<std::string>& DataFrame::getColumnNames() const {
    return columns;
}

size_t DataFrame::getRowCount() const {
    return nRows;
}

DataFrame DataFrame::head(size_t n) const {
    DataFrame df;
    n = std::min(n, nRows);
    df.index = std::vector<std::string>(index.begin(), index.begin() + n);  // Copie l'index

    for (const auto& colName : columns) {
        const std::type_index& colType = csvOptions.columnTypes.at(colName);

        if (colType == typeid(int)) {
            df.addColumn<int>(colName);
            for (size_t i = 0; i < n; ++i) {
                df.addValueToColumn<int>(colName, iat<int>(i, colName));
            }
        } else if (colType == typeid(double)) {
            df.addColumn<double>(colName);
            for (size_t i = 0; i < n; ++i) {
                df.addValueToColumn<double>(colName, iat<double>(i, colName));
            }
        } else if (colType == typeid(std::string)) {
            df.addColumn<std::string>(colName);
            for (size_t i = 0; i < n; ++i) {
                df.addValueToColumn<std::string>(colName, iat<std::string>(i, colName));
            }
        }
    }
    df._set_csv_options(this->csvOptions);
    std::vector<std::string> ind(&this->index[0], &this->index[n]);
    df.setIndex(ind);
    return df;
}

DataFrame DataFrame::tail(size_t n) const {
    DataFrame df;
    n = std::min(n, nRows);
    df.index = std::vector<std::string>(index.end() - n, index.end());  // Copie l'index

    for (const auto& colName : columns) {
        const std::type_index& colType = csvOptions.columnTypes.at(colName);

        if (colType == typeid(int)) {
            df.addColumn<int>(colName);
            for (size_t i = nRows - n; i < nRows; ++i) {
                df.addValueToColumn<int>(colName, iat<int>(i, colName));
            }
        } else if (colType == typeid(double)) {
            df.addColumn<double>(colName);
            for (size_t i = nRows - n; i < nRows; ++i) {
                df.addValueToColumn<double>(colName, iat<double>(i, colName));
            }
        } else if (colType == typeid(std::string)) {
            df.addColumn<std::string>(colName);
            for (size_t i = nRows - n; i < nRows; ++i) {
                df.addValueToColumn<std::string>(colName, iat<std::string>(i, colName));
            }
        }
    }

    df._set_csv_options(this->csvOptions);
    std::vector<std::string> ind(&this->index[0], &this->index[n]);
    df.setIndex(ind);
    return df;
}

void DataFrame::describe() const {
    std::cout << "Description des colonnes numériques :" << std::endl;
    
    for (const auto& colName : columns) {
        // const std::type_index& colType = csvOptions.columnTypes.at(colName);

        if (csvOptions.columnTypes.at(colName) == typeid(int)) {
            describeColumn<int>(colName);
        } else if (csvOptions.columnTypes.at(colName) == typeid(double)) {
            describeColumn<double>(colName);
        } else {
            std::cout << "La colonne '" << colName << "' n'est pas numérique, donc elle est ignorée." << std::endl;
        }
    }
}

void DataFrame::to_csv(const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file");
    }

    const auto& columns = this->getColumnNames();
    size_t nRows = this->getRowCount();

    // Écrire l'en-tête (index + colonnes)
    if (this->csvOptions.hasIndex) {
        file << "Index";
        for (const auto& colName : columns) {
            file << "," << colName;
        }
        file << "\n";
    } else {
        for (const auto& colName : columns) {
            file << colName;
            if (&colName != &columns.back()) {
                file << ",";
            }
        }
        file << "\n";
    }

    // Écrire les données
    for (size_t i = 0; i < nRows; ++i) {
        if (this->csvOptions.hasIndex) {
            // Écrire l'index
            file << i; // Utiliser l'index simple ici ; vous pouvez personnaliser selon vos besoins.
        }
        
        for (size_t j = 0; j < columns.size(); ++j) {
            const std::string& colName = columns[j];

            // Utiliser le type d'index pour déterminer comment écrire chaque colonne
            // Par exemple : si le type est int, double ou string.
            if (this->csvOptions.columnTypes.at(colName) == typeid(int)) {
                file << this->iat<int>(i, colName);
            } else if (this->csvOptions.columnTypes.at(colName) == typeid(double)) {
                file << this->iat<double>(i, colName);
            } else if (this->csvOptions.columnTypes.at(colName) == typeid(std::string)) {
                file << this->iat<std::string>(i, colName);
            } else {
                throw std::runtime_error("Unsupported column type for column: " + colName);
            }

            if (j != columns.size() - 1 || this->csvOptions.hasIndex) {
                file << ",";
            }
        }
        file << "\n";
    }

    file.close();
    std::cout << "DataFrame successfully written to " << filename << std::endl;
}

const std::vector<std::string>& DataFrame::getIndex() const {
    
    if (index.empty()) {
        throw std::runtime_error("Index is not set");
    }
    return index;
}

void DataFrame::setIndex(const std::vector<std::string>& newIndex) {
    if (newIndex.size() != nRows && nRows > 0) {
        throw std::invalid_argument("Index size does not match number of rows");
    }
    index = newIndex;
}