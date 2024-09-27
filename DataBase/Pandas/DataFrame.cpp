#include "DataFrame.h"
#include <iostream>

void DataFrame::print() const {
    // Afficher l'index si disponible
    if (index) {
        try {
            const auto& idx = getIndex<std::string>();  // On assume un index de type int pour la démonstration
            std::cout << "Index\t";
        } catch (...) {
            // Gérer les erreurs de typage
        }
    }

    // Afficher les noms des colonnes
    for (const auto& colName : columns) {
        std::cout << colName << "\t";
    }
    std::cout << std::endl;

    // Afficher les données
    for (size_t i = 0; i < nRows; ++i) {
        // Afficher l'index si disponible
        if (index) {
            try {
                const auto& idx = getIndex<std::string>();  // On assume un index de type int
                std::cout << idx.iat(i) << "\t";
            } catch (...) {
                std::cout << "?\t"; // En cas d'erreur de typage
            }
        }

        // Afficher les valeurs de chaque colonne
        for (const auto& colName : columns) {
            // Supposons que les colonnes sont de types variés
            try {
                // On assume temporairement que les colonnes sont des int
                std::cout << iat<int>(i, colName) << "\t";
            } catch (const std::exception&) {
                std::cout << "NaN\t"; // Valeur manquante ou mauvais type
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
    for (const auto& colName : columns) {
        df.addColumn<int>(colName);
        for (size_t i = 0; i < n; ++i) {
            df.addValueToColumn<int>(colName, iat<int>(i, colName));
        }
    }
    return df;
}

DataFrame DataFrame::tail(size_t n) const {
    DataFrame df;
    n = std::min(n, nRows);
    for (const auto& colName : columns) {
        df.addColumn<int>(colName);
        for (size_t i = nRows - n; i < nRows; ++i) {
            df.addValueToColumn<int>(colName, iat<int>(i, colName));
        }
    }
    return df;
}

void DataFrame::describe() const {
    std::cout << "Describe function is not fully implemented yet." << std::endl;
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