#include "json_parser.h"

std::string trim(const std::string& str) {
    if (str.empty()) {
        return str;
    }
    size_t first = str.find_first_not_of(' ');
    if (first == std::string::npos) {
        return ""; // la chaîne ne contient que des espaces
    }
    size_t last = str.find_last_not_of(' ');
    return str.substr(first, (last - first + 1));
}

std::string remove_quotes(const std::string& str) {
    std::string result = str;
    result.erase(std::remove(result.begin(), result.end(), '"'), result.end());
    result.erase(std::remove(result.begin(), result.end(), ','), result.end());
    return result;
}

std::vector<std::string> split(const std::string& str, char delimiter) {
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(str);
    while (std::getline(tokenStream, token, delimiter)) {
        tokens.push_back(trim(token));
    }
    return tokens;
}


void read_json(const std::string& filename, std::vector<Value>& values, std::vector<Correlation>& correlations) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    std::string line;
    enum class Section { NONE, VALUES, CORRELATIONS };
    Section section = Section::NONE;

    while (std::getline(file, line)) {
        line = trim(line);

        if (line.find("\"values\"") != std::string::npos) {
            section = Section::VALUES;
            continue;
        }

        if (line.find("\"correlations\"") != std::string::npos) {
            section = Section::CORRELATIONS;
            continue;
        }

        if (line.find("{") != std::string::npos || line.find("}") != std::string::npos || line.empty()) {
            continue;
        }

        if (section == Section::VALUES) {
            if (line.find("name") != std::string::npos) {
                Value val;
                std::vector<std::string> parts;
                
                // Lire le champ "name"
                parts = split(line, ':');
                if (parts.size() > 1) {
                    val.name = remove_quotes(parts[1]);
                    std::cout << "Name: " << val.name << std::endl;
                } else {
                    std::cerr << "Error parsing name in values section" << std::endl;
                    continue;
                }

                // Lire le champ "central_value"
                if (std::getline(file, line)) {
                    parts = split(line, ':');
                    if (parts.size() > 1) {
                        val.central_value = std::stod(parts[1]);
                        std::cout << "Central Value: " << val.central_value << std::endl;
                    } else {
                        std::cerr << "Error parsing central_value in values section" << std::endl;
                        continue;
                    }
                }

                // Lire le champ "stat_error"
                if (std::getline(file, line)) {
                    parts = split(line, ':');
                    if (parts.size() > 1) {
                        val.stat_error = std::stod(parts[1]);
                        std::cout << "Stat Error: " << val.stat_error << std::endl;
                    } else {
                        std::cerr << "Error parsing stat_error in values section" << std::endl;
                        continue;
                    }
                }

                // Lire le champ "syst_error"
                if (std::getline(file, line)) {
                    parts = split(line, ':');
                    if (parts.size() > 1) {
                        val.syst_error = std::stod(parts[1]);
                        std::cout << "Syst Error: " << val.syst_error << std::endl;
                    } else {
                        std::cerr << "Error parsing syst_error in values section" << std::endl;
                        continue;
                    }
                }

                values.push_back(val);

                // Lire la ligne de fermeture
                if (std::getline(file, line)) {
                    continue;
                }
            }
        }

        if (section == Section::CORRELATIONS) {
            if (line.find("name1") != std::string::npos) {
                Correlation corr;
                std::vector<std::string> parts;
                
                // Lire le champ "name1"
                parts = split(line, ':');
                if (parts.size() > 1) {
                    corr.name1 = remove_quotes(parts[1]);
                    std::cout << "Name1: " << corr.name1 << std::endl;
                } else {
                    std::cerr << "Error parsing name1 in correlations section" << std::endl;
                    continue;
                }

                // Lire le champ "name2"
                if (std::getline(file, line)) {
                    parts = split(line, ':');
                    if (parts.size() > 1) {
                        corr.name2 = remove_quotes(parts[1]);
                        std::cout << "Name2: " << corr.name2 << std::endl;
                    } else {
                        std::cerr << "Error parsing name2 in correlations section" << std::endl;
                        continue;
                    }
                }

                // Lire le champ "value"
                if (std::getline(file, line)) {
                    parts = split(line, ':');
                    if (parts.size() > 1) {
                        corr.value = std::stod(parts[1]);
                        std::cout << "Value: " << corr.value << std::endl;
                    } else {
                        std::cerr << "Error parsing value in correlations section" << std::endl;
                        continue;
                    }
                }

                correlations.push_back(corr);

                // Lire la ligne de fermeture
                if (std::getline(file, line)) {
                    continue;
                }
            }
        }
    }
}
