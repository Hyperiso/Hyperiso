#include <map>
#include <string>
#include <unordered_map>
#include <fstream>
#include <iostream>
#include <sstream>
#include "ModelModifier.h"
#include "Extractor.h"
#include "Interpreter.h"
#include "SMParamSetter.h"

class NumModelModifier : public ModelModifier {
private:
    std::map<std::string, std::string> paramMap;  // Map pour stocker les paramètres
    std::unordered_map<std::string, double> params;
    bool done = false;  // Indique si la section marquée par "//42" a été atteinte
    bool forceMode = false;  // Mode forcé pour ignorer le "//42"
    int count = 0;  // Compteur de lignes
    std::string wilson;  // Nom spécifique à Wilson
    Extractor extractor;
    Interpreter interpreter;
    SMParamSetter paramSetter; 

public:
    // Constructeur qui extrait et interprète les paramètres dès l'initialisation
    NumModelModifier(const std::string& wilson, bool force = false)
        : wilson(wilson), forceMode(force), paramSetter(params) {
        // Extraire et interpréter les paramètres au niveau du constructeur
        std::string filename = "libs/" + wilson + "_SM/include/params.h";
        auto extractedParams = extractor.extract(filename);
        auto interpretedParams = interpreter.interpret(extractedParams);  // Initialiser la map de paramètres

        // Utilisation du ParamSetter pour initialiser les paramètres
        for (const auto& [name, interpreted] : interpretedParams) {
            paramSetter.setParam(name, interpreted);  // SetParam via SMParamSetter
        }

        // Copie des paramètres dans la map finale
        // params = paramSetter.getParams();
    }

    void addLine(std::ofstream& outputFile, const std::string& currentLine, bool addBefore) override {
        // Si la ligne contient "//42" et que le mode "force" n'est pas activé, marquer comme terminé
        std::cout << "----------------------------------" << std::endl;
        std::cout << currentLine << std::endl;
        if (currentLine.find("//42") != std::string::npos && !forceMode) {
            this->done = true;
        }

        // Si la section est marquée comme terminée et que le mode "force" n'est pas activé, on n'écrit plus
        if (done && !forceMode) {
            outputFile << currentLine << "\n";
            return;
        }

        // Si on trouve un paramètre déjà défini dans le fichier
        if (currentLine.find("param.") != std::string::npos) {
            std::string paramName = extractParamName(currentLine);
            std::string paramValue = extractParamValue(currentLine);

            // Si le paramètre existe dans la map, on le met à jour
            if (params.find(paramName) != params.end()) {
                outputFile << "\tparam." << paramName << " = " << params[paramName] << ";\n";
                params.erase(paramName);  // On supprime le paramètre de la map car il a déjà été traité
            } else {
                // Si le paramètre n’est pas dans la map, on écrit la ligne inchangée
                outputFile << currentLine << "\n";
            }
        } 
        // Si on arrive à "return 0;", on écrit tous les paramètres restants
        else if (currentLine.find("return 0;") != std::string::npos) {
            if (addBefore) {
                outputFile << "\tparam_t param;\n";

                // Ajouter tous les paramètres restants dans la map
                for (const auto& [name, value] : params) {
                    std::string real_name = name;
                    if (name.find("_im") != std::string::npos) {
                        real_name = name.substr(0, name.size() - 3);
                        outputFile << "\tparam." << real_name << " = {"
                                   << params[real_name + "_re"] << ", " << value << "};\n";
                    } else if (name.find("_re") != std::string::npos) {
                        continue;
                    } else {
                        outputFile << "\tparam." << real_name << " = " << value << ";\n";
                    }
                }

                // Ajouter le code spécifique lié à Wilson
                outputFile << "\tauto out = std::ofstream(\"" + wilson + "_SM.txt\");\n";
                outputFile << "\tout << " + wilson + "(param).real() << \" \" << "
                           + wilson + "(param).imag() << std::endl;\n";
                outputFile << "\tout.close();\n";
            }

            // Écrire "return 0;" après avoir ajouté les paramètres restants
            outputFile << currentLine << "\n";
        } 
        // Gestion spéciale pour inclure <fstream> lors de l'utilisation de "using namespace"
        else if (currentLine.find("using namespace") != std::string::npos) {
            outputFile << "#include <fstream>\n";
            outputFile << currentLine << "\n";
        } 
        // Si aucune condition spéciale n'est remplie, on écrit la ligne inchangée
        else {
            outputFile << currentLine << "\n";
        }

        count++;
    }

    void modifyLine(std::string& line) override {}
private:
    // Méthode pour extraire le nom du paramètre à partir de la ligne
    std::string extractParamName(const std::string& line) {
        size_t startPos = line.find("param.") + 6;
        size_t endPos = line.find(" =");
        return line.substr(startPos, endPos - startPos);
    }

    // Méthode pour extraire la valeur du paramètre à partir de la ligne
    std::string extractParamValue(const std::string& line) {
        size_t startPos = line.find("= ") + 2;
        size_t endPos = line.find(";");
        return line.substr(startPos, endPos - startPos);
    }
};