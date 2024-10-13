#ifndef CINEMATIC_EXTRACTOR_H
#define CINEMATIC_EXTRACTOR_H

#include <string>
#include <utility>
#include <regex>
#include <fstream>

int countMatchInRegex(std::string s, std::string re) {
    std::regex words_regex(re);
    auto words_begin = std::sregex_iterator(s.begin(), s.end(), words_regex);
    auto words_end = std::sregex_iterator();

    return std::distance(words_begin, words_end);
}

class CinematicExtractor {
public:

    std::pair<int, int> extract(const std::string& filename) {
        std::pair<int,int> cinematic{};

        std::ifstream file(filename);
        std::string line;
        while (std::getline(file, line)) {
            cinematic.first += countMatchInRegex(line, "Incoming");
            cinematic.second += countMatchInRegex(line, "Outgoing");
        }
        return cinematic;
    }
};
#endif // CINEMATIC_EXTRACTOR_H
