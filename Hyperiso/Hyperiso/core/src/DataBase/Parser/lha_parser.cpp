#include "lha_parser.h"

void LhaParser::tokenize() {
    int cLine = 0;
    int cCol = 0;

    auto rit = std::sregex_iterator(src.begin(), src.end(), analyzer_rx);
    auto rend = std::sregex_iterator();

    while (rit != rend) {
        std::smatch m = *rit;
        size_t group_index = m.size();

        for (size_t idx = 1; idx < m.size(); ++idx) {
            if (m[idx].matched) {
                group_index = idx - 1;
                break;
            }
        }

        auto tokenType = static_cast<TokenType>(group_index);
        auto value = m[group_index + 1].str();

        if (tokenType == TokenType::NEWLINE) {
            ++cLine;
            cCol = 0;
        } else if (tokenType != TokenType::SKIP && value != "") {
            this->tokens.emplace_back(Token{tokenType, value, cLine, cCol});
            ++cCol;
        }

        ++rit;
    }
}