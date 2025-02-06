/**
 * @file CinematicExtractor.h
 * @brief Definition of the CinematicExtractor class and the utility function countMatchInRegex.
 * 
 * Here we defines a class that extracts the cinematic of a given process in Marty (for wilson calculation).
 * It also includes a utility function to count occurrences of a regex pattern in a string.
 */

#ifndef CINEMATIC_EXTRACTOR_H
#define CINEMATIC_EXTRACTOR_H

#include <string>
#include <utility>
#include <regex>
#include <fstream>

/**
 * @brief Counts the number of occurrences of a regex pattern in a given string.
 * 
 * @param s The string in which to search.
 * @param re The regex pattern as a string.
 * @return The number of matches found in the string.
 */
int countMatchInRegex(std::string s, std::string re);

/**
 * @class CinematicExtractor
 * @brief A class to extract cinematic information of a given process
 * 
 * This class analyzes a text file and counts the occurrences of the keywords "Incoming" and "Outgoing".
 */
class CinematicExtractor {
public:

    /**
     * @brief Extracts and counts occurrences of the keywords "Incoming" and "Outgoing" in a file.
     * 
     * @param filename The path to the file to analyze.
     * @return A pair of integers containing the number of occurrences of "Incoming" (first element)
     * and "Outgoing" (second element).
     */
    std::pair<int, int> extract(const std::string& filename);
};
#endif // CINEMATIC_EXTRACTOR_H
