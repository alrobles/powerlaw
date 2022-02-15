#pragma once
#include <string>
#include <vector>
#include "../src/VectorUtilities.h"

/// <summary>
/// Parse csv line to a vector.
/// </summary>
/// <typeparam name="T">Cast type of the line element.</typeparam>
/// <param name="in">Csv formatted line.</param>
/// <param name="ignoreNested">If false, doesn't separate by commas in nested structures.</param>
/// <returns>A list with the parsed elements.</returns>
template<typename T> std::vector<T> parse_csv_line(const std::string& input, const bool ignoreNested = true)
{
    const std::vector<char> beginNestSymbols = { '(', '[' };
    const std::vector<char> endNestSymbols = { ')', ']' };
    int nestLevel = 0;

    std::vector<T> out;
    
    int startPos = 0;

    for (unsigned int i = 0; i < input.length(); i++)
    {
        char currentCharacter = input.at(i);
        if (!ignoreNested && VectorUtilities::VectorContainsQ(beginNestSymbols, currentCharacter))
            nestLevel++;

        if (!ignoreNested && VectorUtilities::VectorContainsQ(endNestSymbols, currentCharacter))
            nestLevel--;

        if (nestLevel == 0 && input.at(i) == ',')
        {
            out.push_back(VectorUtilities::ConvertTo<T>(input.substr(startPos, i - startPos)));
            startPos = i + 1;
        }
        if (i + 1 == input.length())
        {
            if (input[input.length() - 1] == '\n' || input[input.length() - 1] == '\r')
                out.push_back(VectorUtilities::ConvertTo<T>(input.substr(startPos, input.length() - startPos - 1)));
            else
                out.push_back(VectorUtilities::ConvertTo<T>(input.substr(startPos, input.length() - startPos)));
        }
    }
    return out;
}