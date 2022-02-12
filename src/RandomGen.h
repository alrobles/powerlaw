#pragma once
#include <numeric>
#include <random>

/**
* @enum RandomAlgorithm
* @brief Types of random generators
*/
enum class RandomAlgorithm
{
    LCG, MT19937, RANLUX24, RANLUX48
};

/**
* @class RandomGen
* @brief Random number generator
*/
class RandomGen
{
    static RandomAlgorithm m_ra;
    static std::mt19937 mt;
    static std::ranlux24 rl24;
    static std::ranlux48 rl48;
public:
    static void SetAlgorithm(RandomAlgorithm ra);
    static void Seed(int seed = -1);
    static int GetInt(int max);
    static double GetUniform01();
};