#include "RandomGen.h"
#include <chrono>

RandomAlgorithm RandomGen::m_ra = RandomAlgorithm::MT19937;
std::mt19937 RandomGen::mt;
std::ranlux24 RandomGen::rl24;
std::ranlux48 RandomGen::rl48;

void RandomGen::SetAlgorithm(RandomAlgorithm ra)
{
    m_ra = ra;
}
void RandomGen::Seed(int seed)
{
    if (seed == -1)
        seed = static_cast<int>(std::chrono::system_clock::now().time_since_epoch().count());

    switch (m_ra)
    {
        case RandomAlgorithm::MT19937:
            mt.seed(seed);
            break;
        case RandomAlgorithm::LCG:
            srand(seed);
            break;
        case RandomAlgorithm::RANLUX24:
            rl24.seed(seed);
            break;
        case RandomAlgorithm::RANLUX48:
            rl48.seed(seed);
            break;
    };
}
int RandomGen::GetInt(int max)
{
    switch (m_ra)
    {
        case RandomAlgorithm::MT19937:
            return (int) mt() % max;
        case RandomAlgorithm::LCG:
            return rand() % max;
        case RandomAlgorithm::RANLUX24:
            return (int) rl24() % max;
        case RandomAlgorithm::RANLUX48:
            return (int) rl48() % max;
    };
    return 0;
}
double RandomGen::GetUniform01()
{
    switch (m_ra)
    {
        case RandomAlgorithm::MT19937:
            return (double)mt() / (double) mt.max();
        case RandomAlgorithm::LCG:
            return (double)rand() / (double)RAND_MAX;
        case RandomAlgorithm::RANLUX24:
            return (double)rl24() / (double)rl24.max();
        case RandomAlgorithm::RANLUX48:
            return (double)rl48() / (double)rl48.max();
    };
    return 0.0;
}
