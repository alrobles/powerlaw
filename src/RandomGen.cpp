#include "../include/RandomGen.h"
using namespace std;

random_device RandomGen::rd;
mt19937 RandomGen::gen;

void RandomGen::Seed()
{
    gen.seed(rd());
}
int RandomGen::GetInt(int max)
{
    uniform_int_distribution<> uniformIntDistribution(0, max);
    return uniformIntDistribution(gen);
}
double RandomGen::GetUniform01()
{
    uniform_real_distribution<> uniformRealDistribution(0.0, 1.0);
    return uniformRealDistribution(gen);
}
