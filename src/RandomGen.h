#pragma once
#include <numeric>
#include <random>

class RandomGen
{
    static std::random_device rd;
    static std::mt19937 gen;

public:
    static void Seed();
    static int GetInt(int max);
    static double GetUniform01();
};