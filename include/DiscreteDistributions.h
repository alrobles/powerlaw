#pragma once
#include <vector>
#include "RandomGen.h"

enum class DiscreteRandomSampleType
{
    Approximate, Precise
};

class DiscreteEmpiricalDistribution
{
private:
    int _xMin, _xMax;
    std::vector<double> _cdf;

    void PrecalculateTables(const std::vector<int>& sortedTailSample);
public:
    DiscreteEmpiricalDistribution(const std::vector<int>& sampleData, int xMin);
    [[nodiscard]] double GetCDF(int x) const;
    [[nodiscard]] int GetMinElement() const;
    [[nodiscard]] int GetMaxElement() const;
};

class DiscretePowerLawDistribution
{
private:
    double _alpha;
    int _xMin, _xMax;
    int _sampleSize;
    std::vector<double> _cdf;

    [[nodiscard]] int BinarySearch(int l, int r, double x) const;
    [[nodiscard]] int GetRandomNumberApproximate() const;
    [[nodiscard]] int GetRandomNumberPrecise() const;
    [[nodiscard]] double CalculateCDF(int x) const;
    static double AlphaMLEEstimation(const std::vector<int>& data, int xMin);
    void PrecalculateTables();
public:
    DiscretePowerLawDistribution(const std::vector<int>& sampleData, double alpha, int xMin);
    DiscretePowerLawDistribution(const std::vector<int>& sampleData, int xMin);
    explicit DiscretePowerLawDistribution(const std::vector<int>& sampleData);

    [[nodiscard]] std::vector<int> GenerateRandomSequence(int n, DiscreteRandomSampleType sampleType = DiscreteRandomSampleType::Approximate) const;
    [[nodiscard]] int GenerateRandomSample(DiscreteRandomSampleType sampleType = DiscreteRandomSampleType::Approximate) const;
    [[nodiscard]] double GetPDF(int x) const;
    [[nodiscard]] double GetCDF(int x) const;
    [[nodiscard]] double GetAlpha() const;
    [[nodiscard]] double GetStandardError() const;
    [[nodiscard]] double GetStandardError(int sampleSize) const;
    [[nodiscard]] int GetXMin() const;
};

class SyntheticPowerLawGenerator
{
private:
    DiscretePowerLawDistribution _powerLawDistribution;
    DiscreteRandomSampleType _sampleType;
    std::vector<int> _notInTailData;
    double _tailProbability;
    int _sampleDataSize;

    [[nodiscard]] int SampleFromNotInTail() const;
public:
    SyntheticPowerLawGenerator(double alpha, int xMin, const std::vector<int>& sampleData, DiscreteRandomSampleType sampleType);
    [[nodiscard]] std::vector<int> GenerateSynthetic() const;
};