#pragma once
#include <vector>

enum class DiscreteRandomSampleType
{
    Approximate, Precise
};

class DiscreteEmpiricalDistribution
{
private:
    std::vector<double> _cdf;
    int _minElement, _maxElement;

    void PrecalculateTables(std::vector<int> sampleData);
public:
    explicit DiscreteEmpiricalDistribution(const std::vector<int>& sampleData);
    [[nodiscard]] double GetCDF(int x) const;
    [[nodiscard]] int GetMaxElement() const;
};

class DiscretePowerLawDistribution
{
private:
    double _alpha;
    int _xMin;
    int _sampleSize;

    [[nodiscard]] int BinarySearch(int l, int r, double x) const;
    [[nodiscard]] int GetRandomNumberApproximate() const;
    [[nodiscard]] int GetRandomNumberPrecise() const;
    static double AlphaMLEEstimation(const std::vector<int>& data, int xMin);
public:
    DiscretePowerLawDistribution(double alpha, int xMin, int sampleSize = -1);
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
    std::vector<int> _notInTailData;
    double _tailProbability;
    int _sampleDataSize;

    [[nodiscard]] int SampleFromNotInTail() const;
public:
    SyntheticPowerLawGenerator(double alpha, int xMin, const std::vector<int>& sampleData);
    [[nodiscard]] std::vector<int> GenerateSynthetic(DiscreteRandomSampleType sampleType = DiscreteRandomSampleType::Approximate) const;
};