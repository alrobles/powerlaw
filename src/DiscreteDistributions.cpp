#include "../include/DiscreteDistributions.h"
#include "../include/TestStatistics.h"
#include "Zeta.h"
#include "VectorUtilities.h"
#include <iostream>
using namespace std;

/******************************************
*      DiscreteEmpiricalDistribution      *
******************************************/

DiscreteEmpiricalDistribution::DiscreteEmpiricalDistribution(const vector<int>& sampleData, int xMin)
{
    // Select tail
    vector<int> sortedTailSample = sampleData;
    VectorUtilities::RemoveLower(sortedTailSample, xMin);
    VectorUtilities::Sort(sortedTailSample);

    // Assign and precalculate
    _xMin = xMin;
    _xMax = sortedTailSample.back();
    PrecalculateTables(sortedTailSample);
}

void DiscreteEmpiricalDistribution::PrecalculateTables(const std::vector<int>& sortedTailSample)
{
    const auto sortedTailSampleSize = (double) sortedTailSample.size();
    _cdf.reserve(_xMax - _xMin + 1);

    for (int x = _xMin; x <= _xMax; ++x)
    {
        const double foundIndex = VectorUtilities::IndexOf(sortedTailSample, x - 1);
        const double cdfVal = 1.0 - (foundIndex / sortedTailSampleSize);
        _cdf.push_back(cdfVal);
    }
}
double DiscreteEmpiricalDistribution::GetCDF(int x) const
{
    if (x >= _xMin && x <= _xMax)
        return _cdf[x - _xMin];
    else if (x < _xMin)
        return 1.0;
    else
        return 0.0;
}

int DiscreteEmpiricalDistribution::GetMinElement() const
{
    return _xMin;
}
int DiscreteEmpiricalDistribution::GetMaxElement() const
{
    return _xMax;
}

/******************************************
*       DiscretePowerLawDistribution      *
******************************************/

DiscretePowerLawDistribution::DiscretePowerLawDistribution(const std::vector<int>& sampleData, double alpha, int xMin)
{
    _alpha = alpha;
    _xMin = xMin;
    _xMax = *max_element(sampleData.begin(), sampleData.end());
    _sampleSize = (int) sampleData.size();
    PrecalculateTables();
}
DiscretePowerLawDistribution::DiscretePowerLawDistribution(const vector<int> &sampleData, int xMin)
{
    _xMin = xMin;
    _xMax = *max_element(sampleData.begin(), sampleData.end());
    _alpha = AlphaMLEEstimation(sampleData, xMin);
    _sampleSize = (int) sampleData.size();
    PrecalculateTables();
}
DiscretePowerLawDistribution::DiscretePowerLawDistribution(const vector<int> &sampleData)
{
    // Estimate xMin via KS minimization.
    const int minElement = VectorUtilities::Min(sampleData);
    const int maxElement = VectorUtilities::Max(sampleData);

    double minKs = std::numeric_limits<double>::infinity();
    int xMinEstimator = 0;
    for (int x = minElement; x <= maxElement; ++x)
    {
        DiscreteEmpiricalDistribution empiricalDistribution(sampleData, x);
        DiscretePowerLawDistribution model(sampleData, x);
        const double ks = ks_statistic(empiricalDistribution, model);
        if (ks < minKs)
            minKs = ks;
        else
        {
            xMinEstimator = x - 1;
            break;
        }
    }

    _xMin = xMinEstimator;
    _xMax = maxElement;
    _alpha = AlphaMLEEstimation(sampleData, _xMin);
    _sampleSize = VectorUtilities::NumberOfGreaterOrEqual(sampleData, _xMin);
    PrecalculateTables();
}

void DiscretePowerLawDistribution::PrecalculateTables()
{
    _cdf.reserve(_xMax - _xMin + 1);
    for (int x = _xMin; x <= _xMax; ++x)
        _cdf.push_back(CalculateCDF(x));
}

double DiscretePowerLawDistribution::AlphaMLEEstimation(const vector<int> &data, const int xMin)
{
    const auto n = (double) VectorUtilities::NumberOfGreaterOrEqual(data, xMin);
    double sum = 0.0;
    for (double x: data)
    {
        if (x >= xMin)
            sum += log((double) x / (xMin - 0.5));
    }

    return 1.0 + (n / sum);
}

int DiscretePowerLawDistribution::BinarySearch(int l, int r, double x) const
{
    while (l <= r)
    {
        const int mid = l + (r - l) / 2;
        const double cdf = GetCDF(mid);
        const double rCdf = GetCDF(mid + 1);
        const double lCdf = GetCDF(mid - 1);

        // If the value is between mid - 1 and mid + 1, select the one whose cdf is closer to x.
        if (x < lCdf && x > rCdf)
        {
            const double diff = abs(cdf - x);
            const double lDiff = abs(lCdf - x);
            const double rDiff = abs(rCdf - x);

            if (lDiff < diff)
                return mid - 1;
            else if (rDiff < diff)
                return mid + 1;
            else
                return mid;
        }
        // If not, proceed with the binary search.
        else if (cdf < x)
            r = mid - 1;
        else
            l = mid + 1;
    }

    return -1;
}

std::vector<int> DiscretePowerLawDistribution::GenerateRandomSequence(int n) const
{
    vector<int> randomSequence;
    randomSequence.reserve(n);

    for (int i = 0; i < n; ++i)
        randomSequence.push_back(GenerateRandomSample());

    return randomSequence;
}

int DiscretePowerLawDistribution::GenerateRandomSample() const
{
    const double r = RandomGen::GetUniform01();
    const double diff = 1 - r;

    // Find the search interval.
    int x1, x2;
    double cdf;
    x2 = _xMin;
    do
    {
        x1 = x2;
        x2 = 2 * x1;
        cdf = GetCDF(x2);
    }
    while (cdf >= diff);

    // Find exact solution in the interval by binary search
    int randomNumber = BinarySearch(x1, x2, diff);
    return randomNumber % _xMax;
}

double DiscretePowerLawDistribution::GetPDF(int x) const
{
    double numerator = pow(x, -_alpha);
    double denominator = hurwitz_zeta(_alpha, _xMin).real();
    return numerator / denominator;
}

double DiscretePowerLawDistribution::GetCDF(int x) const
{
    if (x >= _xMin && x <= _xMax)
        return _cdf[x - _xMin];
    else if (x < _xMin)
        return 1.0;
    else
        return 0.0;
}

double DiscretePowerLawDistribution::CalculateCDF(int x) const
{
    if (x >= _xMin)
    {
        const double numerator = hurwitz_zeta(_alpha, x).real();
        const double denominator = hurwitz_zeta(_alpha, _xMin).real();
        return numerator / denominator;
    }
    else
        return 1.0;
}

double DiscretePowerLawDistribution::GetAlpha() const
{
    return _alpha;
}

int DiscretePowerLawDistribution::GetXMin() const
{
    return _xMin;
}

double DiscretePowerLawDistribution::GetStandardError() const
{
    return GetStandardError(_sampleSize);
}

double DiscretePowerLawDistribution::GetStandardError(int sampleSize) const
{
    return (_alpha - 1.0) / (double) sampleSize;
}

/******************************************
*       SyntheticPowerLawGenerator        *
******************************************/

SyntheticPowerLawGenerator::SyntheticPowerLawGenerator(double alpha, int xMin, const vector<int>& sampleData)
: _powerLawDistribution(sampleData, alpha, xMin)
{
    _notInTailData = sampleData;
    _sampleDataSize = (int) sampleData.size();

    VectorUtilities::RemoveGreaterOrEqual(_notInTailData, xMin);
    _tailProbability = 1.0 - (double) _notInTailData.size() / (double) _sampleDataSize;
}

int SyntheticPowerLawGenerator::SampleFromNotInTail() const
{
    const int randomIndex = RandomGen::GetInt((int)_notInTailData.size() - 1);
    return _notInTailData[randomIndex];
}

vector<int> SyntheticPowerLawGenerator::GenerateSynthetic() const
{
    vector<int> syntheticDataset;
    syntheticDataset.reserve(_sampleDataSize);
    for (int i = 0; i < _sampleDataSize; ++i)
    {
        if (RandomGen::GetUniform01() < _tailProbability)
            syntheticDataset.push_back(_powerLawDistribution.GenerateRandomSample());
        else
            syntheticDataset.push_back(SampleFromNotInTail());
    }

    return syntheticDataset;
}