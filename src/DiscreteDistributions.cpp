#include "DiscreteDistributions.h"
#include "TestStatistics.h"
#include "Zeta.h"
#include "RandomGen.h"
#include <iostream>
using namespace std;

/******************************************
*      DiscreteEmpiricalDistribution      *
******************************************/

DiscreteEmpiricalDistribution::DiscreteEmpiricalDistribution(const vector<int>& sampleData)
{
    PrecalculateTables(sampleData);
}
void DiscreteEmpiricalDistribution::PrecalculateTables(vector<int> sampleData)
{
    sort(sampleData.begin(), sampleData.end());
    _minElement = sampleData.front();
    _maxElement = sampleData.back();
    _cdf.reserve(_maxElement - _minElement);

    for (int i = _minElement; i < _maxElement; i++)
    {
        double foundIndex = (double)(upper_bound(sampleData.begin(), sampleData.end(), i - 1) - sampleData.begin());
        _cdf.push_back(1.0 - (foundIndex / ((double) sampleData.size())));
    }
}
double DiscreteEmpiricalDistribution::GetCDF(int x) const
{
    if (x > 0 && x < _cdf.size())
        return _cdf[x - 1];
    else if (x <= 0)
        return 1.0;
    else
        return 0.0;
}

int DiscreteEmpiricalDistribution::GetMaxElement() const
{
    return _maxElement;
}

/******************************************
*       DiscretePowerLawDistribution      *
******************************************/

DiscretePowerLawDistribution::DiscretePowerLawDistribution(double alpha, int xMin, int sampleSize)
{
    _alpha = alpha;
    _xMin = xMin;
    _sampleSize = sampleSize;
}
DiscretePowerLawDistribution::DiscretePowerLawDistribution(const vector<int> &sampleData, int xMin)
{
    _xMin = xMin;
    _alpha = AlphaMLEEstimation(sampleData, xMin);
    _sampleSize = (int) sampleData.size();
}
DiscretePowerLawDistribution::DiscretePowerLawDistribution(const vector<int> &sampleData)
{
    // Estimate xMin via KS minimization.
    DiscreteEmpiricalDistribution empiricalDistribution(sampleData);
    const int maxElement = empiricalDistribution.GetMaxElement();
    double minKs = std::numeric_limits<double>::infinity();
    int xMinEstimator = 0;
    for (int x = 1; x < maxElement; x++)
    {
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
    _alpha = AlphaMLEEstimation(sampleData, _xMin);
    _sampleSize = (int) sampleData.size();
}

double DiscretePowerLawDistribution::AlphaMLEEstimation(const vector<int> &data, const int xMin)
{
    const auto n = (double) data.size();
    double sum = 0.0;
    for (double x: data)
    {
        if (x >= xMin)
            sum += log((double) x / (xMin - 0.5));
    }

    return 1.0 + (n / sum);
}

int DiscretePowerLawDistribution::GetRandomNumberApproximate() const
{
    double r = RandomGen::GetUniform01();
    r = clamp(r, 0.0, 0.9999); // Hack to avoid excessively high values.
    const double realRand = (_xMin - 0.5) * pow((1.0 - r), -1.0 / (_alpha - 1.0)) + 0.5;
    return (int) round(realRand);
}

int DiscretePowerLawDistribution::BinarySearch(int l, int r, double x) const
{
    if (l <= r)
    {
        int mid = l + (r - l) / 2;
        const double cdf = GetCDF(mid);
        const double rCdf = GetCDF(mid + 1);
        const double lCdf = GetCDF(mid - 1);
        const double diff = abs(cdf - x);
        const double lDiff = abs(lCdf - x);
        const double rDiff = abs(rCdf - x);

        // If the value is between mid - 1 and mid + 1, select the one whose cdf is closer to x.
        if (x < lCdf && x > rCdf)
        {
            if (lDiff < diff)
                return mid - 1;
            else if (rDiff < diff)
                return mid + 1;
            else
                return mid;
        }
        // If not, proceed with the binary search.
        else if (cdf < x)
            return BinarySearch(l, mid - 1, x);
        else
            return BinarySearch(mid + 1, r, x);
    }
    else
        return -1;
}

int DiscretePowerLawDistribution::GetRandomNumberPrecise() const
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
    return BinarySearch(x1, x2, diff);
}

std::vector<int> DiscretePowerLawDistribution::GenerateRandomSequence(int n, DiscreteRandomSampleType sampleType) const
{
    vector<int> randomSequence;
    randomSequence.reserve(n);

    if (sampleType == DiscreteRandomSampleType::Approximate)
    {
        for (int i = 0; i < n; i++)
            randomSequence.push_back(GetRandomNumberApproximate());
    }
    else if (sampleType == DiscreteRandomSampleType::Precise)
    {
        for (int i = 0; i < n; i++)
            randomSequence.push_back(GetRandomNumberPrecise());
    }
    return randomSequence;
}

int DiscretePowerLawDistribution::GenerateRandomSample(DiscreteRandomSampleType sampleType) const
{
    return (sampleType == DiscreteRandomSampleType::Approximate) ? GetRandomNumberApproximate() : GetRandomNumberPrecise();
}

double DiscretePowerLawDistribution::GetPDF(int x) const
{
    double numerator = pow(x, -_alpha);
    double denominator = hurwitz_zeta(_alpha, _xMin).real();
    return numerator / denominator;
}

double DiscretePowerLawDistribution::GetCDF(int x) const
{
    if (x > 0)
    {
        const double numerator = hurwitz_zeta(_alpha, x).real();
        const double denominator = hurwitz_zeta(_alpha, _xMin).real();
        const double cdf = numerator / denominator;
        return cdf;
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

SyntheticPowerLawGenerator::SyntheticPowerLawGenerator(double alpha, int xMin, const vector<int> &sampleData)
: _powerLawDistribution(alpha, xMin, (int)sampleData.size())
{
    _notInTailData = sampleData;
    _sampleDataSize = (int) sampleData.size();

    // Remove elements of the power-law tail.
    auto it = remove_if(_notInTailData.begin(), _notInTailData.end(),
                        [&](auto const& val){ return val >= xMin; });
    _notInTailData.erase(it, _notInTailData.end());

    _tailProbability = 1.0 - (double) _notInTailData.size() / (double) sampleData.size();
}

int SyntheticPowerLawGenerator::SampleFromNotInTail() const
{
    const int randomIndex = RandomGen::GetInt((int)_notInTailData.size());
    return _notInTailData[randomIndex];
}

vector<int> SyntheticPowerLawGenerator::GenerateSynthetic(DiscreteRandomSampleType sampleType) const
{
    vector<int> syntheticDataset;
    syntheticDataset.reserve(_sampleDataSize);
    for (int i = 0; i < _sampleDataSize; i++)
    {
        if (RandomGen::GetUniform01() < _tailProbability)
            syntheticDataset.push_back(_powerLawDistribution.GenerateRandomSample(sampleType));
        else
            syntheticDataset.push_back(SampleFromNotInTail());
    }

    return syntheticDataset;
}