#include "../include/DiscreteDistributions.h"
#include "../include/TestStatistics.h"
#include "Zeta.h"
#include "VectorUtilities.h"
#include <iostream>
using namespace std;

/******************************************
*      DiscreteEmpiricalDistribution      *
******************************************/

DiscreteEmpiricalDistribution::DiscreteEmpiricalDistribution(const vector<int>& sampleData)
{
    // Select tail
    vector<int> sortedTailSample = sampleData;
    VectorUtilities::Sort(sortedTailSample);

    // Assign and precalculate
    _xMin = sortedTailSample.front();
    _xMax = sortedTailSample.back();
    PrecalculateCDF(sortedTailSample);
}

DiscreteEmpiricalDistribution::DiscreteEmpiricalDistribution(const vector<int>& sampleData, int xMin)
{
    // Select tail
    vector<int> sortedTailSample = sampleData;
    VectorUtilities::RemoveLower(sortedTailSample, xMin);
    VectorUtilities::Sort(sortedTailSample);

    // Assign and precalculate
    _xMin = xMin;
    _xMax = sortedTailSample.back();
    PrecalculateCDF(sortedTailSample);
}

void DiscreteEmpiricalDistribution::PrecalculateCDF(const std::vector<int>& sortedTailSample)
{
    if (_xMax > _xMin)
    {
        const auto sortedTailSampleSize = (double) sortedTailSample.size();
        _cdf.reserve(_xMax - _xMin + 1);

        _cdf.push_back(1.0);
        for (int x = _xMin; x < _xMax; ++x)
        {
            const double foundIndex = VectorUtilities::IndexOf(sortedTailSample, x);
            const double cdfVal = 1.0 - (foundIndex / sortedTailSampleSize);
            _cdf.push_back(cdfVal);
        }
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

DiscretePowerLawDistribution::DiscretePowerLawDistribution(double alpha, int xMin, int xMax)
{
    _alpha = alpha;
    _xMin = xMin;
    _xMax = xMax;
    _stateIsOk = (_xMin < _xMax);
    _sampleSize = -1;

    if (_stateIsOk)
        PrecalculateCDF();
}

DiscretePowerLawDistribution::DiscretePowerLawDistribution(const std::vector<int>& sampleData, double alpha, int xMin)
{
    _alpha = alpha;
    _xMin = xMin;
    _xMax = VectorUtilities::Max(sampleData);
    _stateIsOk = (_xMin < _xMax) && !sampleData.empty();
    _sampleSize = (int) sampleData.size();

    if (_stateIsOk)
        PrecalculateCDF();

    // Calculate KS
    DiscreteEmpiricalDistribution empiricalDistribution(sampleData, xMin);
    _ksStatistic = ks_statistic(empiricalDistribution, *this);
}
DiscretePowerLawDistribution::DiscretePowerLawDistribution(const vector<int> &sampleData, int xMin)
{
    _xMin = xMin;
    if (!sampleData.empty())
    {
        _xMax = VectorUtilities::Max(sampleData);
        _stateIsOk = (_xMin < _xMax);
        _alpha = AlphaMLEEstimation(sampleData, xMin);
        _sampleSize = (int) sampleData.size();
    }
    else
        _stateIsOk = false;

    if (_stateIsOk)
        PrecalculateCDF();

    // Calculate KS
    DiscreteEmpiricalDistribution empiricalDistribution(sampleData, xMin);
    _ksStatistic = ks_statistic(empiricalDistribution, *this);
}
DiscretePowerLawDistribution::DiscretePowerLawDistribution(const vector<int> &sampleData)
{
    _stateIsOk = !sampleData.empty();

    if (_stateIsOk)
    {
        // Estimate xMin via KS minimization.
        const int minElement = VectorUtilities::Min(sampleData);
        const int maxElement = VectorUtilities::Max(sampleData);

        double minKs = std::numeric_limits<double>::infinity();
        int xMinEstimator = 0;
        for (int x = minElement; x < maxElement; ++x)
        {
            DiscretePowerLawDistribution model(sampleData, x);
            const double ks = model.GetKSStatistic();
            if (ks < minKs)
                minKs = ks;
            else
            {
                _ksStatistic = minKs;
                xMinEstimator = x - 1;
                break;
            }
        }

        _xMin = clamp(xMinEstimator, 1, maxElement);
        _xMax = maxElement;
        _alpha = AlphaMLEEstimation(sampleData, _xMin);
        _sampleSize = VectorUtilities::NumberOfGreaterOrEqual(sampleData, _xMin);
        PrecalculateCDF();
    }
}

void DiscretePowerLawDistribution::PrecalculateCDF()
{
    _cdf.reserve(_xMax - _xMin + 1);
    for (int x = _xMin; x <= _xMax; ++x)
        _cdf.push_back(CalculateCDF(x));
}

double DiscretePowerLawDistribution::AlphaMLEEstimationApproximated(const vector<int> &data, int xMin)
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

double DiscretePowerLawDistribution::AlphaMLEEstimation(const vector<int> &data, int xMin)
{
    const int lowerIntAlpha = 150;
    const int upperIntAlpha = 351;
    vector<double> logLikelihoods;
    logLikelihoods.reserve(upperIntAlpha - lowerIntAlpha);

    for (int intAlpha = lowerIntAlpha; intAlpha < upperIntAlpha; intAlpha++)
    {
        const double alpha = (double) intAlpha / 100.0;
        logLikelihoods.push_back(CalculateLogLikelihood(data, alpha, xMin));
    }

    const int maxLikelihoodIntAlpha = VectorUtilities::IndexOfMax(logLikelihoods) + lowerIntAlpha;
    return (double) maxLikelihoodIntAlpha / 100.0;
}

double DiscretePowerLawDistribution::CalculateLogLikelihood(const vector<int> &data, double alpha, int xMin)
{
    const auto n = (double) VectorUtilities::NumberOfGreaterOrEqual(data, xMin);

    double logXSum = 0;
    for (double x: data)
        if (x >= xMin)
            logXSum += log(x);

    return - n * log(hurwitz_zeta(alpha, xMin).real()) - alpha * logXSum;
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

            if (lDiff < diff)
                return mid - 1;
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
    if (_stateIsOk)
    {
        const double r = RandomGen::GetUniform01();

        // Find the search interval.
        int x1, x2;
        double cdf;
        x2 = _xMin;
        do
        {
            x1 = x2;
            x2 = 2 * x1;
            cdf = GetCDF(x2);
        } while (cdf >= r);

        // Find exact solution in the interval by binary search
        const int randomNumber = BinarySearch(x1, x2, r);
        return randomNumber;
    }
    else
        return numeric_limits<int>::quiet_NaN();
}

double DiscretePowerLawDistribution::GetPDF(int x) const
{
    if (_stateIsOk)
    {
        double numerator = pow(x, -_alpha);
        double denominator = hurwitz_zeta(_alpha, _xMin).real();
        return numerator / denominator;
    }
    else
        return numeric_limits<double>::quiet_NaN();
}

double DiscretePowerLawDistribution::GetCDF(int x) const
{
    if (_stateIsOk)
    {
        if (x >= _xMin && x <= _xMax)
            return _cdf[x - _xMin];
        else if (x < _xMin)
            return 1.0;
        else
            return 0.0;
    }
    else
        return numeric_limits<double>::quiet_NaN();
}

double DiscretePowerLawDistribution::GetKSStatistic() const
{
    return _ksStatistic;
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
    if (_stateIsOk)
        return _alpha;
    else
        return numeric_limits<double>::quiet_NaN();
}

int DiscretePowerLawDistribution::GetXMin() const
{
    if (_stateIsOk)
        return _xMin;
    else
        return numeric_limits<int>::quiet_NaN();
}

double DiscretePowerLawDistribution::GetStandardError() const
{
    if (_stateIsOk)
        return GetStandardError(_sampleSize);
    else
        return numeric_limits<double>::quiet_NaN();
}

double DiscretePowerLawDistribution::GetStandardError(int sampleSize) const
{
    return (_alpha - 1.0) / (double) sampleSize;
}

double DiscretePowerLawDistribution::GetLogLikelihood(const vector<int> &data) const
{
    return DiscretePowerLawDistribution::CalculateLogLikelihood(data, _alpha, _xMin);
}

int DiscretePowerLawDistribution::StateIsOk() const
{
    return _stateIsOk;
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

SyntheticPowerLawGenerator::SyntheticPowerLawGenerator(double alpha, int xMin, int xMax, int replicaSize)
: _powerLawDistribution(alpha, xMin, xMax)
{
    _sampleDataSize = replicaSize;
    _tailProbability = 1.0;
}


int SyntheticPowerLawGenerator::SampleFromNotInTail() const
{
    const int randomIndex = RandomGen::GetInt((int)_notInTailData.size() - 1);
    const int randomNumber = _notInTailData[randomIndex];
    return randomNumber;
}

vector<int> SyntheticPowerLawGenerator::GenerateSynthetic() const
{
    vector<int> syntheticDataset;
    syntheticDataset.reserve(_sampleDataSize);
    for (int i = 0; i < _sampleDataSize; ++i)
    {
        const int randomNumber = (RandomGen::GetUniform01() < _tailProbability) ?
                                 _powerLawDistribution.GenerateRandomSample() : SampleFromNotInTail();
        syntheticDataset.push_back(randomNumber);
    }

    return syntheticDataset;
}