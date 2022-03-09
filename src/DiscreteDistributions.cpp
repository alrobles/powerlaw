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

DiscretePowerLawDistribution::DiscretePowerLawDistribution(const DiscretePowerLawDistribution &other)
{
    _alpha = other._alpha;
    _xMin = other._xMin;
    _xMax = other._xMax;
    _stateIsOk = other._stateIsOk;
    _sampleSize = other._sampleSize;
    _alphaPrecision = other._alphaPrecision;
    _testStatisticType = other._testStatisticType;
    _testStatistic = other._testStatistic;
    _cdf = other._cdf;
}

DiscretePowerLawDistribution::DiscretePowerLawDistribution(const vector<int> &sampleData, int xMin, double alphaPrecision,
                                                           TestStatisticType testStatisticType)
{
    _xMin = xMin;
    _alphaPrecision = alphaPrecision;
    _testStatisticType = testStatisticType;

    if (!sampleData.empty())
    {
        _xMax = VectorUtilities::Max(sampleData);
        _stateIsOk = (_xMin < _xMax);
        _alpha = AlphaMLEEstimation(sampleData, xMin, alphaPrecision);
        _sampleSize = (int) sampleData.size();
    }
    else
        _stateIsOk = false;

    if (_stateIsOk)
    {
        PrecalculateCDF();
        _testStatistic = CalculateTestStatistic(sampleData, testStatisticType);
    }
}
DiscretePowerLawDistribution::DiscretePowerLawDistribution(const vector<int> &sampleData, double alphaPrecision,
                                                           TestStatisticType testStatisticType)
{
    _stateIsOk = !sampleData.empty();
    _alphaPrecision = alphaPrecision;
    _testStatisticType = testStatisticType;

    if (_stateIsOk)
    {
        // Estimate xMin via KS minimization.
        const int minElement = VectorUtilities::Min(sampleData);
        const int maxElement = VectorUtilities::Max(sampleData);

        double minKsStatistic = std::numeric_limits<double>::infinity();
        int xMinEstimator = 0;
        for (int x = minElement; x < maxElement; ++x)
        {
            const DiscretePowerLawDistribution model(sampleData, x, alphaPrecision, TestStatisticType::KolmogorovSmirnov);
            const double ksStatistic = model.GetTestStatistic();
            if (ksStatistic < minKsStatistic)
                minKsStatistic = ksStatistic;
            else
            {
                xMinEstimator = x - 1;
                break;
            }
        }

        _xMin = clamp(xMinEstimator, 1, maxElement);
        _xMax = maxElement;
        _alpha = AlphaMLEEstimation(sampleData, _xMin, alphaPrecision);
        _sampleSize = VectorUtilities::NumberOfGreaterOrEqual(sampleData, _xMin);
        PrecalculateCDF();
        _testStatistic = CalculateTestStatistic(sampleData, testStatisticType);
    }
}

void DiscretePowerLawDistribution::PrecalculateCDF()
{
    _cdf.reserve(_xMax - _xMin + 1);
    for (int x = _xMin; x <= _xMax; ++x)
        _cdf.push_back(CalculateCDF(x));
}

double DiscretePowerLawDistribution::AlphaMLEEstimation(const vector<int> &data, int xMin, double precision)
{
    const int div = static_cast<int>(1.0 / precision);
    const int lowerIntAlpha = static_cast<int>(1.50 * div);
    const int upperIntAlpha = static_cast<int>(3.51 * div);

    vector<double> logLikelihoods;
    logLikelihoods.reserve(upperIntAlpha - lowerIntAlpha);

    for (int intAlpha = lowerIntAlpha; intAlpha < upperIntAlpha; intAlpha++)
    {
        const double alpha = (double) intAlpha / div;
        logLikelihoods.push_back(CalculateLogLikelihood(data, alpha, xMin));
    }

    const int maxLikelihoodIntAlpha = VectorUtilities::IndexOfMax(logLikelihoods) + lowerIntAlpha;
    return (double) maxLikelihoodIntAlpha / div;
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
            if (x < lCdf && x > cdf)
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
        return BinarySearch(x1, x2, r);
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

double DiscretePowerLawDistribution::GetTestStatistic() const
{
    if (_stateIsOk)
        return _testStatistic;
    else
        return numeric_limits<double>::infinity();
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

double DiscretePowerLawDistribution::GetAlphaPrecision() const
{
    return _alphaPrecision;
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

double DiscretePowerLawDistribution::CalculateTestStatistic(const vector<int> &data, TestStatisticType type) const
{
    switch (type)
    {
        case TestStatisticType::KolmogorovSmirnov:
            return CalculateKolmogorovSmirnovStatistic(data);
        case TestStatisticType::CramerVonMises:
            return CalculateCramerVonMisesStatistic(data);
        case TestStatisticType::AndersonDarling:
            return CalculateAndersonDarlingStatistic(data);
        default:
            return numeric_limits<double>::infinity();
    }
}

double DiscretePowerLawDistribution::CalculateKolmogorovSmirnovStatistic(const vector<int> &data) const
{
    DiscreteEmpiricalDistribution empirical(data, _xMin);

    const int xMin = _xMin;
    const int xMax = empirical.GetMaxElement();

    // Error handling
    if (xMin >= xMax || !StateIsOk())
        return numeric_limits<double>::infinity();

    vector<double> diffs;
    diffs.reserve(xMax - xMin + 1);
    for (int x = xMin; x <= xMax; ++x)
        diffs.push_back(abs(empirical.GetCDF(x) - GetCDF(x)));

    const double maxDiff = VectorUtilities::Max(diffs);
    return maxDiff;
}

double DiscretePowerLawDistribution::CalculateCramerVonMisesStatistic(const vector<int> &data) const
{
    DiscreteEmpiricalDistribution empirical(data, _xMin);

    const int xMin = _xMin;
    const int xMax = empirical.GetMaxElement();
    const int N = (int) data.size();

    // Error handling
    if (xMin >= xMax || !StateIsOk())
        return numeric_limits<double>::infinity();

    vector<double> sumElements;
    sumElements.reserve(xMax - xMin + 1);
    for (int x = xMin; x <= xMax; ++x)
        sumElements.push_back(pow(empirical.GetCDF(x) - GetCDF(x), 2.0) * GetPDF(x));

    return N * VectorUtilities::Total(sumElements);
}

double DiscretePowerLawDistribution::CalculateAndersonDarlingStatistic(const vector<int> &data) const
{
    DiscreteEmpiricalDistribution empirical(data, _xMin);

    const int xMin = _xMin;
    const int xMax = empirical.GetMaxElement();
    const int N = (int) data.size();

    // Error handling
    if (xMin >= xMax || !StateIsOk())
        return numeric_limits<double>::infinity();

    vector<double> sumElements;
    sumElements.reserve(xMax - xMin + 1);
    for (int x = xMin; x <= xMax; ++x)
    {
        const double numerator = pow(empirical.GetCDF(x) - GetCDF(x), 2.0) * GetPDF(x);
        const double denominator = GetCDF(x) * (1.0 - GetCDF(x));
        if (denominator != 0)
            sumElements.push_back(numerator / denominator);
        else
            sumElements.push_back(0.0);
    }
    return N * VectorUtilities::Total(sumElements);
}

int DiscretePowerLawDistribution::StateIsOk() const
{
    return _stateIsOk;
}

TestStatisticType DiscretePowerLawDistribution::GetTestStatisticType() const
{
    return _testStatisticType;
}

string DiscretePowerLawDistribution::GetTestStatisticTypeStr() const
{
    switch(_testStatisticType)
    {
        case TestStatisticType::KolmogorovSmirnov:
            return "Kolmogorov-Smirnov";
        case TestStatisticType::CramerVonMises:
            return "Cramer-VonMises";
        case TestStatisticType::AndersonDarling:
            return "Anderson-Darling";
        case TestStatisticType::None:
        default:
            return "<unknown>";
    }
}

/******************************************
*       SyntheticPowerLawGenerator        *
******************************************/

SyntheticPowerLawGenerator::SyntheticPowerLawGenerator(const DiscretePowerLawDistribution &model, const vector<int>& sampleData)
: _powerLawDistribution(model)
{
    _bulkData = sampleData;
    VectorUtilities::RemoveGreaterOrEqual(_bulkData, model.GetXMin());

    _sampleDataSize = (int) sampleData.size();
    _tailProbability = 1.0 - (double) _bulkData.size() / (double) _sampleDataSize;
}

SyntheticPowerLawGenerator::SyntheticPowerLawGenerator(const DiscretePowerLawDistribution &model, int replicaSize)
: _powerLawDistribution(model)
{
    _sampleDataSize = replicaSize;
    _tailProbability = 1.0;
}

int SyntheticPowerLawGenerator::SampleFromBulk() const
{
    const int randomIndex = RandomGen::GetInt((int)_bulkData.size() - 1);
    const int randomNumber = _bulkData[randomIndex];
    return randomNumber;
}

std::vector<int> SyntheticPowerLawGenerator::SampleFromBulk(int n) const
{
    vector<int> bulkSamples;
    bulkSamples.reserve(n);
    for (int i = 0; i < n; ++i)
        bulkSamples.push_back(SampleFromBulk());

    return bulkSamples;
}

vector<int> SyntheticPowerLawGenerator::GenerateSynthetic() const
{
    vector<int> syntheticDataset;
    syntheticDataset.reserve(_sampleDataSize);

    const int tailSampleSize = floor(_tailProbability * _sampleDataSize);
    VectorUtilities::Insert(syntheticDataset, _powerLawDistribution.GenerateRandomSequence(tailSampleSize));
    VectorUtilities::Insert(syntheticDataset, SampleFromBulk(_sampleDataSize - tailSampleSize));

    return syntheticDataset;
}

const DiscretePowerLawDistribution& SyntheticPowerLawGenerator::GetModel() const
{
    return _powerLawDistribution;
}