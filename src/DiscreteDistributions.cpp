#include "../include/DiscreteDistributions.h"
#include "../include/TestStatistics.h"
#include "Zeta.h"
#include "VectorUtilities.h"
#include <iostream>
using namespace std;

/******************************************
*      DiscreteEmpiricalDistribution      *
******************************************/

DiscreteEmpiricalDistribution::DiscreteEmpiricalDistribution(const vector<int>& sampleData, int xMin, int xMax)
{
    // Select tail
    vector<int> sortedTailSample = sampleData;
    VectorUtilities::RemoveLower(sortedTailSample, xMin);
    VectorUtilities::RemoveGreater(sortedTailSample, xMax);
    VectorUtilities::Sort(sortedTailSample);

    // Assign and precalculate
    _xMin = xMin;
    _xMax = xMax;
    PrecalculateCDF(sortedTailSample);
}

void DiscreteEmpiricalDistribution::PrecalculateCDF(const vector<int>& sortedTailSample)
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
    _state = other._state;
    _sampleSize = other._sampleSize;
    _alphaPrecision = other._alphaPrecision;
    _testStatisticType = other._testStatisticType;
    _testStatistic = other._testStatistic;
    _distributionType = other._distributionType;
    _cdf = other._cdf;
}

DiscretePowerLawDistribution::DiscretePowerLawDistribution(const vector<int> &sampleData, int xMin, double alphaPrecision,
                                                           TestStatisticType testStatisticType)
{
    _state = InputValidator(sampleData, xMin);
    _xMin = xMin;
    _alphaPrecision = alphaPrecision;
    _testStatisticType = testStatisticType;
    _distributionType = DistributionType::LeftBounded;

    if (_state == DistributionState::Valid)
    {
        _xMax = VectorUtilities::Max(sampleData);
        _alpha = EstimateAlpha(sampleData, xMin, alphaPrecision);
        _sampleSize = VectorUtilities::NumberOfGreaterOrEqual(sampleData, xMin);

        PrecalculateCDF();
        _testStatistic = CalculateTestStatistic(sampleData, testStatisticType);
    }
}

DiscretePowerLawDistribution::DiscretePowerLawDistribution(const vector<int> &sampleData, int xMin, int xMax,
                                                           double alphaPrecision, TestStatisticType testStatisticType)
{
    _state = InputValidator(sampleData, xMin, xMax);
    _xMin = xMin;
    _xMax = xMax;
    _alphaPrecision = alphaPrecision;
    _testStatisticType = testStatisticType;
    _distributionType = DistributionType::LeftAndRightBounded;

    if (_state == DistributionState::Valid)
    {
        _alpha = EstimateAlpha(sampleData, xMin, xMax,alphaPrecision);
        _sampleSize = VectorUtilities::NumberInInterval(sampleData, xMin, xMax);

        PrecalculateCDF();
        _testStatistic = CalculateTestStatistic(sampleData, testStatisticType);
    }
}

DiscretePowerLawDistribution::DiscretePowerLawDistribution(const vector<int> &sampleData, double alphaPrecision,
                                                           TestStatisticType testStatisticType, DistributionType distributionType,
                                                           int smallestInterval)
{
    _state = InputValidator(sampleData);
    _alphaPrecision = alphaPrecision;
    _testStatisticType = testStatisticType;
    _distributionType = distributionType;

    if (_state == DistributionState::Valid)
    {
        _xMin = EstimateLowerBound(sampleData, alphaPrecision);

        if (distributionType == DistributionType::LeftBounded)
        {
            _xMax = VectorUtilities::Max(sampleData);
            _alpha = EstimateAlpha(sampleData, _xMin, alphaPrecision);
            _sampleSize = VectorUtilities::NumberOfGreaterOrEqual(sampleData, _xMin);
        }
        else if (distributionType == DistributionType::LeftAndRightBounded)
        {
            _xMax = EstimateUpperBound(sampleData, _xMin, alphaPrecision, smallestInterval);
            _alpha = EstimateAlpha(sampleData, _xMin, _xMax,alphaPrecision);
            _sampleSize = VectorUtilities::NumberInInterval(sampleData, _xMin, _xMax);
        }

        PrecalculateCDF();
        _testStatistic = CalculateTestStatistic(sampleData, testStatisticType);
    }
}

DistributionState DiscretePowerLawDistribution::InputValidator(const vector<int> &data)
{
    return !data.empty() ? DistributionState::Valid : DistributionState::NoInput;
}

DistributionState DiscretePowerLawDistribution::InputValidator(const vector<int> &data, int xMin)
{
    if (data.empty())
        return DistributionState::NoInput;

    int maxElement = VectorUtilities::Max(data);
    if (xMin >= maxElement)
        return DistributionState::InvalidInput;

    return DistributionState::Valid;
}

DistributionState DiscretePowerLawDistribution::InputValidator(const vector<int> &data, int xMin, int xMax)
{
    if (data.empty())
        return DistributionState::NoInput;

    int minElement = VectorUtilities::Min(data);
    int maxElement = VectorUtilities::Max(data);
    if (xMin >= maxElement || xMax <= minElement)
        return DistributionState::InvalidInput;

    return DistributionState::Valid;
}

void DiscretePowerLawDistribution::PrecalculateCDF()
{
    _cdf.reserve(_xMax - _xMin + 1);
    for (int x = _xMin; x <= _xMax; ++x)
    {
        double cdfValue = (_distributionType == DistributionType::LeftBounded) ?
                CalculateCDF(x, _alpha, _xMin) : CalculateCDF(x, _alpha, _xMin, _xMax);
        _cdf.push_back(cdfValue);
    }
}

double DiscretePowerLawDistribution::EstimateAlpha(const vector<int> &data, int xMin, double precision)
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

double DiscretePowerLawDistribution::EstimateAlpha(const vector<int> &data, int xMin, int xMax, double precision)
{
    const int div = static_cast<int>(1.0 / precision);
    const int lowerIntAlpha = static_cast<int>(1.50 * div);
    const int upperIntAlpha = static_cast<int>(3.51 * div);

    vector<double> logLikelihoods;
    logLikelihoods.reserve(upperIntAlpha - lowerIntAlpha);

    for (int intAlpha = lowerIntAlpha; intAlpha < upperIntAlpha; intAlpha++)
    {
        const double alpha = (double) intAlpha / div;
        logLikelihoods.push_back(CalculateLogLikelihood(data, alpha, xMin, xMax));
    }

    const int maxLikelihoodIntAlpha = VectorUtilities::IndexOfMax(logLikelihoods) + lowerIntAlpha;
    return (double) maxLikelihoodIntAlpha / div;
}

int DiscretePowerLawDistribution::EstimateLowerBound(const vector<int> &data, double precision)
{
    // Estimate xMin via finding the first local minima of KS test-statistic
    const int minElement = VectorUtilities::Min(data);
    const int maxElement = VectorUtilities::Max(data);

    double minKsStatistic = numeric_limits<double>::infinity();
    int xMinEstimator = 0;
    for (int x = minElement; x < maxElement; ++x)
    {
        const DiscretePowerLawDistribution model(data, x, precision, TestStatisticType::KolmogorovSmirnov);
        const double ksStatistic = model.GetTestStatistic();
        if (ksStatistic < minKsStatistic)
            minKsStatistic = ksStatistic;
        else
        {
            xMinEstimator = x - 1;
            break;
        }
    }

    return clamp(xMinEstimator, 1, maxElement);
}
int DiscretePowerLawDistribution::EstimateUpperBound(const vector<int> &data, int xMin, double precision, int smallestInterval)
{
    // Estimate xMin via KS minimization.
    const int minElement = xMin + smallestInterval;
    const int maxElement = VectorUtilities::Max(data);

    vector<double> ksValues;
    ksValues.reserve(maxElement - minElement);
    for (int x = minElement; x < maxElement; ++x)
    {
        const DiscretePowerLawDistribution model(data, xMin, x, precision, TestStatisticType::KolmogorovSmirnov);
        ksValues.push_back(model.GetTestStatistic());
    }

    const int xMax = VectorUtilities::IndexOfMin(ksValues) + minElement;
    return xMax;
}

double DiscretePowerLawDistribution::CalculateLogLikelihood(const vector<int> &data, double alpha, int xMin)
{
    const auto n = (double) VectorUtilities::NumberOfGreaterOrEqual(data, xMin);

    double logXSum = 0;
    for (double x: data)
        if (x >= xMin)
            logXSum += log(x);

    return - n * log(real_hurwitz_zeta(alpha, xMin)) - alpha * logXSum;
}

double DiscretePowerLawDistribution::CalculateLogLikelihood(const vector<int> &data, double alpha, int xMin, int xMax)
{
    const auto n = (double) VectorUtilities::NumberInInterval(data, xMin, xMax);

    double logXSum = 0;
    for (double x: data)
        if (x >= xMin && x <= xMax)
            logXSum += log(x);

    return - n * log(real_hurwitz_zeta(alpha, xMin) - real_hurwitz_zeta(alpha, 1 + xMax)) - alpha * logXSum;
}

double DiscretePowerLawDistribution::CalculateCDF(int x, double alpha, int xMin)
{
    if (x >= xMin)
    {
        const double numerator = real_hurwitz_zeta(alpha, x);
        const double denominator = real_hurwitz_zeta(alpha, xMin);
        return numerator / denominator;
    }
    else
        return 1.0;
}

double DiscretePowerLawDistribution::CalculateCDF(int x, double alpha, int xMin, int xMax)
{
    if (x >= xMin && x <= xMax)
    {
        const double numerator = real_hurwitz_zeta(alpha, x) - real_hurwitz_zeta(alpha, 1 + xMax);
        const double denominator = real_hurwitz_zeta(alpha, xMin) - real_hurwitz_zeta(alpha, 1 + xMax);
        return numerator / denominator;
    }
    else if (x < xMin)
        return 1.0;
    else
        return 0.0;
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

vector<int> DiscretePowerLawDistribution::GenerateRandomSequence(int n) const
{
    vector<int> randomSequence;
    randomSequence.reserve(n);

    for (int i = 0; i < n; ++i)
        randomSequence.push_back(GenerateRandomSample());

    return randomSequence;
}

int DiscretePowerLawDistribution::GenerateRandomSample() const
{
    if (_state == DistributionState::Valid)
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
    if (_state == DistributionState::Valid)
    {
        double numerator = pow(x, -_alpha);
        double denominator = real_hurwitz_zeta(_alpha, _xMin);
        return numerator / denominator;
    }
    else
        return numeric_limits<double>::quiet_NaN();
}

double DiscretePowerLawDistribution::GetCDF(int x) const
{
    if (_state == DistributionState::Valid)
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
    if (_state == DistributionState::Valid)
        return _testStatistic;
    else
        return numeric_limits<double>::infinity();
}

double DiscretePowerLawDistribution::GetAlpha() const
{
    if (_state == DistributionState::Valid)
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
    if (_state == DistributionState::Valid)
        return _xMin;
    else
        return numeric_limits<int>::quiet_NaN();
}

int DiscretePowerLawDistribution::GetXMax() const
{
    if (_state == DistributionState::Valid)
        return _xMax;
    else
        return numeric_limits<int>::quiet_NaN();
}

double DiscretePowerLawDistribution::GetStandardError() const
{
    if (_state == DistributionState::Valid)
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
    DiscreteEmpiricalDistribution empirical(data, _xMin, _xMax);

    const int xMin = _xMin;
    const int xMax = empirical.GetMaxElement();

    // Error handling
    if (xMin >= xMax || !StateIsValid())
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
    DiscreteEmpiricalDistribution empirical(data, _xMin, _xMax);

    const int xMin = _xMin;
    const int xMax = empirical.GetMaxElement();
    const int N = (int) data.size();

    // Error handling
    if (xMin >= xMax || !StateIsValid())
        return numeric_limits<double>::infinity();

    vector<double> sumElements;
    sumElements.reserve(xMax - xMin + 1);
    for (int x = xMin; x <= xMax; ++x)
        sumElements.push_back(pow(empirical.GetCDF(x) - GetCDF(x), 2.0) * GetPDF(x));

    return N * VectorUtilities::Total(sumElements);
}

double DiscretePowerLawDistribution::CalculateAndersonDarlingStatistic(const vector<int> &data) const
{
    DiscreteEmpiricalDistribution empirical(data, _xMin, _xMax);

    const int xMin = _xMin;
    const int xMax = empirical.GetMaxElement();
    const int N = (int) data.size();

    // Error handling
    if (xMin >= xMax || !StateIsValid())
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

bool DiscretePowerLawDistribution::StateIsValid() const
{
    return _state == DistributionState::Valid;
}

DistributionState DiscretePowerLawDistribution::GetState() const
{
    return _state;
}

DistributionType DiscretePowerLawDistribution::GetDistributionType() const
{
    return _distributionType;
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

vector<int> SyntheticPowerLawGenerator::SampleFromBulk(int n) const
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