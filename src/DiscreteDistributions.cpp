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
    _ksStatistic = other._ksStatistic;
    _distributionType = other._distributionType;
    _cdf = other._cdf;
}

DiscretePowerLawDistribution::DiscretePowerLawDistribution(const vector<int> &sampleData, int xParameter, double alphaPrecision,
                                                           DistributionType distributionType)
{
    _state = InputValidator(sampleData, xParameter, distributionType);
    _alphaPrecision = alphaPrecision;
    _distributionType = distributionType;

    if (_state == DistributionState::Valid)
    {
        if (distributionType == DistributionType::LeftBounded)
        {
            _xMin = xParameter;
            _xMax = VectorUtilities::Max(sampleData);
            _alpha = EstimateAlpha(sampleData, _xMin, alphaPrecision);
            _sampleSize = VectorUtilities::NumberOfGreaterOrEqual(sampleData, _xMin);
        }
        else if (distributionType == DistributionType::RightBounded)
        {
            _xMin = 1;
            _xMax = xParameter;
            _alpha = EstimateAlpha(sampleData, _xMin, _xMax, alphaPrecision);
            _sampleSize = VectorUtilities::NumberOfLowerOrEqual(sampleData, _xMax);
        }

        PrecalculateCDF();
        _ksStatistic = CalculateKSStatistic(sampleData);
    }
}

DiscretePowerLawDistribution::DiscretePowerLawDistribution(const vector<int> &sampleData, double alphaPrecision,
                                                           DistributionType distributionType, int smallestInterval)
{
    _state = InputValidator(sampleData);
    _alphaPrecision = alphaPrecision;
    _distributionType = distributionType;

    if (_state == DistributionState::Valid)
    {
        if (distributionType == DistributionType::LeftBounded)
        {
            _xMin = EstimateLowerBound(sampleData, alphaPrecision);
            _xMax = VectorUtilities::Max(sampleData);
            _alpha = EstimateAlpha(sampleData, _xMin, alphaPrecision);
            _sampleSize = VectorUtilities::NumberOfGreaterOrEqual(sampleData, _xMin);
        }
        else if (distributionType == DistributionType::RightBounded)
        {
            _xMin = 1;
            _xMax = EstimateUpperBound(sampleData, alphaPrecision, smallestInterval);
            _alpha = EstimateAlpha(sampleData, _xMin, _xMax,alphaPrecision);
            _sampleSize = VectorUtilities::NumberOfLowerOrEqual(sampleData, _xMax);
        }

        PrecalculateCDF();
        _ksStatistic = CalculateKSStatistic(sampleData);
    }
}

DistributionState DiscretePowerLawDistribution::InputValidator(const vector<int> &data)
{
    return !data.empty() ? DistributionState::Valid : DistributionState::NoInput;
}

DistributionState DiscretePowerLawDistribution::InputValidator(const vector<int> &data, int xParameter, DistributionType distributionType)
{
    if (data.empty())
        return DistributionState::NoInput;

    if (distributionType == DistributionType::LeftBounded)
    {
        int maxElement = VectorUtilities::Max(data);
        if (xParameter >= maxElement)
            return DistributionState::InvalidInput;
    }
    else if (distributionType == DistributionType::RightBounded)
    {
        int minElement = VectorUtilities::Min(data);
        if (xParameter <= minElement)
            return DistributionState::InvalidInput;
    }

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
        logLikelihoods.push_back(CalculateLogLikelihoodLeftBounded(data, alpha, xMin));
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
        logLikelihoods.push_back(CalculateLogLikelihoodRightBounded(data, alpha, xMax));
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
        const DiscretePowerLawDistribution model(data, x, precision,
                                                 DistributionType::LeftBounded);
        const double ksStatistic = model.GetKSStatistic();
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
int DiscretePowerLawDistribution::EstimateUpperBound(const vector<int> &data, double precision, int smallestInterval)
{
    // Estimate xMin via KS minimization.
    const int minElement = 1 + smallestInterval;
    const int maxElement = VectorUtilities::Max(data);

    vector<double> ksValues;
    ksValues.reserve(maxElement - minElement);
    for (int x = minElement; x < maxElement; ++x)
    {
        const DiscretePowerLawDistribution model(data, x, precision,
                                                 DistributionType::RightBounded);
        ksValues.push_back(model.GetKSStatistic());
    }

    const int xMax = VectorUtilities::IndexOfMin(ksValues) + minElement;
    return xMax;
}

double DiscretePowerLawDistribution::CalculateLogLikelihoodLeftBounded(const vector<int> &data, double alpha, int xMin)
{
    const auto n = (double) VectorUtilities::NumberOfGreaterOrEqual(data, xMin);

    double logXSum = 0;
    for (double x: data)
        if (x >= xMin)
            logXSum += log(x);

    return - n * log(real_hurwitz_zeta(alpha, xMin)) - alpha * logXSum;
}

double DiscretePowerLawDistribution::CalculateLogLikelihoodRightBounded(const vector<int> &data, double alpha, int xMax)
{
    const auto n = (double) VectorUtilities::NumberOfLowerOrEqual(data, xMax);

    double logXSum = 0;
    for (double x: data)
        if (x >= 1 && x <= xMax)
            logXSum += log(x);

    return - n * log(real_hurwitz_zeta(alpha, 1) - real_hurwitz_zeta(alpha, 1 + xMax)) - alpha * logXSum;
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

double DiscretePowerLawDistribution::GetKSStatistic() const
{
    if (_state == DistributionState::Valid)
        return _ksStatistic;
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
    return DiscretePowerLawDistribution::CalculateLogLikelihoodLeftBounded(data, _alpha, _xMin);
}

double DiscretePowerLawDistribution::CalculateKSStatistic(const vector<int> &data) const
{
    DiscreteEmpiricalDistribution empirical(data, _xMin, _xMax);

    // Error handling
    if (!StateIsValid())
        return numeric_limits<double>::infinity();

    vector<double> diffs;
    diffs.reserve(_xMax - _xMin + 1);
    for (int x = _xMin; x <= _xMax; ++x)
        diffs.push_back(abs(empirical.GetCDF(x) - GetCDF(x)));

    const double maxDiff = VectorUtilities::Max(diffs);
    return maxDiff;
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

std::string DiscretePowerLawDistribution::GetDistributionTypeStr() const
{
    switch (_distributionType)
    {
        case DistributionType::LeftBounded:
            return "Left bounded";
        case DistributionType::RightBounded:
            return "Right bounded";
        default:
            return "<unknown>";
    }
}

/******************************************
*       SyntheticPowerLawGenerator        *
******************************************/

SyntheticPowerLawGenerator::SyntheticPowerLawGenerator(const DiscretePowerLawDistribution &model, const vector<int>& sampleData,
                                                       SyntheticGeneratorMode mode)
: _powerLawDistribution(model)
{
    _sampleDataSize = (int) sampleData.size();
    _mode = mode;

    if (mode == SyntheticGeneratorMode::SemiParametric)
    {
        _nonModelData = sampleData;
        if (model.GetDistributionType() == DistributionType::LeftBounded)
            VectorUtilities::RemoveGreaterOrEqual(_nonModelData, model.GetXMin());
        else
            VectorUtilities::RemoveLowerOrEqual(_nonModelData, model.GetXMax());

        _modelSampleProbability = 1.0 - (double) _nonModelData.size() / (double) _sampleDataSize;
    }
    else if (mode == SyntheticGeneratorMode::FullParametric)
        _modelSampleProbability = 1.0;
}

int SyntheticPowerLawGenerator::SampleFromData() const
{
    const int randomIndex = RandomGen::GetInt((int)_nonModelData.size() - 1);
    const int randomNumber = _nonModelData[randomIndex];
    return randomNumber;
}

vector<int> SyntheticPowerLawGenerator::SampleFromData(int n) const
{
    vector<int> dataSamples;
    dataSamples.reserve(n);
    for (int i = 0; i < n; ++i)
        dataSamples.push_back(SampleFromData());

    return dataSamples;
}

vector<int> SyntheticPowerLawGenerator::GenerateSynthetic() const
{
    vector<int> syntheticDataset;
    syntheticDataset.reserve(_sampleDataSize);

    const int modelSampleSize = floor(_modelSampleProbability * _sampleDataSize);
    VectorUtilities::Insert(syntheticDataset, _powerLawDistribution.GenerateRandomSequence(modelSampleSize));
    VectorUtilities::Insert(syntheticDataset, SampleFromData(_sampleDataSize - modelSampleSize));

    return syntheticDataset;
}

double SyntheticPowerLawGenerator::MeasureKsStatisticOfReplica() const
{
    const vector<int> &syntheticSample = GenerateSynthetic();
    const DistributionType distributionType = _powerLawDistribution.GetDistributionType();
    const double alphaPrecision = _powerLawDistribution.GetAlphaPrecision();

    if (_mode == SyntheticGeneratorMode::SemiParametric)
    {
        const DiscretePowerLawDistribution model(syntheticSample, alphaPrecision, distributionType);
        return model.GetKSStatistic();
    }
    else // _mode == SyntheticGeneratorMode::FullParametric
    {
        const int xParameter = (distributionType == DistributionType::LeftBounded) ?
                _powerLawDistribution.GetXMin() : _powerLawDistribution.GetXMax();
        const DiscretePowerLawDistribution model(syntheticSample, xParameter, alphaPrecision, distributionType);
        return model.GetKSStatistic();
    }
}