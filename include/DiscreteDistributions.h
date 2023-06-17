#pragma once
#include <vector>
#include "RandomGen.h"

/**
 * Storage for a discrete empirical distribution with truncated xMin.
 * CDF values are calculated on initialization for fast runtime access.
 */
class DiscreteEmpiricalDistribution
{
private:
    int _xMin, _xMax;
    std::vector<double> _cdf;

    void PrecalculateCDF(const std::vector<int>& sortedTailSample);
public:
    /**
     * Power-law discrete empirical distribution with known xMin.
     * @param sampleData Power-law distributed sample data.
     * @param xMin Known cut-off value of xMin.
     */
    DiscreteEmpiricalDistribution(const std::vector<int>& sampleData, int xMin, int xMax);

    /// Obtain the cumulative density function at value x.
    [[nodiscard]] double GetCDF(int x) const;
};

enum class DistributionType
{
    LeftBounded, // Type I
    RightBounded // Type II
};

enum class DistributionState
{
    Valid, NoInput, InvalidInput
};

/**
 * Implementation of a discrete power law distribution as described in https://arxiv.org/abs/0706.1062
 * Can be used for parameter estimation, generating a power-law distributed sample and obtaining PDF and CDF values.
 */
class DiscretePowerLawDistribution
{
private:
    DistributionType _distributionType;
    DistributionState _state;
    double _alpha;
    double _ksStatistic;
    double _alphaPrecision;
    int _xMin, _xMax;
    int _sampleSize;
    std::vector<double> _cdf;

    static DistributionState InputValidator(const std::vector<int>& data);
    static DistributionState InputValidator(const std::vector<int>& data, int xParameter, DistributionType distributionType);

    /**
     * Estimate Alpha for model type I
     * @param data Sample data
     * @param xMin Known xMin
     * @param precision Multiple of the desired alpha precision.
     * @return The estimated value for alpha
     */
    static double EstimateAlpha(const std::vector<int>& data, int xMin, double precision = 0.01);

    /**
     * Estimate Alpha for model type II
     * @param data Sample data
     * @param xMin Known xMin
     * @param xMax Known xMax
     * @param precision Multiple of the desired alpha precision
     * @return The estimated value for alpha
     */
    static double EstimateAlpha(const std::vector<int>& data, int xMin, int xMax, double precision = 0.01);

    /**
     * Calculate the estimated value for xMin
     * @param data Sample data.
     * @param precision Multiple of the desired alpha precision
     * @return xMin value
     */
    static int EstimateLowerBound(const std::vector<int>& data, double precision = 0.01);

    /**
     * Calculate the estimated value for xMax
     * @param data Sample data
     * @param precision Multiple of the desired alpha precision
     * @param smallestInterval Minimum xMax-xMin interval
     * @return xMax value
     */
    static int EstimateUpperBound(const std::vector<int>& data, double precision = 0.01, int smallestInterval = 20);

    /// Log-likelihood for model type I
    static double CalculateLogLikelihoodLeftBounded(const std::vector<int>& data, double alpha, int xMin);

    /// Log-likelihood for model type II
    static double CalculateLogLikelihoodRightBounded(const std::vector<int>& data, double alpha, int xMax);

    /// Calculates the CDF for the model type I
    static double CalculateCDF(int x, double alpha, int xMin);

    /// Calculates the CDF for the model type II
    static double CalculateCDF(int x, double alpha, int xMin, int xMax);

    [[nodiscard]] double CalculateKSStatistic(const std::vector<int>& data) const;
    [[nodiscard]] int BinarySearch(int l, int r, double x) const;
    [[nodiscard]] double GetStandardError(int sampleSize) const;

    /// Precomputes the cumulative distribution function for fast access
    void PrecalculateCDF();

public:
    /**
     * Copy constructor
     */
    DiscretePowerLawDistribution(const DiscretePowerLawDistribution& other);

    /**
     * Constructor for a type I distribution with known xParameter. Estimates alpha from the sample.
     * @param sampleData Sample data to estimate alpha from
     * @param xParameter Known value for the xParameter parameter
     * @param alphaPrecision Multiple of the desired alpha precision
     * @param testStatisticType The type of test statistic that will be used for gof estimation
     */
    DiscretePowerLawDistribution(const std::vector<int>& sampleData, int xParameter, double alphaPrecision = 0.01,
                                 DistributionType distributionType = DistributionType::LeftBounded);

    /**
     * Constructor for a distribution with no known parameters. Estimates alpha and xMin from the sample data.
     * @param sampleData Data for the parameter estimation.
     */
    explicit DiscretePowerLawDistribution(const std::vector<int>& sampleData, double alphaPrecision = 0.01,
                                          DistributionType distributionType = DistributionType::LeftBounded,
                                          int smallestInterval = 20);

    /**
     * Generates a sequence of n power-law distributed random numbers.
     * @param n The length of the sequence to generate.
     * @return A vector sequence of power-law distributed integers.
     */
    [[nodiscard]] std::vector<int> GenerateRandomSequence(int n) const;

    /**
     * Generates one power-law distributed sample.
     * @return An integer drawn from a power-law distributed random variable.
     */
    [[nodiscard]] int GenerateRandomSample() const;

    /// Calculate the probability density function at value x.
    [[nodiscard]] double GetPDF(int x) const;

    /// Obtain the cumulative density function at value x.
    [[nodiscard]] double GetCDF(int x) const;

    [[nodiscard]] double GetKSStatistic() const;

    /// Obtain the estimated alpha value.
    [[nodiscard]] double GetAlpha() const;

    /// Obtain the precision in which estimate alpha values.
    [[nodiscard]] double GetAlphaPrecision() const;

    /// Obtain the estimated standard error for alpha.
    [[nodiscard]] double GetStandardError() const;

    // Obtain the log-likelihood that the sample was drawn from the model.
    [[nodiscard]] double GetLogLikelihood(const std::vector<int>& data) const;

    /// Obtain the estimated xMin value.
    [[nodiscard]] int GetXMin() const;

    /// Obtain the estimated xMax value.
    [[nodiscard]] int GetXMax() const;

    /// Check that there arent any error conditions.
    [[nodiscard]] bool StateIsValid() const;

    [[nodiscard]] DistributionState GetState() const;

    [[nodiscard]] DistributionType GetDistributionType() const;

    [[nodiscard]] std::string GetDistributionTypeStr() const;
};

enum class SyntheticGeneratorMode
{
    SemiParametric, FullParametric
};

/**
 * Generator of power-law distributed synthetic replicas.
 */
class SyntheticPowerLawGenerator
{
private:
    DiscretePowerLawDistribution _powerLawDistribution;
    SyntheticGeneratorMode _mode;
    std::vector<int> _nonModelData;
    double _modelSampleProbability;
    int _sampleDataSize;

    [[nodiscard]] int SampleFromData() const;
    [[nodiscard]] std::vector<int> SampleFromData(int n) const;
public:
    /**
     * Default constructor that takes the parameters of a fitted model.
     * @param model Fitted model.
     * @param sampleData Data to extract the non-powerlaw part of the sample.
     */
    SyntheticPowerLawGenerator(const DiscretePowerLawDistribution& model, const std::vector<int>& sampleData,
                               SyntheticGeneratorMode mode = SyntheticGeneratorMode::SemiParametric);

    /// Generates a synthetic replica of the sample data.
    [[nodiscard]] std::vector<int> GenerateSynthetic() const;
    [[nodiscard]] double MeasureKsStatisticOfReplica() const;
};