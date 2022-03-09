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
    DiscreteEmpiricalDistribution(const std::vector<int>& sampleData, int xMin);

    /// Obtain the cumulative density function at value x.
    [[nodiscard]] double GetCDF(int x) const;

    /// Returns the minimum element, either specified as xMin on the constructor or the minimum element of the sample data.
    [[nodiscard]] int GetMinElement() const;

    /// Returns the maximum element from the sample data.
    [[nodiscard]] int GetMaxElement() const;
};

enum class TestStatisticType
{
    None, KolmogorovSmirnov, CramerVonMises, AndersonDarling
};

/**
 * Implementation of a discrete power law distribution as described in https://arxiv.org/abs/0706.1062
 * Can be used for parameter estimation, generating a power-law distributed sample and obtaining PDF and CDF values.
 */
class DiscretePowerLawDistribution
{
private:
    TestStatisticType _testStatisticType;
    double _alpha;
    double _testStatistic;
    double _alphaPrecision;
    int _xMin, _xMax;
    int _sampleSize;
    bool _stateIsOk;
    std::vector<double> _cdf;

    [[nodiscard]] double CalculateCDF(int x) const;
    [[nodiscard]] double CalculateKolmogorovSmirnovStatistic(const std::vector<int>& data) const;
    [[nodiscard]] double CalculateCramerVonMisesStatistic(const std::vector<int>& data) const;
    [[nodiscard]] double CalculateAndersonDarlingStatistic(const std::vector<int>& data) const;
    [[nodiscard]] double CalculateTestStatistic(const std::vector<int>& data, TestStatisticType type) const;
    [[nodiscard]] int BinarySearch(int l, int r, double x) const;
    [[nodiscard]] double GetStandardError(int sampleSize) const;
    void PrecalculateCDF();

public:
    static double AlphaMLEEstimation(const std::vector<int>& data, int xMin, double precision = 0.01);
    static double CalculateLogLikelihood(const std::vector<int>& data, double alpha, int xMin);

    /**
     * Copy constructor
     */
    DiscretePowerLawDistribution(const DiscretePowerLawDistribution& other);

    /**
     * Constructor for a distribution with known xMin. Estimates alpha from the sample.
     * @param sampleData Sample data to estimate alpha from.
     * @param xMin Known value for the xMin parameter.
     */
    DiscretePowerLawDistribution(const std::vector<int>& sampleData, int xMin, double alphaPrecision = 0.01,
                                 TestStatisticType testStatisticType = TestStatisticType::KolmogorovSmirnov);

    /**
     * Constructor for a distribution with no known parameters. Estimates alpha and xMin from the sample data.
     * @param sampleData Data for the parameter estimation.
     */
    explicit DiscretePowerLawDistribution(const std::vector<int>& sampleData, double alphaPrecision = 0.01,
                                          TestStatisticType testStatisticType = TestStatisticType::KolmogorovSmirnov);

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

    [[nodiscard]] double GetTestStatistic() const;

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

    /// Check that there arent any error conditions.
    [[nodiscard]] int StateIsOk() const;

    [[nodiscard]] TestStatisticType GetTestStatisticType() const;

    [[nodiscard]] std::string GetTestStatisticTypeStr() const;
};

/**
 * Generator of power-law distributed synthetic replicas.
 */
class SyntheticPowerLawGenerator
{
private:
    DiscretePowerLawDistribution _powerLawDistribution;
    std::vector<int> _bulkData;
    double _tailProbability;
    int _sampleDataSize;

    [[nodiscard]] int SampleFromBulk() const;
    [[nodiscard]] std::vector<int> SampleFromBulk(int n) const;
public:
    /**
     * Default constructor that takes the parameters of a fitted model.
     * @param model Fitted model.
     * @param sampleData Data to extract the non-powerlaw part of the sample.
     */
    SyntheticPowerLawGenerator(const DiscretePowerLawDistribution& model, const std::vector<int>& sampleData);

    /**
     * Pure power-law synthetic generator
     * @param model Fitted model.
     * @param replicaSize Size of the synthetic replica.
     */
    SyntheticPowerLawGenerator(const DiscretePowerLawDistribution& model, int replicaSize);

    /// Generates a synthetic replica of the sample data.
    [[nodiscard]] std::vector<int> GenerateSynthetic() const;

    [[nodiscard]] const DiscretePowerLawDistribution& GetModel() const;
};