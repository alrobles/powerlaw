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
     * Generic discrete empirical distribution.
     * @param sampleData Any type of sample data.
     */
    explicit DiscreteEmpiricalDistribution(const std::vector<int>& sampleData);

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

/**
 * Implementation of a discrete power law distribution as described in https://arxiv.org/abs/0706.1062
 * Can be used for parameter estimation, generating a power-law distributed sample and obtaining PDF and CDF values.
 */
class DiscretePowerLawDistribution
{
private:
    double _alpha;
    double _ksStatistic;
    int _xMin, _xMax;
    int _sampleSize;
    bool _stateIsOk;
    std::vector<double> _cdf;

    [[nodiscard]] double CalculateCDF(int x) const;
    [[nodiscard]] double CalculateKSStatistic(const std::vector<int>& data) const;
    [[nodiscard]] int BinarySearch(int l, int r, double x) const;
    [[nodiscard]] double GetStandardError(int sampleSize) const;
    void PrecalculateCDF();

public:
    static double AlphaMLEEstimationApproximated(const std::vector<int>& data, int xMin);
    static double AlphaMLEEstimation(const std::vector<int>& data, int xMin);
    static double CalculateLogLikelihood(const std::vector<int>& data, double alpha, int xMin);

    /**
     * Pure parametric constructor.
     * @param alpha Known value for the alpha parameters.
     * @param xMin Known value for the xMin parameter.
     * @param xMax Maximum sample value.
     */
    DiscretePowerLawDistribution(double alpha, int xMin, int xMax);

    /**
     * Constructor for a distribution with known parameters.
     * @param sampleData Used for estimation of xMax and sampleSize. Needed to calculate StandardError.
     * @param alpha Known value for the alpha parameters.
     * @param xMin Known value for the xMin parameter.
     */
    DiscretePowerLawDistribution(const std::vector<int>& sampleData, double alpha, int xMin);

    /**
     * Constructor for a distribution with known xMin. Estimates alpha from the sample.
     * @param sampleData Sample data to estimate alpha from.
     * @param xMin Known value for the xMin parameter.
     */
    DiscretePowerLawDistribution(const std::vector<int>& sampleData, int xMin);

    /**
     * Constructor for a distribution with no known parameters. Estimates alpha and xMin from the sample data.
     * @param sampleData Data for the parameter estimation.
     */
    explicit DiscretePowerLawDistribution(const std::vector<int>& sampleData);

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

    /// Obtain the estimated standard error for alpha.
    [[nodiscard]] double GetStandardError() const;

    [[nodiscard]] double GetLogLikelihood(const std::vector<int>& data) const;

    /// Obtain the estimated xMin value.
    [[nodiscard]] int GetXMin() const;

    /// Check that there arent any error conditions.
    [[nodiscard]] int StateIsOk() const;
};

/**
 * Generator of power-law distributed synthetic replicas.
 */
class SyntheticPowerLawGenerator
{
private:
    DiscretePowerLawDistribution _powerLawDistribution;
    std::vector<int> _notInTailData;
    double _tailProbability;
    int _sampleDataSize;

    [[nodiscard]] int SampleFromNotInTail() const;
public:
    /**
     * Default constructor that takes the parameters of a fitted model.
     * @param alpha Alpha value of the fitted model.
     * @param xMin xMin value of the fitted model.
     * @param sampleData Data to extract the non-powerlaw part of the sample.
     */
    SyntheticPowerLawGenerator(double alpha, int xMin, const std::vector<int>& sampleData);

    /**
     * Pure power-law synthetic generator
     * @param alpha Alpha value of the fitted model.
     * @param xMin xMin value of the fitted model.
     * @param xMax Maximum sample value.
     * @param replicaSize Size of the synthetic replica.
     */
    SyntheticPowerLawGenerator(double alpha, int xMin, int xMax, int replicaSize);

    /// Generates a synthetic replica of the sample data.
    [[nodiscard]] std::vector<int> GenerateSynthetic() const;
};