#include "../include/TestStatistics.h"
#include "ThreadPool.h"
#include "VectorUtilities.h"
using namespace std;

/// Threadpool used in the MultiThread runtime mode.
static thread_pool pool;

double ks_statistic(const DiscreteEmpiricalDistribution& empirical, const DiscretePowerLawDistribution& model)
{
    const int xMin = model.GetXMin();
    const int xMax = empirical.GetMaxElement();

    // Error handling
    if (xMin >= xMax || !model.StateIsOk())
        return numeric_limits<double>::infinity();

    vector<double> diffs;
    diffs.reserve(xMax - xMin + 1);
    for (int x = xMin; x <= xMax; ++x)
    {
        const double empiricalCDF = empirical.GetCDF(x);
        const double modelCDF = model.GetCDF(x);
        diffs.push_back(abs(empiricalCDF - modelCDF));
    }

    const double maxDiff = VectorUtilities::Max(diffs);
    return maxDiff;
}

double measure_ks_of_replica(const SyntheticPowerLawGenerator& syntheticGenerator)
{
    const vector<int> &syntheticSample = syntheticGenerator.GenerateSynthetic();
    return fit_model(syntheticSample).GetKSStatistic();
}

double measure_ks_of_replica_fixed_min(const SyntheticPowerLawGenerator& syntheticGenerator, int xMin)
{
    const vector<int> &syntheticSample = syntheticGenerator.GenerateSynthetic();
    return fit_model(syntheticSample, xMin).GetKSStatistic();
}

double calculate_gof(const DiscretePowerLawDistribution &fittedModel, const vector<int> &sampleData, int replicas, RuntimeMode mode)
{
    RandomGen::Seed();

    // Error handling
    if (!fittedModel.StateIsOk())
        return 0.0;

    const double testKsValue = fittedModel.GetKSStatistic();

    // Create KS-Statistic distribution from synthetic replicas.
    SyntheticPowerLawGenerator syntheticGenerator(fittedModel.GetAlpha(), fittedModel.GetXMin(), sampleData);
    vector<double> ksDistribution;
    ksDistribution.reserve(replicas);

    if (mode == RuntimeMode::SingleThread)
    {
        for (int i = 0; i < replicas; ++i)
            ksDistribution.push_back(measure_ks_of_replica(syntheticGenerator));
    }
    else if (mode == RuntimeMode::MultiThread)
    {
        // Launch threads
        vector<future<double>> futures;
        futures.reserve(replicas);
        for (int i = 0; i < replicas; ++i)
            futures.push_back(pool.submit(measure_ks_of_replica, syntheticGenerator));

        // Get results
        for (future<double>& result : futures)
            ksDistribution.push_back(result.get());
    }

    int syntheticLargerThanEmpirical = VectorUtilities::NumberOfGreater(ksDistribution, testKsValue);
    return (double) syntheticLargerThanEmpirical / (double) ksDistribution.size();
}

double calculate_fixed_min_gof(const DiscretePowerLawDistribution &fittedModel, const vector<int> &sampleData,
                               int replicas, RuntimeMode mode)
{
    RandomGen::Seed();

    // Error handling
    if (!fittedModel.StateIsOk())
        return 0.0;

    const double testKsValue = fittedModel.GetKSStatistic();

    // Create KS-Statistic distribution from synthetic replicas.
    const int n = VectorUtilities::NumberOfGreaterOrEqual(sampleData, fittedModel.GetXMin());
    const int xMax = VectorUtilities::Max(sampleData);
    SyntheticPowerLawGenerator syntheticGenerator(fittedModel.GetAlpha(), fittedModel.GetXMin(), xMax, n);
    vector<double> ksDistribution;
    ksDistribution.reserve(replicas);

    if (mode == RuntimeMode::SingleThread)
    {
        for (int i = 0; i < replicas; ++i)
            ksDistribution.push_back(measure_ks_of_replica_fixed_min(syntheticGenerator, fittedModel.GetXMin()));
    }
    else if (mode == RuntimeMode::MultiThread)
    {
        // Launch threads
        vector<future<double>> futures;
        futures.reserve(replicas);
        for (int i = 0; i < replicas; ++i)
            futures.push_back(pool.submit(measure_ks_of_replica_fixed_min, syntheticGenerator, fittedModel.GetXMin()));

        // Get results
        for (future<double>& result : futures)
            ksDistribution.push_back(result.get());
    }

    int syntheticLargerThanEmpirical = VectorUtilities::NumberOfGreater(ksDistribution, testKsValue);
    return (double) syntheticLargerThanEmpirical / (double) ksDistribution.size();
}


DiscretePowerLawDistribution fit_model(const vector<int> &sampleData)
{
    return DiscretePowerLawDistribution(sampleData);
}

DiscretePowerLawDistribution fit_model(const vector<int> &sampleData, int xMin)
{
    return DiscretePowerLawDistribution(sampleData, xMin);
}