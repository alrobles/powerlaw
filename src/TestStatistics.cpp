#include "../include/TestStatistics.h"
#include "ThreadPool.h"
#include "VectorUtilities.h"
using namespace std;

/// Threadpool used in the MultiThread runtime mode.
static thread_pool pool;

DiscretePowerLawDistribution fit_model(const vector<int> &sampleData, double alphaPrecision)
{
    return DiscretePowerLawDistribution(sampleData, alphaPrecision);
}

DiscretePowerLawDistribution fit_model(const vector<int> &sampleData, int xMin, double alphaPrecision)
{
    return DiscretePowerLawDistribution(sampleData, xMin, alphaPrecision);
}

double measure_ks_of_replica(const SyntheticPowerLawGenerator& syntheticGenerator, double alphaPrecision)
{
    const vector<int> &syntheticSample = syntheticGenerator.GenerateSynthetic();
    const DiscretePowerLawDistribution model = fit_model(syntheticSample, alphaPrecision);
    return model.GetKSStatistic();
}

double measure_ks_of_replica_fixed_min(const SyntheticPowerLawGenerator& syntheticGenerator, int xMin, double alphaPrecision)
{
    const vector<int> &syntheticSample = syntheticGenerator.GenerateSynthetic();
    const DiscretePowerLawDistribution model = fit_model(syntheticSample, xMin, alphaPrecision);
    return model.GetKSStatistic();
}

vector<double> measure_bootstrap_ks(const SyntheticPowerLawGenerator& syntheticGenerator, int replicas, RuntimeMode mode,
                                    double alphaPrecision)
{
    vector<double> ksDistribution;
    ksDistribution.reserve(replicas);

    if (mode == RuntimeMode::SingleThread)
    {
        for (int i = 0; i < replicas; ++i)
            ksDistribution.push_back(measure_ks_of_replica(syntheticGenerator, alphaPrecision));
    }
    else if (mode == RuntimeMode::MultiThread)
    {
        // Launch threads
        vector<future<double>> futures;
        futures.reserve(replicas);
        for (int i = 0; i < replicas; ++i)
            futures.push_back(pool.submit(measure_ks_of_replica, syntheticGenerator, alphaPrecision));

        // Get results
        for (future<double>& result : futures)
            ksDistribution.push_back(result.get());
    }
    return ksDistribution;
}

vector<double> measure_bootstrap_ks_fixed_min(const SyntheticPowerLawGenerator& syntheticGenerator, int xMin,
                                              int replicas, RuntimeMode mode, double alphaPrecision)
{
    vector<double> ksDistribution;
    ksDistribution.reserve(replicas);

    if (mode == RuntimeMode::SingleThread)
    {
        for (int i = 0; i < replicas; ++i)
            ksDistribution.push_back(measure_ks_of_replica_fixed_min(syntheticGenerator, xMin, alphaPrecision));
    }
    else if (mode == RuntimeMode::MultiThread)
    {
        // Launch threads
        vector<future<double>> futures;
        futures.reserve(replicas);
        for (int i = 0; i < replicas; ++i)
            futures.push_back(pool.submit(measure_ks_of_replica_fixed_min, syntheticGenerator, xMin, alphaPrecision));

        // Get results
        for (future<double>& result : futures)
            ksDistribution.push_back(result.get());
    }
    return ksDistribution;
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
    vector<double> ksDistribution = measure_bootstrap_ks(syntheticGenerator, replicas, mode, fittedModel.GetAlphaPrecision());

    // Measure p-value
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
    const int xMin = fittedModel.GetXMin();
    const int xMax = VectorUtilities::Max(sampleData);
    SyntheticPowerLawGenerator syntheticGenerator(fittedModel.GetAlpha(), xMin, xMax, n);
    vector<double> ksDistribution = measure_bootstrap_ks_fixed_min(syntheticGenerator, xMin, replicas, mode,
                                                                   fittedModel.GetAlphaPrecision());

    // Measure p-value
    int syntheticLargerThanEmpirical = VectorUtilities::NumberOfGreater(ksDistribution, testKsValue);
    return (double) syntheticLargerThanEmpirical / (double) ksDistribution.size();
}