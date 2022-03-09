#include "../include/TestStatistics.h"
#include "ThreadPool.h"
#include "VectorUtilities.h"
using namespace std;

/// Threadpool used in the MultiThread runtime mode.
static thread_pool pool;

DiscretePowerLawDistribution fit_model(const vector<int> &sampleData, double alphaPrecision, TestStatisticType testStatisticType)
{
    return DiscretePowerLawDistribution(sampleData, alphaPrecision, testStatisticType);
}

DiscretePowerLawDistribution fit_model(const vector<int> &sampleData, int xMin, double alphaPrecision, TestStatisticType testStatisticType)
{
    return DiscretePowerLawDistribution(sampleData, xMin, alphaPrecision, testStatisticType);
}

double measure_test_statistic_of_replica(const SyntheticPowerLawGenerator& syntheticGenerator)
{
    const vector<int> &syntheticSample = syntheticGenerator.GenerateSynthetic();
    const TestStatisticType testStatisticType = syntheticGenerator.GetModel().GetTestStatisticType();
    const double alphaPrecision = syntheticGenerator.GetModel().GetAlphaPrecision();
    const DiscretePowerLawDistribution model = fit_model(syntheticSample, alphaPrecision, testStatisticType);
    return model.GetTestStatistic();
}

double measure_test_statistic_of_replica_fixed_min(const SyntheticPowerLawGenerator& syntheticGenerator)
{
    const vector<int> &syntheticSample = syntheticGenerator.GenerateSynthetic();
    const TestStatisticType testStatisticType = syntheticGenerator.GetModel().GetTestStatisticType();
    const double alphaPrecision = syntheticGenerator.GetModel().GetAlphaPrecision();
    const int xMin = syntheticGenerator.GetModel().GetXMin();
    const DiscretePowerLawDistribution model = fit_model(syntheticSample, xMin, alphaPrecision, testStatisticType);
    return model.GetTestStatistic();
}

vector<double> measure_bootstrap_test_statistic(const SyntheticPowerLawGenerator& syntheticGenerator, int replicas, RuntimeMode mode)
{
    vector<double> tsDistribution;
    tsDistribution.reserve(replicas);

    if (mode == RuntimeMode::SingleThread)
    {
        for (int i = 0; i < replicas; ++i)
            tsDistribution.push_back(measure_test_statistic_of_replica(syntheticGenerator));
    }
    else if (mode == RuntimeMode::MultiThread)
    {
        // Launch threads
        vector<future<double>> futures;
        futures.reserve(replicas);
        for (int i = 0; i < replicas; ++i)
            futures.push_back(pool.submit(measure_test_statistic_of_replica, syntheticGenerator));

        // Get results
        for (future<double>& result : futures)
            tsDistribution.push_back(result.get());
    }
    return tsDistribution;
}

vector<double> measure_bootstrap_test_statistic_fixed_min(const SyntheticPowerLawGenerator& syntheticGenerator,
                                                          int replicas, RuntimeMode mode)
{
    vector<double> tsDistribution;
    tsDistribution.reserve(replicas);

    if (mode == RuntimeMode::SingleThread)
    {
        for (int i = 0; i < replicas; ++i)
            tsDistribution.push_back(
                    measure_test_statistic_of_replica_fixed_min(syntheticGenerator));
    }
    else if (mode == RuntimeMode::MultiThread)
    {
        // Launch threads
        vector<future<double>> futures;
        futures.reserve(replicas);
        for (int i = 0; i < replicas; ++i)
            futures.push_back(
                    pool.submit(measure_test_statistic_of_replica_fixed_min, syntheticGenerator));

        // Get results
        for (future<double>& result : futures)
            tsDistribution.push_back(result.get());
    }
    return tsDistribution;
}

double calculate_gof(const DiscretePowerLawDistribution &fittedModel, const vector<int> &sampleData, int replicas, RuntimeMode mode)
{
    RandomGen::Seed();

    // Error handling
    if (!fittedModel.StateIsOk())
        return 0.0;

    const double testKsValue = fittedModel.GetTestStatistic();

    // Create KS-Statistic distribution from synthetic replicas.
    SyntheticPowerLawGenerator syntheticGenerator(fittedModel, sampleData);
    vector<double> ksDistribution = measure_bootstrap_test_statistic(syntheticGenerator, replicas, mode);

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

    const double testKsValue = fittedModel.GetTestStatistic();

    // Create KS-Statistic distribution from synthetic replicas.
    const int n = VectorUtilities::NumberOfGreaterOrEqual(sampleData, fittedModel.GetXMin());
    SyntheticPowerLawGenerator syntheticGenerator(fittedModel, n);
    vector<double> ksDistribution = measure_bootstrap_test_statistic_fixed_min(syntheticGenerator, replicas, mode);

    // Measure p-value
    int syntheticLargerThanEmpirical = VectorUtilities::NumberOfGreater(ksDistribution, testKsValue);
    return (double) syntheticLargerThanEmpirical / (double) ksDistribution.size();
}