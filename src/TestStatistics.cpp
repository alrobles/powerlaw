#include "../include/TestStatistics.h"
#include "ThreadPool.h"
#include "VectorUtilities.h"
using namespace std;

/// Threadpool used in the MultiThread runtime mode.
static thread_pool pool;

double measure_test_statistic_of_replica(const SyntheticPowerLawGenerator& syntheticGenerator)
{
    const vector<int> &syntheticSample = syntheticGenerator.GenerateSynthetic();
    const TestStatisticType testStatisticType = syntheticGenerator.GetModel().GetTestStatisticType();
    const DistributionType distributionType = syntheticGenerator.GetModel().GetDistributionType();
    const double alphaPrecision = syntheticGenerator.GetModel().GetAlphaPrecision();

    const DiscretePowerLawDistribution model(syntheticSample, alphaPrecision, testStatisticType, distributionType);
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

double calculate_gof(const DiscretePowerLawDistribution &fittedModel, const vector<int> &sampleData, int replicas, RuntimeMode mode)
{
    RandomGen::Seed();

    // Error handling
    if (!fittedModel.StateIsValid())
        return 0.0;

    const double testKsValue = fittedModel.GetTestStatistic();

    // Create KS-Statistic distribution from synthetic replicas.
    SyntheticPowerLawGenerator syntheticGenerator(fittedModel, sampleData);
    vector<double> ksDistribution = measure_bootstrap_test_statistic(syntheticGenerator, replicas, mode);

    // Measure p-value
    int syntheticLargerThanEmpirical = VectorUtilities::NumberOfGreater(ksDistribution, testKsValue);
    return (double) syntheticLargerThanEmpirical / (double) ksDistribution.size();
}