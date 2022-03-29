#include "../include/TestStatistics.h"
#include "ThreadPool.h"
#include "VectorUtilities.h"
#include "ProgressBar.h"
using namespace std;

/// Thread pool used in the MultiThread runtime mode.
static thread_pool pool;

vector<double> measure_bootstrap_ks_statistic(const SyntheticPowerLawGenerator& syntheticGenerator, int replicas, RuntimeMode mode)
{
    vector<double> tsDistribution;
    tsDistribution.reserve(replicas);

    if (mode == RuntimeMode::SingleThread)
    {
        for (int i = 0; i < replicas; ++i)
        {
            progress_bar(i, 0, replicas, 1);
            tsDistribution.push_back(syntheticGenerator.MeasureKsStatisticOfReplica());
        }
        progress_bar(1.0);
    }
    else if (mode == RuntimeMode::MultiThread)
    {
        // Launch threads
        vector<future<double>> futures;
        futures.reserve(replicas);
        for (int i = 0; i < replicas; ++i)
            futures.push_back(pool.submit([&syntheticGenerator]{ return syntheticGenerator.MeasureKsStatisticOfReplica();}));

        // Get results
        for (future<double>& result : futures)
        {
            tsDistribution.push_back(result.get());
            progress_bar((double)(&result - &futures[0]), 0.0, (double) futures.size(), 1.0);
        }
        progress_bar(1.0);
    }
    return tsDistribution;
}

double calculate_gof(const DiscretePowerLawDistribution &fittedModel, const vector<int> &sampleData, int replicas,
                     SyntheticGeneratorMode syntheticGeneratorMode, RuntimeMode runtimeMode)
{
    RandomGen::Seed();

    // Error handling
    if (!fittedModel.StateIsValid())
        return 0.0;

    const double testKsValue = fittedModel.GetKSStatistic();

    // Create KS-Statistic distribution from synthetic replicas.
    SyntheticPowerLawGenerator syntheticGenerator(fittedModel, sampleData, syntheticGeneratorMode);
    vector<double> ksDistribution = measure_bootstrap_ks_statistic(syntheticGenerator, replicas, runtimeMode);

    // Measure p-value
    int syntheticLargerThanEmpirical = VectorUtilities::NumberOfGreater(ksDistribution, testKsValue);
    return (double) syntheticLargerThanEmpirical / (double) ksDistribution.size();
}