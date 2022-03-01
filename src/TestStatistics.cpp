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

double measure_ks_of_replica(const SyntheticPowerLawGenerator& syntheticGenerator, double alphaPrecision)
{
    const vector<int> &syntheticSample = syntheticGenerator.GenerateSynthetic();
    const DiscretePowerLawDistribution model = fit_model(syntheticSample, alphaPrecision);
    return model.GetKSStatistic();
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
            ksDistribution.push_back(measure_ks_of_replica(syntheticGenerator, fittedModel.GetAlphaPrecision()));
    }
    else if (mode == RuntimeMode::MultiThread)
    {
        // Launch threads
        vector<future<double>> futures;
        futures.reserve(replicas);
        for (int i = 0; i < replicas; ++i)
            futures.push_back(pool.submit(measure_ks_of_replica, syntheticGenerator, fittedModel.GetAlphaPrecision()));

        // Get results
        for (future<double>& result : futures)
            ksDistribution.push_back(result.get());
    }

    int syntheticLargerThanEmpirical = VectorUtilities::NumberOfGreater(ksDistribution, testKsValue);
    return (double) syntheticLargerThanEmpirical / (double) ksDistribution.size();
}