#include "TestStatistics.h"
#include "ThreadPool.h"
using namespace std;

double ks_statistic(const DiscreteEmpiricalDistribution& empirical, const DiscretePowerLawDistribution& model)
{
    const int xMin = model.GetXMin();
    const int xMax = empirical.GetMaxElement();

    vector<double> diffs;
    diffs.reserve(xMax - xMin);
    for (int x = xMin; x < xMax; x++)
        diffs.push_back(abs(empirical.GetCDF(x) - model.GetCDF(x)));

    const double maxDiff = *max_element(diffs.begin(), diffs.end());
    return maxDiff;
}

double measure_ks_of_replica(const DiscretePowerLawDistribution& fittedModel, const SyntheticPowerLawGenerator& syntheticGenerator)
{
    const vector<int>& syntheticSample = syntheticGenerator.GenerateSynthetic(DiscreteRandomSampleType::Precise);
    DiscreteEmpiricalDistribution empiricalSynthetic(syntheticSample);
    return ks_statistic(empiricalSynthetic, fittedModel);
}

double calculate_gof(const DiscretePowerLawDistribution &fittedModel, const vector<int> &sampleData, int replicas, RuntimeMode mode)
{
    // Calculate KS-Statistic value of the fitted model.
    DiscreteEmpiricalDistribution empirical(sampleData);
    double testKsValue = ks_statistic(empirical, fittedModel);

    // Create KS-Statistic distribution from synthetic replicas.
    SyntheticPowerLawGenerator syntheticGenerator(fittedModel.GetAlpha(), fittedModel.GetXMin(), sampleData);
    vector<double> ksDistribution;
    ksDistribution.reserve(replicas);

    if (mode == RuntimeMode::SingleThread)
    {
        for (int i = 0; i < replicas; i++)
            ksDistribution.push_back(measure_ks_of_replica(fittedModel, syntheticGenerator));
    }
    else if (mode == RuntimeMode::MultiThread)
    {
        // Launch threads
        thread_pool pool;
        vector<future<double>> futures;
        futures.reserve(replicas);
        for (int i = 0; i < replicas; i++)
            futures.push_back(pool.submit(measure_ks_of_replica, fittedModel, syntheticGenerator));

        // Get results
        for (future<double>& result : futures)
            ksDistribution.push_back(result.get());
    }

    int syntheticLargerThanEmpirical = (int) count_if(ksDistribution.begin(), ksDistribution.end(),
                                                      [&](auto const& val){ return val > testKsValue; });
    return (double) syntheticLargerThanEmpirical / (double) ksDistribution.size();
}

DiscretePowerLawDistribution fit_model(const vector<int> &sampleData)
{
    return DiscretePowerLawDistribution(sampleData);
}
