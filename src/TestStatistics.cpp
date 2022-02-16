#include "../include/TestStatistics.h"
#include "ThreadPool.h"
#include "VectorUtilities.h"
using namespace std;

static thread_pool pool;

double ks_statistic(const DiscreteEmpiricalDistribution& empirical, const DiscretePowerLawDistribution& model)
{
    const int xMin = model.GetXMin();
    const int xMax = empirical.GetMaxElement();

    if (xMin >= xMax)
        return numeric_limits<double>::infinity();

    vector<double> diffs;
    diffs.reserve(xMax - xMin + 1);
    for (int x = xMin; x <= xMax; ++x)
        diffs.push_back(abs(empirical.GetCDF(x) - model.GetCDF(x)));

    const double maxDiff = VectorUtilities::Max(diffs);
    return maxDiff;
}

double measure_ks_of_replica(const DiscretePowerLawDistribution& fittedModel, const SyntheticPowerLawGenerator& syntheticGenerator)
{
    const vector<int>& syntheticSample = syntheticGenerator.GenerateSynthetic();
    DiscreteEmpiricalDistribution empiricalSynthetic(syntheticSample, fittedModel.GetXMin());
    return ks_statistic(empiricalSynthetic, fittedModel);
}

double calculate_ks_statistic_of_fit(const vector<int> &sampleData)
{
    return calculate_ks_statistic_of_fit(fit_model(sampleData), sampleData);
}
double calculate_ks_statistic_of_fit(const DiscretePowerLawDistribution &fittedModel, const vector<int> &sampleData)
{
    DiscreteEmpiricalDistribution empirical(sampleData, fittedModel.GetXMin());
    return ks_statistic(empirical, fittedModel);
}

double calculate_gof(const vector<int>& sampleData, int replicas, RuntimeMode mode, DiscreteRandomSampleType sampleType)
{
    return calculate_gof(fit_model(sampleData), sampleData, replicas, mode, sampleType);
}
double calculate_gof(const DiscretePowerLawDistribution &fittedModel, const vector<int> &sampleData, int replicas,
                     RuntimeMode mode, DiscreteRandomSampleType sampleType)
{
    RandomGen::Seed();

    // Calculate KS-Statistic value of the fitted model.
    DiscreteEmpiricalDistribution empirical(sampleData, fittedModel.GetXMin());
    double testKsValue = ks_statistic(empirical, fittedModel);

    // Create KS-Statistic distribution from synthetic replicas.
    SyntheticPowerLawGenerator syntheticGenerator(fittedModel.GetAlpha(), fittedModel.GetXMin(), sampleData, sampleType);
    vector<double> ksDistribution;
    ksDistribution.reserve(replicas);

    if (mode == RuntimeMode::SingleThread)
    {
        for (int i = 0; i < replicas; ++i)
            ksDistribution.push_back(measure_ks_of_replica(fittedModel, syntheticGenerator));
    }
    else if (mode == RuntimeMode::MultiThread)
    {
        // Launch threads
        vector<future<double>> futures;
        futures.reserve(replicas);
        for (int i = 0; i < replicas; ++i)
            futures.push_back(pool.submit(measure_ks_of_replica, fittedModel, syntheticGenerator));

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
