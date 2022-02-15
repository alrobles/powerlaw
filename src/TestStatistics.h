#pragma once
#include <vector>
#include "DiscreteDistributions.h"

enum class RuntimeMode
{
    SingleThread, MultiThread
};

DiscretePowerLawDistribution fit_model(const std::vector<int>& sampleData);

double ks_statistic(const DiscreteEmpiricalDistribution& empirical, const DiscretePowerLawDistribution& model);

double calculate_ks_statistic_of_fit(const std::vector<int>& sampleData);
double calculate_ks_statistic_of_fit(const DiscretePowerLawDistribution &fittedModel, const std::vector<int>& sampleData);

double calculate_gof(const std::vector<int>& sampleData, int replicas = 1000, RuntimeMode mode = RuntimeMode::SingleThread,
                     DiscreteRandomSampleType sampleType = DiscreteRandomSampleType::Approximate);
double calculate_gof(const DiscretePowerLawDistribution& fittedModel, const std::vector<int>& sampleData,
                     int replicas = 1000, RuntimeMode mode = RuntimeMode::SingleThread,
                     DiscreteRandomSampleType sampleType = DiscreteRandomSampleType::Approximate);