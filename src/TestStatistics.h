#pragma once
#include <vector>
#include "DiscreteDistributions.h"

double ks_statistic(const DiscreteEmpiricalDistribution& empirical, const DiscretePowerLawDistribution& model);
double calculate_gof(const DiscretePowerLawDistribution& fittedModel, const std::vector<int> &sampleData, int replicas = 1000);
DiscretePowerLawDistribution fit_model(const std::vector<int>& sampleData);