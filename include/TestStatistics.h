#pragma once
#include <vector>
#include "DiscreteDistributions.h"

enum class RuntimeMode
{
    SingleThread, MultiThread
};

/// Generates a fitted model from sample data.
DiscretePowerLawDistribution fit_model(const std::vector<int>& sampleData);

/// Generates a fitted model from sample data with known xMin.
DiscretePowerLawDistribution fit_model(const std::vector<int>& sampleData, int xMin);

/**
 * Calculates the goodness of fit of a power-law model.
 * @param fittedModel Reference to the fitted power-law model.
 * @param sampleData Power-law distributed sample data.
 * @param replicas Number of bootstrap replicas.
 * @param mode Whether run the process as a single thread or multi thread.
 * @return A p-value that represents the goodness of fit.
 */
double calculate_gof(const DiscretePowerLawDistribution& fittedModel, const std::vector<int>& sampleData,
                     int replicas = 1000, RuntimeMode mode = RuntimeMode::MultiThread);

double calculate_fixed_min_gof(const DiscretePowerLawDistribution& fittedModel, const std::vector<int>& sampleData,
                               int replicas = 1000, RuntimeMode mode = RuntimeMode::MultiThread);