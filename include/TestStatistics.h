#pragma once
#include <vector>
#include "DiscreteDistributions.h"

enum class RuntimeMode
{
    SingleThread, MultiThread
};

/**
 * Calculates the goodness of fit of a power-law model.
 * @param fittedModel Reference to the fitted power-law model.
 * @param sampleData Power-law distributed sample data.
 * @param replicas Number of bootstrap replicas.
 * @param runtimeMode Whether run the process as a single thread or multi thread.
 * @return A p-value that represents the goodness of fit.
 */
double calculate_gof(const DiscretePowerLawDistribution& fittedModel, const std::vector<int>& sampleData,
                     int replicas = 1000, SyntheticGeneratorMode syntheticGeneratorMode = SyntheticGeneratorMode::SemiParametric,
                     RuntimeMode runtimeMode = RuntimeMode::MultiThread);