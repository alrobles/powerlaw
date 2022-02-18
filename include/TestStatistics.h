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
 * Calculates the value of the Kolmogorov-Smirnov test statistic.
 * @param empirical Reference to an empirical distribution.
 * @param model Reference to a fitted power-law model.
 * @return The numeric value of the KS statistic.
 */
double ks_statistic(const DiscreteEmpiricalDistribution& empirical, const DiscretePowerLawDistribution& model);

/**
 * Calculates the value of the Kolmogorov-Smirnov test statistic of a fitted model using the sample data.
 * @param sampleData Power-law distributed sample data.
 * @return The numeric value of the KS statistic.
 */
double calculate_ks_statistic_of_fit(const std::vector<int>& sampleData);

/**
 * Calculates the value of the Kolmogorov-Smirnov test statistic of a fitted model using an existing fitted model.
 * @param fittedModel Reference to a fitted model.
 * @param sampleData Power-law distributed sample data.
 * @return The numeric value of the KS statistic.
 */
double calculate_ks_statistic_of_fit(const DiscretePowerLawDistribution &fittedModel, const std::vector<int>& sampleData);

/**
 * Calculates the goodness of fit of the best fit model obtained from the sample data.
 * @param sampleData Power-law distributed sample data.
 * @param replicas Number of bootstrap replicas.
 * @param mode Whether run the process as a single thread or multi thread.
 * @return A p-value that represents the goodness of fit.
 */
double calculate_gof(const std::vector<int>& sampleData, int replicas = 1000, RuntimeMode mode = RuntimeMode::MultiThread);

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