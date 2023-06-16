#include <iostream>
#include <string>
#include <chrono>
#include "CsvParser.h"
#include "OptionParser.h"
#include "../include/TestStatistics.h"
using namespace std;

/****************************
*       Option Parser       *
****************************/

struct Arg : public option::Arg
{
    static void PrintError(const char* msg1, const option::Option& opt, const char* msg2)
    {
#if defined(_WIN32)
        fprintf(stderr, "%s", msg1);
        fwrite(opt.name, opt.namelen, 1, stderr);
        fprintf(stderr, "%s", msg2);
#endif
    }

    static option::ArgStatus Required(const option::Option& option, bool msg)
    {
        if (option.arg != 0)
            return option::ARG_OK;

        if (msg) PrintError("Option '", option, "' requires an argument.\n");
        return option::ARG_ILLEGAL;
    }
};

enum optionIndex
{
    UNKNOWN, DATA, BOOTSTRAP_REPLICAS, ALPHA_PRECISION, MODEL_TYPE, FULL_PARAMETRIC, X_PARAMETER, SINGLE_THREAD, HELP
};

const option::Descriptor usage[] =
{
        {UNKNOWN,             0, "",  "",                Arg::None,     "INSTRUCTIONS: PowerLawFitterCpp [options]\n"},
        {DATA,                0, "d", "data",            Arg::Required, "  -d <data_to_test>, \t--data=<data_to_test>  \tSample data as a list of comma-separated integers." },
        {BOOTSTRAP_REPLICAS,  0, "r", "replicas",        Arg::Required, "  -r <number_of_replicas>, \t--replicas=<number_of_replicas>  \tNumber of bootstrap replicas. Default is 2000." },
        {ALPHA_PRECISION,     0, "a", "alpha_precision", Arg::Required, "  -a <least_significant>, \t--alpha_precision=<least_significant>  \tPrecision for alpha estimation. Default is 0.01." },
        {X_PARAMETER,         0, "x", "x_parameter",     Arg::Required, "  -x <value>, \t--x_parameter=<value>  \tKnown value of the x parameter if there is any." },
        {MODEL_TYPE,          0, "m", "model_type",      Arg::Required, "  -m <type>, \t--model_type=<type>  \tType of model. Can be LeftBounded or RightBounded. Default is LeftBounded." },
        {FULL_PARAMETRIC,     0, "f", "full_parametric", Arg::None,     "  -f, \t--full_parametric  \tWhether to bootstrap using a full parametric approach. Default is semi-parametric." },
        {SINGLE_THREAD,       0, "s", "single_thread",   Arg::None,     "  -s, \t--single_thread  \tUse only one thread for the boot-strapping." },
        {HELP,                0, "",  "help",            Arg::None,     "  \t--help  \tShow instructions." },
        {0,                   0, 0,   0,                 0,             0}
};

int main(int argc, char* argv[])
{
    vector<int> data;
    int bootstrapReplicas = 2000;
    int xParameter = -1;
    double alphaPrecision = 0.01;
    RuntimeMode runtimeMode = RuntimeMode::MultiThread;
    SyntheticGeneratorMode syntheticGeneratorMode = SyntheticGeneratorMode::SemiParametric;
    DistributionType distributionType = DistributionType::LeftBounded;

    // Argument parser
    argc -= (argc > 0); argv += (argc > 0);
    const option::Stats  stats(static_cast<const option::Descriptor*>(usage), argc, argv);
    vector<option::Option> options(stats.options_max);
    vector<option::Option> buffer(stats.buffer_max);
    option::Parser parse(static_cast<const option::Descriptor*>(usage), argc, argv, &options.at(0), &buffer.at(0));

    if (parse.error())
        return 1;

    if (options.at(HELP) || argc == 0)
    {
        option::printUsage(cout, static_cast<const option::Descriptor*>(usage));
        return 0;
    }

    for (option::Option* opt = options.at(UNKNOWN); opt; opt = opt->next())
        cout << "Unknown option: " << string(opt->name, opt->namelen) << "\n";

    for (int i = 0; i < parse.optionsCount(); ++i)
    {
        const option::Option& opt = buffer.at(i);
        switch (opt.index())
        {
            case DATA:
                data = parse_csv_line<int>(opt.arg);
                break;
            case BOOTSTRAP_REPLICAS:
                bootstrapReplicas = stoi(opt.arg);
                break;
            case ALPHA_PRECISION:
                alphaPrecision = stod(opt.arg);
                break;
            case X_PARAMETER:
                xParameter = stoi(opt.arg);
                break;
            case MODEL_TYPE:
                distributionType = (string(opt.arg) == "RightBounded") ? DistributionType::RightBounded : DistributionType::LeftBounded;
                break;
            case SINGLE_THREAD:
                runtimeMode = RuntimeMode::SingleThread;
                break;
            case FULL_PARAMETRIC:
                syntheticGeneratorMode = SyntheticGeneratorMode::FullParametric;
                break;
            default:
                break;
        }
    }

    chrono::steady_clock::time_point beginTime, endTime; // Used for benchmark.
    DiscretePowerLawDistribution* model;

    if (xParameter == -1)
        model = new DiscretePowerLawDistribution(data, alphaPrecision, distributionType);
    else
        model = new DiscretePowerLawDistribution(data, xParameter, alphaPrecision, distributionType);

    cout << "Fitted model:" << endl;
    cout << "Type: " << model->GetDistributionTypeStr() << endl;
    cout << "Alpha: " << model->GetAlpha() << "±" << model->GetStandardError() << endl;

    if (model->GetDistributionType() == DistributionType::LeftBounded)
        cout << "xMin: " << model->GetXMin() << endl;
    else if (model->GetDistributionType() == DistributionType::RightBounded)
        cout << "xMax: " << model->GetXMax() << endl;

    cout << "Fit KS statistic: " << model->GetKSStatistic() << endl;
    cout << "Log-likelihood: " << model->GetLogLikelihood(data) << endl;

    beginTime = std::chrono::steady_clock::now();
    const double gof = calculate_gof(*model, data, bootstrapReplicas, syntheticGeneratorMode, runtimeMode);
    cout << "GoodnessOfFit: " << gof << endl;
    endTime = std::chrono::steady_clock::now();

    auto timePerReplica = chrono::duration_cast<chrono::microseconds>(endTime - beginTime).count();
    cout << "Benchmark: " << timePerReplica / bootstrapReplicas << " [µs] per replica" << endl;

    delete model;

    return 0;
}
