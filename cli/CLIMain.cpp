#include <iostream>
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
    UNKNOWN, DATA, BOOTSTRAP_REPLICAS, XMIN, ALPHA_MLE_APPROXIMATE, SINGLE_THREAD, HELP
};

const option::Descriptor usage[] =
{
        {UNKNOWN, 0, "", "",Arg::None, "INSTRUCTIONS: PowerLawFitterCpp [options]\n"},
        {DATA, 0,"d", "data", Arg::Required, "  -d <data_to_test>, \t--data=<data_to_test>  \tSample data as a list of comma-separated integers." },
        {BOOTSTRAP_REPLICAS, 0,"r", "replicas", Arg::Required, "  -r <number_of_replicas>, \t--replicas=<number_of_replicas>  \tNumber of bootstrap replicas. Default is 2000." },
        {XMIN, 0,"x", "x_min", Arg::Required, "  -x <xMin>, \t--x_min=<xMin>  \tKnown value of xMin if there is any." },
        {ALPHA_MLE_APPROXIMATE,  0, "A", "--alpha_mle_approximate", Arg::None, "  -a, \t--alpha_mle_approximate  \tEstimate alpha with the approximate formula (faster but less accurate)." },
        {SINGLE_THREAD,  0, "s", "single_thread", Arg::None, "  -s, \t--single_thread  \tUse only one thread for the boot-strapping." },
        {HELP, 0,"", "help", Arg::None,    "  \t--help  \tShow instructions." },
        {0,0,0,0,0,0}
};

int main(int argc, char* argv[])
{
    vector<int> data;
    int bootstrapReplicas = 2000;
    int xMin = -1;
    RuntimeMode runtimeMode = RuntimeMode::MultiThread;

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
            case XMIN:
                xMin = stoi(opt.arg);
                break;
            case SINGLE_THREAD:
                runtimeMode = RuntimeMode::SingleThread;
                break;
            default:
                break;
        }
    }

    if (xMin == -1)
    {
        DiscretePowerLawDistribution model = fit_model(data);

        cout << "Fitted model:" << endl;
        cout << "Alpha: " << model.GetAlpha() << "±" << model.GetStandardError() << " xMin: " << model.GetXMin() << endl;
        cout << "Fit KS statistic: " << model.GetKSStatistic() << endl;
        cout << "Log-likelihood: " << model.GetLogLikelihood(data) << endl;
        cout << "GoodnessOfFit: " << calculate_gof(model, data, bootstrapReplicas, runtimeMode) << endl;
    }
    else
    {
        DiscretePowerLawDistribution model = fit_model(data, xMin);

        cout << "Fitted model:" << endl;
        cout << "Alpha: " << model.GetAlpha() << "±" << model.GetStandardError() << endl;
        cout << "Fit KS statistic: " << model.GetKSStatistic() << endl;
        cout << "Log-likelihood: " << model.GetLogLikelihood(data) << endl;
        cout << "GoodnessOfFit: " << calculate_fixed_min_gof(model, data, bootstrapReplicas, runtimeMode) << endl;
    }

    return 0;
}
