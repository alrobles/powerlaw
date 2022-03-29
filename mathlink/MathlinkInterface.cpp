#include <vector>
#include <string>
#include <map>
#include "mathlink.h"
#include "../include/TestStatistics.h"
using namespace std;

double alphaPrecision = 0.01;

#if defined(WINDOWS_MATHLINK)
#include <windows.h>
#endif

void set_alpha_precision(double precision)
{
    alphaPrecision = precision;
    MLPutSymbol(stdlink, "Null");
    MLEndPacket(stdlink);
}

void fit_model(int* data, long dataLength, const char* distributionType)
{
    vector<int> dataVec(data, data + dataLength);
    DistributionType distTypeEnum = (string(distributionType) == "RightBounded") ?
                                    DistributionType::RightBounded : DistributionType::LeftBounded;

    DiscretePowerLawDistribution model(dataVec, alphaPrecision, distTypeEnum);
    const double alpha = model.GetAlpha();
    const double stdError = model.GetStandardError();
    const int xMin = model.GetXMin();
    const int xMax = model.GetXMax();

    MLPutFunction(stdlink, "Association", 4);

    MLPutFunction(stdlink, "Rule", 2);
    MLPutString(stdlink, "Alpha");
    MLPutReal(stdlink, alpha);

    MLPutFunction(stdlink, "Rule", 2);
    MLPutString(stdlink, "AlphaStandardError");
    MLPutReal(stdlink, stdError);

    MLPutFunction(stdlink, "Rule", 2);
    MLPutString(stdlink, "xMin");
    MLPutInteger(stdlink, xMin);

    MLPutFunction(stdlink, "Rule", 2);
    MLPutString(stdlink, "xMax");
    MLPutInteger(stdlink, xMax);

    MLEndPacket(stdlink);
}

void fit_model(int* data, long dataLength, int xParameter, const char* distributionType)
{
    vector<int> dataVec(data, data + dataLength);
    DistributionType distTypeEnum = (string(distributionType) == "RightBounded") ?
            DistributionType::RightBounded : DistributionType::LeftBounded;

    DiscretePowerLawDistribution model(dataVec, xParameter, alphaPrecision, distTypeEnum);
    const double alpha = model.GetAlpha();
    const double stdError = model.GetStandardError();
    const int xMin = model.GetXMin();
    const int xMax = model.GetXMax();

    MLPutFunction(stdlink, "Association", 4);

    MLPutFunction(stdlink, "Rule", 2);
    MLPutString(stdlink, "Alpha");
    MLPutReal(stdlink, alpha);

    MLPutFunction(stdlink, "Rule", 2);
    MLPutString(stdlink, "AlphaStandardError");
    MLPutReal(stdlink, stdError);

    MLPutFunction(stdlink, "Rule", 2);
    MLPutString(stdlink, "xMin");
    MLPutInteger(stdlink, xMin);

    MLPutFunction(stdlink, "Rule", 2);
    MLPutString(stdlink, "xMax");
    MLPutInteger(stdlink, xMax);

    MLEndPacket(stdlink);
}

void calculate_gof(int* data, long dataLength, int replicas)
{
    vector<int> dataVec(data, data + dataLength);

    DiscretePowerLawDistribution model(dataVec, alphaPrecision);
    const double ksStatistic = model.GetKSStatistic();
    const double pValue = calculate_gof(model, dataVec, replicas);

    MLPutFunction(stdlink, "Association", 2);

    MLPutFunction(stdlink, "Rule", 2);
    MLPutString(stdlink, "p-value");
    MLPutReal(stdlink, pValue);

    MLPutFunction(stdlink, "Rule", 2);
    MLPutString(stdlink, "KS-Statistic");
    MLPutReal(stdlink, ksStatistic);

    MLEndPacket(stdlink);
}

void calculate_gof(int* data, long dataLength, int replicas, int xParameter, const char* distributionType,
                   const char* syntheticGeneratorMode)
{
    vector<int> dataVec(data, data + dataLength);
    DistributionType distTypeEnum = (string(distributionType) == "RightBounded") ?
                                    DistributionType::RightBounded : DistributionType::LeftBounded;
    SyntheticGeneratorMode syntheticGeneratorModeEnum = (string(syntheticGeneratorMode) == "FullParametric") ?
                                                        SyntheticGeneratorMode::FullParametric :
                                                        SyntheticGeneratorMode::SemiParametric;

    DiscretePowerLawDistribution model(dataVec, xParameter, alphaPrecision, distTypeEnum);
    const double ksStatistic = model.GetKSStatistic();
    const double pValue = calculate_gof(model, dataVec, replicas, syntheticGeneratorModeEnum);

    MLPutFunction(stdlink, "Association", 2);

    MLPutFunction(stdlink, "Rule", 2);
    MLPutString(stdlink, "p-value");
    MLPutReal(stdlink, pValue);

    MLPutFunction(stdlink, "Rule", 2);
    MLPutString(stdlink, "KS-Statistic");
    MLPutReal(stdlink, ksStatistic);

    MLEndPacket(stdlink);
}

/****************************
*   Mathlink entry point    *
****************************/

#if defined(WINDOWS_MATHLINK)

extern HWND MLInitializeIcon(HINSTANCE hinstCurrent, int nCmdShow);

int WINAPI WinMain(HINSTANCE hinstCurrent, HINSTANCE hinstPrevious, LPSTR lpszCmdLine, int nCmdShow)
{
    char buff[512];
    char* buff_start = buff;
    char* argv[32];
    char** argv_end = argv + 32;

    hinstPrevious = hinstPrevious; /* suppress warning */

    if (!MLInitializeIcon(hinstCurrent, nCmdShow)) return 1;
    MLScanString(argv, &argv_end, &lpszCmdLine, &buff_start);
    return MLMain((int)(argv_end - argv), argv);
}

#else

int main(int argc, char* argv[])
{
    return MLMain(argc, argv);
}

#endif