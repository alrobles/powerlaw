#include <vector>
#include "mathlink.h"
#include "../include/TestStatistics.h"
using namespace std;

#if defined(WINDOWS_MATHLINK)
#include <windows.h>
#endif

void fit_model(int* data, long dataLength)
{
    vector<int> dataVec(data, data + dataLength);

    DiscretePowerLawDistribution model = fit_model(dataVec);
    const double alpha = model.GetAlpha();
    const double stdError = model.GetStandardError();
    const int xMin = model.GetXMin();

    MLPutFunction(stdlink, "Association", 3);

    MLPutFunction(stdlink, "Rule", 2);
    MLPutString(stdlink, "Alpha");
    MLPutReal(stdlink, alpha);

    MLPutFunction(stdlink, "Rule", 2);
    MLPutString(stdlink, "AlphaStandardError");
    MLPutReal(stdlink, stdError);

    MLPutFunction(stdlink, "Rule", 2);
    MLPutString(stdlink, "xMin");
    MLPutReal(stdlink, xMin);

    MLEndPacket(stdlink);
}

void fit_model(int* data, long dataLength, int xMin)
{
    vector<int> dataVec(data, data + dataLength);

    DiscretePowerLawDistribution model = fit_model(dataVec, xMin);
    const double alpha = model.GetAlpha();
    const double stdError = model.GetStandardError();

    MLPutFunction(stdlink, "Association", 2);

    MLPutFunction(stdlink, "Rule", 2);
    MLPutString(stdlink, "Alpha");
    MLPutReal(stdlink, alpha);

    MLPutFunction(stdlink, "Rule", 2);
    MLPutString(stdlink, "AlphaStandardError");
    MLPutReal(stdlink, stdError);

    MLEndPacket(stdlink);
}

void calculate_gof(int* data, long dataLength, int replicas)
{
    vector<int> dataVec(data, data + dataLength);

    DiscretePowerLawDistribution model = fit_model(dataVec);
    const double ksStatistic = calculate_ks_statistic_of_fit(model, dataVec);
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

void calculate_gof(int* data, long dataLength, int xMin, int replicas)
{
    vector<int> dataVec(data, data + dataLength);

    DiscretePowerLawDistribution model = fit_model(dataVec, xMin);
    const double ksStatistic = calculate_ks_statistic_of_fit(model, dataVec);
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