#include <vector>
#include <map>
#include "mathlink.h"
#include "../include/TestStatistics.h"
using namespace std;

double alphaPrecision = 0.01;
TestStatisticType testStatisticType = TestStatisticType::KolmogorovSmirnov;

map<string, TestStatisticType> argToType = {
        { "KolmogorovSmirnov", TestStatisticType::KolmogorovSmirnov },
        { "CramerVonMises", TestStatisticType::CramerVonMises },
        { "AndersonDarling", TestStatisticType::AndersonDarling }
};

#if defined(WINDOWS_MATHLINK)
#include <windows.h>
#endif

void set_alpha_precision(double precision)
{
    alphaPrecision = precision;
    MLPutSymbol(stdlink, "Null");
    MLEndPacket(stdlink);
}

void set_test_statistic(const char* testStatistic)
{
    string testStatisticStr = testStatistic;
    testStatisticType = argToType[testStatisticStr];
    MLPutSymbol(stdlink, "Null");
    MLEndPacket(stdlink);
}

void fit_model(int* data, long dataLength)
{
    vector<int> dataVec(data, data + dataLength);

    DiscretePowerLawDistribution model = fit_model(dataVec, alphaPrecision, testStatisticType);
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

    DiscretePowerLawDistribution model = fit_model(dataVec, xMin, alphaPrecision, testStatisticType);
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

    DiscretePowerLawDistribution model = fit_model(dataVec, alphaPrecision, testStatisticType);
    const double testStatistic = model.GetTestStatistic();
    const string testStatisticTypeName = model.GetTestStatisticTypeStr();
    const double pValue = calculate_gof(model, dataVec, replicas);

    MLPutFunction(stdlink, "Association", 3);

    MLPutFunction(stdlink, "Rule", 2);
    MLPutString(stdlink, "p-value");
    MLPutReal(stdlink, pValue);

    MLPutFunction(stdlink, "Rule", 2);
    MLPutString(stdlink, "Test-Statistic");
    MLPutReal(stdlink, testStatistic);

    MLPutFunction(stdlink, "Rule", 2);
    MLPutString(stdlink, "Test-Statistic Type");
    MLPutString(stdlink, testStatisticTypeName.c_str());

    MLEndPacket(stdlink);
}

void calculate_gof(int* data, long dataLength, int xMin, int replicas)
{
    vector<int> dataVec(data, data + dataLength);

    DiscretePowerLawDistribution model = fit_model(dataVec, xMin, alphaPrecision, testStatisticType);
    const double testStatistic = model.GetTestStatistic();
    const string testStatisticTypeName = model.GetTestStatisticTypeStr();
    const double pValue = calculate_fixed_min_gof(model, dataVec, replicas);

    MLPutFunction(stdlink, "Association", 3);

    MLPutFunction(stdlink, "Rule", 2);
    MLPutString(stdlink, "p-value");
    MLPutReal(stdlink, pValue);

    MLPutFunction(stdlink, "Rule", 2);
    MLPutString(stdlink, "Test-Statistic");
    MLPutReal(stdlink, testStatistic);

    MLPutFunction(stdlink, "Rule", 2);
    MLPutString(stdlink, "Test-Statistic Type");
    MLPutString(stdlink, testStatisticTypeName.c_str());

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