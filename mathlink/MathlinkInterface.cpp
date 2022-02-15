#include "mathlink.h"

#if defined(WINDOWS_MATHLINK)
#include <windows.h>
#endif

void test()
{

}

/****************************
*                           *
*   Mathlink entry point    *
*                           *
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