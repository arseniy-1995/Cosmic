#include "pch.h"
 
#include <TApplication.h>

//#include <TRandom.h>

#include <TF1.h>

using namespace System;

int main(int argc, char** argv)
{
    TApplication theApp("App", &argc, argv);


    TF1* f1 = new TF1("f1", "sin(x)/x", 0, 1);

    theApp.Run();

    return 0;
}
