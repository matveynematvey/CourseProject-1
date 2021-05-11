#pragma once
#include <iomanip>
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

class Func
{
public:
    Func() {}

    //double lambda(double u)
    //{
    //    return u + 1;
    //}

    double u(double x, double y, double t)
    {
        return x + y;
    }

    double func(double x, double y)
    {
        return 0;
    }
};
