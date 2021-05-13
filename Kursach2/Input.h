#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include "Vector.h"


using namespace std;

class Input
{
public:
    double h, x, y, t0, ht;
    size_t Nx, Ny, Nt;
    string path = "test2";
    ifstream f;

    Input() {}

    Input(string test)
    {
        //cout << "Enter the folder: "; cin >> path;
        f.open(test + "/info.txt");
        f >> path >> h >> path >> Nx >> Ny >> path >> x >> y >> path >> t0 >> ht >> Nt;
        f.close();

    }
};
