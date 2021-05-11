#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include "Vector.h"

using namespace std;

class Input
{
public:
    double h, Nx, Ny, x, y;
    string path = "test2";
    ifstream f;

    Input()
    {
        //cout << "Enter the folder: "; cin >> path;
        f.open(path + "/info.txt");
        f >> path >> h >> path >> Nx >> Ny >> path >> x >> y;
        f.close();
    }
};
