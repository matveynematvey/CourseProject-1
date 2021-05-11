#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include "Vector.h"
#include "Input.h"

using namespace std;


class Matrix
{
public:
    vector<vector<double>> D, localM, localG, constM = { {2, 1 , 1}, {1, 2, 1}, {1, 1, 2} };
    vector<int> ig, jg, index;
    vector<double> x, y, hx, hy, globalF, globalM, globalG, globalMdi, globalGdi, globalDI, globalTR;
    size_t Count, CountN;
    double detD, lambda = 1, gamma = 1;
    Input inf;

    Matrix()
    {
        Input inf = Input();
        Resize();
        Step();
        AssembleGlobalMatrices();
    }

    void Resize()
    {
        Count = ((inf.Nx + 1) * (inf.Ny + 1)); //кол-во узлов
        CountN = inf.Nx * inf.Ny * 2;          //кол-во элементов
        index.resize(3);
        x.resize(inf.Nx + 1);
        y.resize(inf.Nx + 1);
        hx.resize(3);
        hy.resize(3);
        D.resize(3, vector<double>(3));
        localG.resize(3, vector<double>(3));
        localM.resize(3, vector<double>(3));
        globalMdi.resize(Count);
        globalGdi.resize(Count);
        globalF.resize(Count);
        ig.resize(Count + 1);
        globalDI.resize(Count);
    }

    void Step()
    {
        x[0] = inf.x; y[0] = inf.y;
        for (size_t i = 1; i < inf.Nx + 1; i++)
            x[i] = x[i - 1] + inf.h;
        for (size_t i = 1; i < inf.Ny + 1; i++) 
            y[i] = y[i - 1] + inf.h;
    }

    void AssembleLocalD(int elem)
    {
        size_t floored = floor(elem / (inf.Nx * 2));
        if (elem % 2)   //нечетные элементы
        {
            index[0] = (elem + 1) / 2 + floored; index[1] = index[0] + inf.Nx; index[2] = index[1] + 1;    //номера узлов
            hx[0] = x[index[0] % 3], hx[2] = hx[0], hx[1] = hx[0] - inf.h, hy[0] = y[floored], hy[1] = hy[0] + inf.h, hy[2] = hy[1];
        }
        else            //четные эелменты
        {
            index[0] = elem / 2 + floored; index[1] = index[0] + 1; index[2] = index[0] + inf.Nx + 1;
            hy[0] = y[floored];
            hx[0] = x[index[0] % 3], hx[1] = hx[0] + inf.h, hx[2] = hx[0], hy[0] = y[floored], hy[1] = hy[0], hy[2] = hy[0] + inf.h;
        }
        detD = (hx[1] - hx[0]) * (hy[2] - hy[0]) - (hx[2] - hx[0]) * (hy[1] - hy[0]);
        D[0][0] = hx[1] * hy[2] - hx[2] * hy[1]; D[0][1] = hy[1] - hy[2]; D[0][2] = hx[2] - hx[1];
        D[1][0] = hx[2] * hy[0] - hx[0] * hy[2]; D[1][1] = hy[2] - hy[0]; D[1][2] = hx[0] - hx[2];
        D[2][0] = hx[0] * hy[1] - hx[1] * hy[0]; D[2][1] = hy[0] - hy[1]; D[2][2] = hx[1] - hx[0];
        D = (1 / detD) * D;
    }


    void AssembleLocalG()  //матрица G
    {
        for (size_t i = 0; i < 3; i++)
            for (size_t j = 0; j < 3; j++)
                localG[i][j] = (abs(detD) / 2) * lambda * (D[i][1] * D[j][1] + D[i][2] * D[j][2]);
    }

    void AssembleLocalM()
    {
        localM = (gamma * detD / 24) * constM;
    }


    void AssembleGlobalF()
    {
        for (size_t i = 0; i < 3; i++)
            globalF[index[i]] += (gamma * detD / 24) * func(hx[i],hy[i]) * (constM[i][0] + constM[i][1] + constM[i][2]);
    }


    void AssembleGlobalMatrices()
    {
        for (size_t i = 0; i < CountN; i++)
        {
            AssembleLocalD(i);
            AssembleLocalM();
            AssembleLocalG();
            AssembleGlobalF();
            FormPortrait();
        }
        globalTR.resize(ig[Count]);
        globalDI = globalMdi + globalGdi;
        globalTR = globalM + globalG;
    }

    void FormPortrait()
    {
        for (size_t i = 0; i < 3; i++)      //заполнение диагоналей
        {
            globalMdi[index[i]] += localM[i][i];
            globalGdi[index[i]] += localG[i][i];
        }
        if (!(ig[index[1] + 1] - ig[index[1]]))        //проверка на пустоту строки
        {
            PlusIG(index[1] + 1);
            jg.insert(jg.begin() + ig[index[1]], index[0]);
            globalM.insert(globalM.begin() + ig[index[1]], localM[1][0]);
            globalG.insert(globalG.begin() + ig[index[1]], localM[1][0]);
        }
        else
        {
            size_t k = ig[index[1]], i, flag = 1, eqflag = 0;
            for (i = ig[index[1]]; i < ig[index[1] + 1]; i++)
            {
                if (jg[i] == index[0])          //элемент уже есть
                {
                    globalM[i] += localM[1][0];
                    globalG[i] += localG[1][0];
                    eqflag = 1;
                    break;
                }
                if (jg[i] < index[0])
                    k = i;
                else
                    flag = 0;
            }
            if (!eqflag)
            {
                if (!flag)        //элемент нужно вставить подальше... 
                {
                    PlusIG(index[1] + 1);
                    jg.insert(jg.begin() + k + 1, index[0]);
                    globalM.insert(globalM.begin() + k + 1, localM[1][0]);
                    globalG.insert(globalG.begin() + k + 1, localG[1][0]);
                }
                else                //элемент вставить поближе....
                {
                    PlusIG(index[1] + 1);
                    jg.insert(jg.begin() + k, index[0]);
                    globalM.insert(globalM.begin() + k, localM[1][0]);
                    globalG.insert(globalG.begin() + k, localG[1][0]);
                }
            }
        }
        if (!(ig[index[2] + 1] - ig[index[2]]))
        {
            PlusIG(index[2] + 1); PlusIG(index[2] + 1);
            jg.insert(jg.begin() + ig[index[2]], index[0]);
            jg.insert(jg.begin() + ig[index[2]] + 1, index[1]);
            globalM.insert(globalM.begin() + ig[index[2]], localM[2][0]);
            globalM.insert(globalM.begin() + ig[index[2]] + 1, localM[2][1]);
            globalG.insert(globalG.begin() + ig[index[2]], localG[2][0]);
            globalG.insert(globalG.begin() + ig[index[2]] + 1, localG[2][1]);
        }
        else
        {
            for (size_t m = 0; m <= 1; m++)
            {
                size_t k = ig[index[2]], i, flag = 1, eqflag = 0;
                for (i = ig[index[2]]; i < ig[index[2] + 1]; i++)
                {
                    if (jg[i] == index[m])          //элемент уже есть
                    {
                        globalM[i] += localM[2][m];
                        globalG[i] += localG[2][m];
                        eqflag = 1;
                        break;
                    }
                    if (jg[i] < index[m])
                        k = i;
                    else
                        flag = 0;
                }
                if (!eqflag)
                {
                    if (!flag)        //элемент нужно вставить подальше... 
                    {
                        PlusIG(index[2] + 1);
                        jg.insert(jg.begin() + k + 1, index[m]);
                        globalM.insert(globalM.begin() + k + 1, localM[2][m]);
                        globalG.insert(globalG.begin() + k + 1, localG[2][m]);
                    }
                    else                //элемент вставить поближе....
                    {
                        PlusIG(index[2] + 1);
                        jg.insert(jg.begin() + k, index[m]);
                        globalM.insert(globalM.begin() + k, localM[2][m]);
                        globalG.insert(globalG.begin() + k, localG[2][m]);
                    }
                }
            }
        }

    }

    void PlusIG(int num)
    {
        for (size_t i = num; i < Count + 1; i++)
            ig[i] += 1;
    }
};

