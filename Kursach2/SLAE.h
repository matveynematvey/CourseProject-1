#pragma once
#include "Vector.h"
#include "Matrix.h"
#include <iomanip>

using namespace std;

class SLAE
{
public:
	int maxiter = 100;
	double eps = 1e-20, dt, dt1, dt0;
	vector<double> xk, r, p, z, res, AU, AL, q_2, q_1, mem;
	Matrix mat;
	string test1;

	SLAE() {}

	SLAE(string test)
	{
		this->test1 = test;
		mat = Matrix(test, 0);
		xk.resize(mat.Count, 0);
		q_2.resize(mat.Count, 0);
		mem.resize(mat.Count, 0);
		q_1.resize(mat.Count, 0);
		r.resize(mat.Count, 0);
		p.resize(mat.Count);
		z.resize(mat.Count);
		IterateTime();
	}

	void IterateTime()
	{
		q_2 = CalculateTrueU(0);
		q_1 = CalculateTrueU(1);
		for (size_t i = 2; i < mat.inf.Nt; i++)
		{
			mat = Matrix(test1, i);
			mem = mat.globalF;
			dt = mat.t[i] - mat.t[i - 2], dt1 = mat.t[i - 1] - mat.t[i - 2], dt0 = mat.t[i] - mat.t[i - 1];
			mat.globalDI = (((dt + dt0) / (dt * dt0)) * mat.globalMdi) + mat.globalGdi;
			mat.globalTR = (((dt + dt0) / (dt * dt0)) * mat.globalM) + mat.globalG;
			AU = AL = mat.globalTR;
			mat.globalF = mat.globalF - ((dt0 / (dt * dt1)) * MultM(q_2)) + ((dt / (dt1 * dt0)) * MultM(q_1));
			Borders(i);
			LOS();
			q_2 = q_1, q_1 = xk, mat.globalF = mem;
		}
		q_1 = CalculateTrueU(mat.inf.Nt - 1);
		//for (int i = 0; i < mat.Count; i++)
		//	cout << "Count[" << setw(log10(mat.Count)+1) << i << "] " << scientific << setw(25) <<  xk[i] << " " << setw(10) << q_1[i] << " " << setw(20) << abs(xk[i] - q_1[i]) << endl;
	}

	vector<double> Mult(const vector<double>& vec)
	{
		res.resize(mat.Count, 0);
		for (int i = 0; i < mat.Count; ++i)
		{
			int gi = mat.ig[i], gi_1 = mat.ig[i + 1];
			res[i] = mat.globalDI[i] * vec[i];
			for (int j = gi; j < gi_1; ++j)
			{
				int column = mat.jg[j];
				res[i] += AL[j] * vec[column];
				res[column] += AU[j] * vec[i];
			}
		}
		return res;
	}

	vector<double> MultM(const vector<double>& vec)
	{
		res.resize(mat.Count, 0);
		for (int i = 0; i < mat.Count; ++i)
		{
			int gi = mat.ig[i], gi_1 = mat.ig[i + 1];
			res[i] = mat.globalMdi[i] * vec[i];
			for (int j = gi; j < gi_1; ++j)
			{
				int column = mat.jg[j];
				res[i] += mat.globalM[j] * vec[column];
				res[column] += mat.globalM[j] * vec[i];
			}
		}
		return res;
	}

	void Borders(size_t time)
	{
		size_t o = 0;
		for (size_t i = 0; i <= mat.inf.Ny; i++)
			for (size_t j = 0; j <= mat.inf.Nx; j++, o++)
			{
				if (((i == 0) || (i == mat.inf.Ny)) || ((j == 0) || (j == mat.inf.Nx)))
				{
					mat.globalDI[o] = 1;
					mat.globalF[o] = mat.u(mat.x[j], mat.y[i], mat.t[time]);
					for (size_t l = mat.ig[o]; l < mat.ig[o + 1]; l++) AL[l] = 0;
					for (size_t l = o +	1; l < mat.Count; l++)
					{
						for (size_t p = mat.ig[l]; p < mat.ig[l + 1]; p++)
						{
							if (mat.jg[p] == o) AU[p] = 0;
						}
					}
				}
			}
	}


	void LOS()
	{
		xk.resize(mat.Count, 0);
		r =  Mult(xk);
		r = mat.globalF - r;
		z = r;
		p = Mult(z);
		int k = 0;
		do
		{
			double alpha = (p * r) / (p * p);
			xk = xk + (alpha * z);
			r = r - (alpha * p);
			double beta = -(p * Mult(r)) / (p * p);
			z = r + (beta * z);
			p = Mult(r) + (beta * p);
			k++;
		} while ((r * r) > eps && k < maxiter);
	}

	vector<double> CalculateTrueU(size_t index)
	{
		res.resize(mat.Count, 0);
		size_t o = 0;
		for (size_t i = 0; i <= mat.inf.Ny; i++)
			for (size_t j = 0; j <= mat.inf.Nx; j++, o++) 
				res[o] = mat.u(mat.x[j], mat.y[i], mat.t[index]);
		return res;
	}
};

