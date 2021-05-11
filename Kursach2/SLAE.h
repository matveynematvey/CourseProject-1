#pragma once
#include "Vector.h"
#include "Matrix.h"
#include <iomanip>

using namespace std;

class SLAE
{
public:
	int maxiter = 100;
	double eps = 1e-20;
	vector<double> xk, r, p, z, res, AU, AL;
	Matrix mat;

	SLAE()
	{
		Matrix mat = Matrix();
		xk.resize(mat.Count, 0);
		r.resize(mat.Count, 0);
		p.resize(mat.Count);
		z.resize(mat.Count);
		AU = AL = mat.globalTR;
		Borders();
		LOS();
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

	void Borders()
	{
		size_t o = 0;
		for (size_t i = 0; i <= mat.inf.Ny; i++)
			for (size_t j = 0; j <= mat.inf.Nx; j++, o++)
			{
				if (((i == 0) || (i == mat.inf.Ny)) || ((j == 0) || (j == mat.inf.Nx)))
				{
					mat.globalDI[o] = 1;
					mat.globalF[o] = u(mat.x[j], mat.y[i]);
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
		} while (sqrt(r * r) > eps && k < maxiter);
		cout << "Iterations: " << k << " " << (r * r) << " " << norm(mat.globalF - Mult(xk)) / norm(mat.globalF) << endl;
		for (int i = 0; i < mat.Count; i++)
			cout << "Count[" << i << "] " <<  setprecision(10) << xk[i] << endl;
	}
};

