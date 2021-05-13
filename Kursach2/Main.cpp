#include <iostream>
#include <fstream>
#include <vector>
#include "vector.h"
#include "Input.h"
#include "Matrix.h"
#include "SLAE.h"

using namespace std;

int main()
{
	SLAE test2 = SLAE("test2");
	SLAE test3 = SLAE("test3");
	SLAE test4 = SLAE("test4");
	SLAE test5 = SLAE("test5");
	cout << scientific << abs(test2.xk[test2.mat.Count / 2] - test2.q_1[test2.mat.Count / 2]) << " " << abs(test3.xk[test3.mat.Count / 2] - test3.q_1[test3.mat.Count / 2]) << " " << (abs(test2.xk[test2.mat.Count / 2] - test2.q_1[test2.mat.Count / 2]) / abs(test3.xk[test3.mat.Count / 2] - test3.q_1[test3.mat.Count / 2])) << endl;
	cout << scientific << abs(test3.xk[test3.mat.Count / 2] - test3.q_1[test3.mat.Count / 2]) << " " << abs(test4.xk[test4.mat.Count / 2] - test4.q_1[test4.mat.Count / 2]) << " " << (abs(test3.xk[test3.mat.Count / 2] - test3.q_1[test3.mat.Count / 2]) / abs(test4.xk[test4.mat.Count / 2] - test4.q_1[test4.mat.Count / 2])) << endl;
	cout << scientific << abs(test4.xk[test4.mat.Count / 2] - test4.q_1[test4.mat.Count / 2]) << " " << abs(test5.xk[test5.mat.Count / 2] - test5.q_1[test5.mat.Count / 2]) << " " << (abs(test4.xk[test4.mat.Count / 2] - test4.q_1[test4.mat.Count / 2]) / abs(test5.xk[test5.mat.Count / 2] - test5.q_1[test5.mat.Count / 2])) << endl;

}