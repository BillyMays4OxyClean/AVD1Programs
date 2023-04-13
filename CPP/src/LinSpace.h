#ifndef LINSPACE_H_
#define LINSPACE_H_
#include <vector>

vector<double> LinSpace(double x1, double x2, const int nseg)
{
	vector<double> double_array{};

	for (int i{}; i < nseg; i++)
	{
		double_array.push_back(i * (x2 - x1) / (nseg - 1));
	}

	return double_array;
};

#endif
