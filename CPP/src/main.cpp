#include <iostream>
#include "IntegrateSimpson.h"
#include "LinSpace.h"
#include <math.h>
using namespace std;

int main(int argc, char** argv)
{
	const int nseg{128};
	vector<double> x = LinSpace((double) 0,(double) M_PI/2, nseg);

	vector<double> y{};
	for (size_t i{}; i < x.size(); i++)
	{
		y.push_back(sin(x.at(i)));
//		cout << y.at(i) << endl;
	}

	printf("The integral from 0 to pi/2 of sinx is: %.2f\n", IntegrateSimpson(y,(double) 0, (double) M_PI/2));
};
