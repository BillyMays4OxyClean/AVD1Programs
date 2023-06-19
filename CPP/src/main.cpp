#include <iostream>
#include <math.h>
#include <string>
#include <vector>


#include "IntegrateSimpson.h"
#include "NACA4.h"

std::vector<double> LinSpace(double x1, double x2, int nseg);

using namespace std;

int main(int argc, char** argv)
{
    const int nseg{128};
    vector<double> x = LinSpace((double) 0,(double) M_PI/2, nseg);

    vector<double> y{};
    for (size_t i{}; i < x.size(); i++)
    {
        y.push_back(sin(x.at(i)));
    }

    printf("The integral from 0 to pi/2 of sinx is: %.2f\n", IntegrateSimpson(y,(double) 0, (double) M_PI/2));

    NACA4 foil(string("2312"), nseg); // NACA 2312

    foil.ExportFile("Upper.dat",foil.x,foil.yu);
    foil.ExportFile("Lower.dat",foil.x,foil.yl);
    foil.ExportFile("MeanChord.dat",foil.x,foil.yc);

};
