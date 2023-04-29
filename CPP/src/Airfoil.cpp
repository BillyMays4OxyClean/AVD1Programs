#include <fstream>
#include <iomanip>

#include "Airfoil.h"
#include "LinSpace.h"


Airfoil::Airfoil()
:
nseg(128),
c(1.0),
x(LinSpace(0.0,c,nseg)),
yu(LinSpace(0.0,0.0,nseg)),
yl(LinSpace(0.0,0.0,nseg))
{
}

int Airfoil::ExportFile(std::string filename)
{
    ofstream dataPoints(filename);
    dataPoints << std::fixed << std::setprecision(6);

    for (size_t i{}; i < nseg; i++)
    {
        dataPoints << yl.at(i) << "\t" << yu.at(i) << std::endl;
    }
    dataPoints.close();
    
    return 0;
}