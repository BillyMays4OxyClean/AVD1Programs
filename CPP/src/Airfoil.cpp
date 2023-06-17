#include <fstream>
#include <iomanip>

#include "Airfoil.h"
#include "LinSpace.h"


Airfoil::Airfoil()
{
}

int Airfoil::ExportFile(std::string filename)
{
    ofstream dataPoints(filename);
    dataPoints << std::fixed << std::setprecision(6);

    for (size_t i{}; i < nseg; i++)
    {
        dataPoints << x.at(i) << "\t" << yu.at(i) << std::endl;
        dataPoints << x.at(i) << "\t" << yl.at(i) << std::endl;
        dataPoints << x.at(i) << "\t" << yc.at(i) << std::endl;
    }
    dataPoints.close();
    
    return 0;
}
