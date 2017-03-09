#ifndef scatter_plotter_h 
#define scatter_plotter_h

#include "maps.h"
#include <fstream>

using namespace std;

namespace plotter
{

class plotter
{
  public:
    plotter();
    ~plotter();

    void good_bad(vector<double>&, vector<double>&, vector<double>&, string, string, string);

};
}

#endif
