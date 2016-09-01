#include <iostream>
#include <fstream>
#include <vector>
#include "LHEF_lite.h"
#include "analyzers_2p2.h"

using namespace std;

int main()
{
    LHEF_lite reader("ttbb_h_bbbbdue.lhe");
    unsigned long neve = 0;
    //const unsigned long max_ev = -1;
    const unsigned long max_ev = 600;

    vector<pair<int, double *>> p; // particles
    vector<pair<int, int>> moth_ID;

    ofstream outfile;
    outfile.open ( "./results/MEM_output.txt" );

    ex1::analyzer2p2 a;
    while (reader.Get_event_GM(p, moth_ID) && neve < max_ev) {
        ++neve;

        a.analz(p, moth_ID, outfile);
    }

    outfile.close();

    for (auto i : p)
        delete i.second;

    return 0;
}
