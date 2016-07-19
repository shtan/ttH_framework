#include <iostream>
#include <vector>
#include "LHEF_lite.h"
#include "analyzers_2p2.h"

using namespace std;

int main()
{
    LHEF_lite reader("ttbb_h_bbbbdue.lhe");
    unsigned long neve = 0;
    const unsigned long max_ev = -1;

    vector<pair<int, double *>> p; // particles
    vector<pair<int, int>> moth_ID;

    ex1::analyzer2p2 a;
    while (reader.Get_event_GM(p, moth_ID) && neve < max_ev) {
        ++neve;

        a.analz(p, moth_ID);
    }

    for (auto i : p)
        delete i.second;

    return 0;
}
