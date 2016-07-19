#include "analyzers.h"
#include "analyzers_1.h"
// #include "smearer.h" // in "analyzers_1.h"
// #include "identifier.h"  // in "analyzers_1.h"
// #include "MEKD.h"    // in "analyzers_1.h"
// #include "ttHReco/ttbar_reco_C.h"

using namespace std;

namespace ex1
{

class analyzer2p1
{
  public:
    analyzer2p1();
    ~analyzer2p1();

    void analz(const pvec &, const movec &);
    inline void eval_MEs_ttbb(const ttbarX &ttH, vector<double> &MEs);

    smearer smr;
    vector<pair<int, double *>> ps; // smeared particles
    ex1::identifier_GM ident;

  private:
    mekd::MEKD *mes;
    mekd::input in;
    vector<int> event_id;
    vector<double *> event_p;
    vector<double> MEs;
};
}