#include "analyzers.h"
#include "smearer.h"
#include "identifier.h"
#include "MEKD.h"

using namespace std;

namespace ex1
{

vector<mekd::process_description> Make_desc();
void My_flags(mekd::flags &);

class analyzer1
{
  public:
    analyzer1();
    ~analyzer1();

    void analz(const pvec &, const movec &);
    inline void eval_MEs_ttbb_bjjbbarlnbbar(const pvec &);
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