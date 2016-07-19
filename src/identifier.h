#include <vector>

using namespace std;

namespace ex1
{

typedef vector<pair<int, double *>> pvec;
typedef vector<pair<int, int>> movec;

enum exp_type {
    ttH_SL_bx // leptop + bbar from tbar
};

struct ttbarX {
    //     pvec p_ttbar;
    pvec p_t1;
    pvec p_t2;
    pvec p_X;
    
    void Print_contents();
    double M_t1();
    double M_t2();
    double M_X();
    double M();
};

class top_GM
{
  public:
    int id[3]; // pb, d1, d2

    double *pb;
    double *d1;
    double *d2; // neutrino

    top_GM() : id{0}
    {
        pb = nullptr;
        d1 = nullptr;
        d2 = nullptr;
    };
};

inline void top_GM_to_ttbarX(const top_GM &, const top_GM &, ttbarX &);

class identifier_GM
{
  public:
    ttbarX identify(const exp_type, const pvec &p, const movec &mo);

  private:
    inline ttbarX id_ttH_SL_bx(const pvec &p, const movec &mo);
    inline ttbarX id_H_bbx(int *id, double **ptr, pair<int, int> *mo,
                           const unsigned int sz);
    inline top_GM id_ttH_SL_bx_lep(int *id, double **ptr, pair<int, int> *mo,
                                   const unsigned int sz);
    inline void id_ttH_SL_bx_tlep(int *id, double **ptr, pair<int, int> *mo,
                                  const unsigned int sz, top_GM &t);
    inline void id_ttH_SL_bx_thad(int *id, double **ptr, pair<int, int> *mo,
                                  const unsigned int sz, top_GM &t);
};
}