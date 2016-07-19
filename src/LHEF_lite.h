#include <iostream>
#include <vector>
#include <algorithm> // for sort

using namespace std;

namespace LHEF
{
class Reader; // forward declaration
}

class LHEF_lite
{
  public:
    LHEF_lite(const string &lhe);
    ~LHEF_lite();

    // closer to experiment with a perfect particle ID
    bool Get_event(vector<pair<int, double *>> &p);
    bool Get_event_sorted(vector<pair<int, double *>> &p);

    // God-mode methods
    bool Get_event_GM(vector<pair<int, double *>> &p,
                      vector<pair<int, int>> &mth_IDs);

  private:
    LHEF::Reader *reader;

    inline void Get_event_(const LHEF::Reader &r,
                           vector<pair<int, double *>> &p, const bool sort);
    inline void Get_event_GM_(const LHEF::Reader &r,
                              vector<pair<int, double *>> &p,
                              vector<pair<int, int>> &mth_IDs);
    inline unsigned int Get_event_sz_prescan(const vector<int> &status);

    inline void Get_event_sorted_(vector<pair<int, double *>> &p);
    inline void resize(vector<pair<int, double *>> &p, const unsigned int sz);
};
