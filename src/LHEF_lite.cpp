#ifndef LHEF_lite_cpp
#define LHEF_lite_cpp

#include "LHEF_lite.h"
#include "LHEF/LHEF.h"

LHEF_lite::LHEF_lite(const string &lhe) { reader = new LHEF::Reader(lhe); }

LHEF_lite::~LHEF_lite() { delete reader; }

bool LHEF_lite::Get_event(vector<pair<int, double *>> &p)
{
    if (reader->readEvent()) {
        Get_event_(*reader, p, false);

        return true;
    }

    return false;
}

bool LHEF_lite::Get_event_sorted(vector<pair<int, double *>> &p)
{
    if (reader->readEvent()) {
        Get_event_(*reader, p, true);

        return true;
    }

    return false;
}

bool LHEF_lite::Get_event_GM(vector<pair<int, double *>> &p,
                             vector<pair<int, int>> &mth_IDs)
{
    if (reader->readEvent()) {
        Get_event_GM_(*reader, p, mth_IDs);

        return true;
    }

    return false;
}

inline void LHEF_lite::Get_event_(const LHEF::Reader &r,
                                  vector<pair<int, double *>> &p,
                                  const bool sort)
{
    const auto ev = &r.hepeup;
    // The PDG id's for the particle entries in this event.
    const auto IDs = &(ev->IDUP);
    // The status codes for the particle entries in this event.
    const auto sts = &(ev->ISTUP);
    // Lab frame momentum (Px, Py, Pz, E and M in GeV) for the particle entries
    // in this event. vector<vector<double>>
    const auto ps = &(ev->PUP);

    const auto sz = IDs->size();
    if (sz != sts->size() || sz != ps->size()) {
        cerr << "Problems with vector sizes\n";
        return;
    }

    // relevant size prescan for memory reallocation
    const unsigned int sz_ = Get_event_sz_prescan(*sts);
    resize(p, sz_);

    // fill vector with final particles
    {
        unsigned int i1 = 0;
        unsigned int i2 = 0;
        while (i1 < sz && i2 < sz_) {
            if ((*sts)[i1] == 1) {
                p[i2].first = (*IDs)[i1];
                p[i2].second[0] = (*ps)[i1][3];
                p[i2].second[1] = (*ps)[i1][0];
                p[i2].second[2] = (*ps)[i1][1];
                p[i2].second[3] = (*ps)[i1][2];
                ++i2;
            }
            ++i1;
        }
    }

    if (!sort)
        return;

    Get_event_sorted_(p);
}

inline void LHEF_lite::Get_event_GM_(const LHEF::Reader &r,
                                     vector<pair<int, double *>> &p,
                                     vector<pair<int, int>> &mth_IDs)
{
    const auto ev = &r.hepeup;
    // The PDG id's for the particle entries in this event.
    const auto IDs = &(ev->IDUP);
    // The status codes for the particle entries in this event.
    const auto sts = &(ev->ISTUP);
    // Lab frame momentum (Px, Py, Pz, E and M in GeV) for the particle entries
    // in this event. vector<vector<double>>
    const auto ps = &(ev->PUP);
    const auto mth = &(ev->MOTHUP);

    const auto sz = IDs->size();
    if (sz != sts->size() || sz != ps->size()) {
        cerr << "Problems with vector sizes\n";
        return;
    }

    // relevant size prescan for memory reallocation
    const unsigned int sz_ = Get_event_sz_prescan(*sts);
    resize(p, sz_);
    mth_IDs.resize(sz_);

    // fill vector with final particles
    {
        unsigned int i1 = 0;
        unsigned int i2 = 0;
        while (i1 < sz && i2 < sz_) {
            if ((*sts)[i1] == 1) {
                p[i2].first = (*IDs)[i1];
                mth_IDs[i2].first = (*IDs)[(*mth)[i1].first - 1];
                mth_IDs[i2].second = (*IDs)[(*mth)[i1].second - 1];
                p[i2].second[0] = (*ps)[i1][3];
                p[i2].second[1] = (*ps)[i1][0];
                p[i2].second[2] = (*ps)[i1][1];
                p[i2].second[3] = (*ps)[i1][2];
                ++i2;
            }
            ++i1;
        }
    }
}

inline unsigned int LHEF_lite::Get_event_sz_prescan(const vector<int> &status)
{
    unsigned int sz_ = 0;
    for (auto i : status) {
        if (i == 1)
            ++sz_;
    }

    return sz_;
}

inline void LHEF_lite::resize(vector<pair<int, double *>> &p,
                              const unsigned int sz)
{
    const auto old_sz = p.size();
    if (old_sz < sz) {
        p.reserve(sz);
        for (auto i = old_sz; i < sz; ++i)
            p.push_back(pair<int, double *>(0, new double[4]));

        return;
    }
    if (old_sz > sz) {
        for (auto i = (old_sz - 1); i >= sz; --i)
            delete p[i].second;
        p.resize(sz);

        return;
    }
}

inline void LHEF_lite::Get_event_sorted_(vector<pair<int, double *>> &p)
{
    // sort by absolute IDs + 0.5 for neg ID
    // thus sort: 1 1 -1 -1 2 2 -2 -2, etc.
    typedef pair<int, double *> cell_t;
    vector<pair<float, cell_t>> vmap;
    vmap.reserve(p.size());
    for (auto i : p) {
        float abs_id;
        const auto id = i.first;
        if (id > 0)
            abs_id = id;
        else
            abs_id = -id + 0.5;

        vmap.push_back(pair<float, cell_t>(abs_id, i));
    }
    sort(vmap.begin(), vmap.end());

    {
        auto ip = p.begin();
        auto ivmap = vmap.begin();
        while (ip != p.end()) {
            *ip = (*ivmap).second;
            ++ip;
            ++ivmap;
        }
    }
}

#endif