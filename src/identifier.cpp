#ifndef identifier_cpp
#define identifier_cpp

#include "identifier.h"
#include <cmath>
#include <iostream>

using namespace std;

namespace ex1
{

inline double Invariant_mass(double p[4])
{
    return sqrt(p[0] * p[0] - p[1] * p[1] - p[2] * p[2] - p[3] * p[3]);
}

inline void Add_p_sum(const pvec &pv, double p_out[4])
{
    for (auto p : pv) {
        p_out[0] += p.second[0];
        p_out[1] += p.second[1];
        p_out[2] += p.second[2];
        p_out[3] += p.second[3];
    }
}

/*
 * ttbarX container
 */
void ttbarX::Print_contents()
{
    cout << "p_t1:\n";
    for (auto p : p_t1)
        cout << p.first << " " << p.second[0] << " " << p.second[1] << " "
             << p.second[2] << " " << p.second[3] << endl;

    cout << "p_t2:\n";
    for (auto p : p_t2)
        cout << p.first << " " << p.second[0] << " " << p.second[1] << " "
             << p.second[2] << " " << p.second[3] << endl;

    cout << "p_X:\n";
    for (auto p : p_X)
        cout << p.first << " " << p.second[0] << " " << p.second[1] << " "
             << p.second[2] << " " << p.second[3] << endl;
    
    cout << endl;
}

double ttbarX::M_t1()
{
    double p_tot[4] = {0};
    Add_p_sum(p_t1, p_tot);
    
    return Invariant_mass(p_tot);
}

double ttbarX::M_t2()
{
    double p_tot[4] = {0};
    Add_p_sum(p_t2, p_tot);
    
    return Invariant_mass(p_tot);
}

double ttbarX::M_X()
{
    double p_tot[4] = {0};
    Add_p_sum(p_X, p_tot);
    
    return Invariant_mass(p_tot);
}

double ttbarX::M()
{
    double p_tot[4] = {0};
    Add_p_sum(p_t1, p_tot);
    Add_p_sum(p_t2, p_tot);
    Add_p_sum(p_X, p_tot);
    
    return Invariant_mass(p_tot);
}

/*
 * Helpers
 */
inline void top_GM_to_ttbarX(const top_GM &t1, const top_GM &t2, ttbarX &X)
{
    if (t1.d2 != nullptr) {
        X.p_t1.resize(3);
        X.p_t1[0].first = t1.id[0];
        X.p_t1[0].second = t1.pb;

        X.p_t1[1].first = t1.id[1];
        X.p_t1[1].second = t1.d1;

        X.p_t1[2].first = t1.id[2];
        X.p_t1[2].second = t1.d2;
    } else {
        X.p_t1.resize(2);
        X.p_t1[0].first = t1.id[0];
        X.p_t1[0].second = t1.pb;

        X.p_t1[1].first = t1.id[1];
        X.p_t1[1].second = t1.d1;
    }

    if (t1.d2 != nullptr) {
        X.p_t2.resize(3);
        X.p_t2[0].first = t2.id[0];
        X.p_t2[0].second = t2.pb;

        X.p_t2[1].first = t2.id[1];
        X.p_t2[1].second = t2.d1;

        X.p_t2[2].first = t2.id[2];
        X.p_t2[2].second = t2.d2;
    } else {
        X.p_t2.resize(2);
        X.p_t2[0].first = t2.id[0];
        X.p_t2[0].second = t2.pb;

        X.p_t2[1].first = t2.id[1];
        X.p_t2[1].second = t2.d1;
    }
}

/*
 * identifier_GM
 */
ttbarX identifier_GM::identify(const exp_type type, const pvec &p,
                               const movec &mo)
{
    if (type == ttH_SL_bx)
        return id_ttH_SL_bx(p, mo);

    return ttbarX();
}

inline ttbarX identifier_GM::id_ttH_SL_bx(const pvec &p, const movec &mo)
{
    // local links to manipulate with
    const auto sz = p.size();
    double **ptr = new double *[sz];
    pair<int, int> *amo = new pair<int, int>[sz];
    int *id = new int[sz];

    for (unsigned int i = 0; i < sz; ++i) {
        id[i] = p[i].first;
        ptr[i] = p[i].second;
        amo[i] = mo[i];
    }

    // start with Higgs
    ttbarX result = id_H_bbx(id, ptr, amo, sz);

    // start leptonic top
    top_GM t1 = id_ttH_SL_bx_lep(id, ptr, amo, sz);

    // identify leptonic comps
    id_ttH_SL_bx_tlep(id, ptr, amo, sz, t1);

    // second top
    top_GM t2;
    id_ttH_SL_bx_thad(id, ptr, amo, sz, t2);

    //     // leftovers
    //     for (unsigned int i = 0; i < sz; ++i)
    //         cout << id[i] << " " << ptr[i] << " " << amo[i].first << endl;

    // place tops into container
    top_GM_to_ttbarX(t1, t2, result);

    return result;
}

inline ttbarX identifier_GM::id_H_bbx(int *id, double **ptr, pair<int, int> *mo,
                                      const unsigned int sz)
{
    ttbarX result;
    result.p_X.reserve(2);
    for (unsigned int i = 0; i < sz; ++i) {
        if (id[i] == 5 || id[i] == -5) {
            if (mo[i].first != 25)
                continue;

            result.p_X.push_back(pair<int, double *>(id[i], ptr[i]));
            id[i] = 0;
            ptr[i] = nullptr;
        }
    }

    return result;
}

inline top_GM identifier_GM::id_ttH_SL_bx_lep(int *id, double **ptr,
                                              pair<int, int> *mo,
                                              const unsigned int sz)
{
    top_GM t;
    for (unsigned int i = 0; i < sz; ++i) {
        if (id[i] == -11 || id[i] == -13) {
            if (mo[i].first != 24)
                continue;

            t.d1 = ptr[i];
            t.id[1] = id[i];
            id[i] = 0;
            ptr[i] = nullptr;

            break;
        }
        if (id[i] == 11 || id[i] == 13) {
            if (mo[i].first != -24)
                continue;

            t.d1 = ptr[i];
            t.id[1] = id[i];
            id[i] = 0;
            ptr[i] = nullptr;

            break;
        }
    }

    return t;
}

inline void identifier_GM::id_ttH_SL_bx_tlep(int *id, double **ptr,
                                             pair<int, int> *mo,
                                             const unsigned int sz, top_GM &t)
{
    if (t.id[1] == -11 || t.id[1] == -13) {
        for (unsigned int i = 0; i < sz; ++i) {
            if (id[i] == 5) {
                if (mo[i].first != 6)
                    continue;
                t.pb = ptr[i];
                t.id[0] = id[i];
                id[i] = 0;
                ptr[i] = nullptr;
                continue;
            }

            if (id[i] == 12 || id[i] == 14) {
                if (mo[i].first != 24)
                    continue;
                t.d2 = ptr[i];
                t.id[2] = id[i];
                id[i] = 0;
                ptr[i] = nullptr;
                continue;
            }
        }
        return;
    }

    if (t.id[1] == 11 || t.id[1] == 13) {
        for (unsigned int i = 0; i < sz; ++i) {
            if (id[i] == -5) {
                if (mo[i].first != -6)
                    continue;
                t.pb = ptr[i];
                t.id[0] = id[i];
                id[i] = 0;
                ptr[i] = nullptr;
                continue;
            }

            if (id[i] == -12 || id[i] == -14) {
                if (mo[i].first != -24)
                    continue;
                t.d2 = ptr[i];
                t.id[2] = id[i];
                id[i] = 0;
                ptr[i] = nullptr;
                continue;
            }
        }
        return;
    }
}

inline void identifier_GM::id_ttH_SL_bx_thad(int *id, double **ptr,
                                             pair<int, int> *mo,
                                             const unsigned int sz, top_GM &t)
{
    // tbar part
    for (unsigned int i = 0; i < sz; ++i) {
        if (id[i] == -5) {
            if (mo[i].first != -6)
                continue;

            t.pb = ptr[i];
            t.id[0] = id[i];
            id[i] = 0;
            ptr[i] = nullptr;
            continue;
        }

        if (id[i] <= 4 && id[i] >= -4 && id[i] != 0) {
            if (mo[i].first != -24)
                continue;
            if (t.d1 == nullptr) {
                t.d1 = ptr[i];
                t.id[1] = id[i];
            } else {
                t.d2 = ptr[i];
                t.id[2] = id[i];
            }
            id[i] = 0;
            ptr[i] = nullptr;
        }
    }
    if (t.pb != nullptr && t.d1 != nullptr)
        return;

    // t part
    for (unsigned int i = 0; i < sz; ++i) {
        if (id[i] == 5) {
            if (mo[i].first != 6)
                continue;

            t.pb = ptr[i];
            t.id[0] = id[i];
            id[i] = 0;
            ptr[i] = nullptr;
            continue;
        }

        if (id[i] <= 4 && id[i] >= -4 && id[i] != 0) {
            if (mo[i].first != 24)
                continue;
            if (t.d1 == nullptr) {
                t.d1 = ptr[i];
                t.id[1] = id[i];
            } else {
                t.d2 = ptr[i];
                t.id[2] = id[i];
            }
            id[i] = 0;
            ptr[i] = nullptr;
        }
    }
}
}

#endif