#ifndef smearer_cpp
#define smearer_cpp

#include "smearer.h"
#include <iostream>

#define PI acos(-1.0)

smearer::smearer() : normal(0, 1)
{
    j_res.pT = 0.1;
    j_res.phi = 0.01;
    j_res.eta = 0.01;

    l_res.pT = 0.01;
    l_res.phi = 0.001;
    l_res.eta = 0.001;

    j_max_id = 5;
}

void smearer::smear(const vector<pair<int, double *>> &p,
                    vector<pair<int, double *>> &op)
{
    resize(op, p.size());
    auto it1 = p.begin();
    auto it2 = op.begin();

    while (it1 != p.end()) {
        const auto id = it1->first;
        const auto in_p = it1->second;
        const auto out_p = it2->second;
        it2->first = id;
        ++it1;
        ++it2;

        if (id <= j_max_id && id >= -j_max_id) {
            smear_(j_res, in_p, 0.003, out_p);
            continue;
        }
        if (id == 11 || id == -11) {
            smear_(l_res, in_p, 0.000511, out_p);
            continue;
        }
        if (id == 13 || id == -13) {
            smear_(l_res, in_p, 0.101, out_p);
            continue;
        }
        out_p[0] = in_p[0];
        out_p[1] = in_p[1];
        out_p[2] = in_p[2];
        out_p[3] = in_p[3];
    }
}

inline void smearer::smear_(const resolution &res, const double *p,
                            const double m, double *op)
{
    const double pT = sqrt(p[1] * p[1] + p[2] * p[2]);
    const double eta = asinh(p[3] / pT);
    double phi = asin(p[2] / pT);
    if (p[1] < 0)
        phi = PI - phi;

    const double rel_pT_dev = 1 + res.pT * normal(generator);
    const double pT_n = pT * rel_pT_dev;
    const double phi_n = phi + res.phi * normal(generator);
    const double eta_n = eta + res.eta * normal(generator);

    const double px = cos(phi_n) * pT_n;
    const double py = sin(phi_n) * pT_n;
    const double pz = sinh(eta_n) * pT_n;

    op[0] = sqrt(m * m + px * px + py * py + pz * pz);
    op[1] = px;
    op[2] = py;
    op[3] = pz;
}

inline void smearer::resize(vector<pair<int, double *>> &p,
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

#endif