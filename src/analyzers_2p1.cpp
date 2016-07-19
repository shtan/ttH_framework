#ifndef analyzers_2p1_cpp
#define analyzers_2p1_cpp

#include "analyzers_2p1.h"
#include <iostream>
#include <cmath>

namespace ex1
{

/* 
 * analyzer2p1. For single-leptonic case
 */
analyzer2p1::analyzer2p1()
{
    mes = new mekd::MEKD(Make_desc());
    My_flags(mes->flag);

    in.id = &event_id;
    in.p = &event_p;

    event_id.resize(4, 0);
    event_p.resize(4, nullptr);
    for (auto i = event_p.begin(); i != event_p.end(); ++i)
        (*i) = new double[4];
    
    // adjust your smearer's resolution
    smr.j_res.pT = 0.1;
    smr.j_res.phi = 0.01;
    smr.j_res.eta = 0.01;

    smr.l_res.pT = 0.01;
    smr.l_res.phi = 0.001;
    smr.l_res.eta = 0.001;
}

analyzer2p1::~analyzer2p1()
{
    for (auto i : ps)
        delete i.second;
    for (auto i : event_p)
        delete i;

    delete mes;
}

inline void Replace_nu_wpT_(const pvec &pv, double p[4])
{
    for (auto pi : pv) {
        if (pi.first == 12 || pi.first == -12 || pi.first == -14
            || pi.first == -14)
            continue;
        p[0] += pi.second[0];
        p[1] += pi.second[1];
        p[2] += pi.second[2];
        p[3] += pi.second[3];
    }
}

inline void Replace_nu_wpT_2(const double p_MET[4], pvec &pv)
{
    for (auto pi : pv) {
        if (pi.first == 12 || pi.first == -12 || pi.first == -14
            || pi.first == -14) {
            pi.second[0] = p_MET[0];
            pi.second[1] = p_MET[1];
            pi.second[2] = p_MET[2];
            pi.second[3] = p_MET[3];
        }
    }
}

inline void Replace_nu_wpT(ttbarX &ttH)
{
    double p[4] = {0};
    Replace_nu_wpT_(ttH.p_t1, p);
    Replace_nu_wpT_(ttH.p_t2, p);
    Replace_nu_wpT_(ttH.p_X, p);
    
    const double MET_x = -p[1];
    const double MET_y = -p[2];
    const double MET_E = sqrt(MET_x * MET_x + MET_y * MET_y);
    const double p_MET[4] = {MET_E, MET_x, MET_y, 0};
    
    Replace_nu_wpT_2(p_MET, ttH.p_t1);
    Replace_nu_wpT_2(p_MET, ttH.p_t2);
}

void analyzer2p1::analz(const pvec &p, const movec &moth_ID)
{
    smr.smear(p, ps);
    ttbarX ttX = ident.identify(ex1::ttH_SL_bx, p, moth_ID);
    ttbarX ttX_smr = ident.identify(ex1::ttH_SL_bx, ps, moth_ID);
    // ttX_smr.Print_contents();
    Replace_nu_wpT(ttX_smr);    // smeared and MET substituted
    // ttX_smr.Print_contents();

    eval_MEs_ttbb(ttX, MEs);
    const double lgME_org[2] = {log10(MEs[0]), log10(MEs[1])};
    
    eval_MEs_ttbb(ttX_smr, MEs);
    const double lgME_smr[2] = {log10(MEs[0]), log10(MEs[1])};
    
    const double d_ME_sig = lgME_smr[0] - lgME_org[0];
    const double d_ME_bkg = lgME_smr[1] - lgME_org[1];
    const double d_m_sys = ttX_smr.M() - ttX.M();
    const double d_m_t1 = ttX_smr.M_t1() - ttX.M_t1();
    const double d_m_t2 = ttX_smr.M_t2() - ttX.M_t2();
    const double d_m_X = ttX_smr.M_X() - ttX.M_X();
    cout << lgME_org[0] << " " << lgME_org[1] << " "
         << d_ME_sig << " " << d_ME_bkg << " " << d_m_sys << " "
         << d_m_t1 << " " << d_m_t2 << " " << d_m_X
         << endl;
}

inline void Fill_4momenta(const pvec &pv, double *p)
{
    p[0] = 0;
    p[1] = 0;
    p[2] = 0;
    p[3] = 0;
    
    for (auto pi : pv) {
        p[0] += pi.second[0];
        p[1] += pi.second[1];
        p[2] += pi.second[2];
        p[3] += pi.second[3];
    }
}

inline void analyzer2p1::eval_MEs_ttbb(const ttbarX &ttH, vector<double> &MEs)
{
    // inside "in" object
    event_id[0] = 6;
    event_id[1] = -6;
    event_id[2] = 5;
    event_id[3] = -5;
    
    bool t1 = false;
    for (auto p : ttH.p_t1)
        if (p.first == 5)
            t1 = true;
    if (t1) {
        Fill_4momenta(ttH.p_t1, event_p[0]);
        Fill_4momenta(ttH.p_t2, event_p[1]);
    } else {
        Fill_4momenta(ttH.p_t1, event_p[1]);
        Fill_4momenta(ttH.p_t2, event_p[0]);
    }
    
    for (auto p : ttH.p_X) {
        if (p.first == 5) {
            event_p[2][0] = p.second[0];
            event_p[2][1] = p.second[1];
            event_p[2][2] = p.second[2];
            event_p[2][3] = p.second[3];
            continue;
        }
        if (p.first == -5) {
            event_p[3][0] = p.second[0];
            event_p[3][1] = p.second[1];
            event_p[3][2] = p.second[2];
            event_p[3][3] = p.second[3];
            continue;
        }
    }

    mes->eval_MEs(in, MEs);
}
}

#endif