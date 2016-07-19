#ifndef analyzers_1_cpp
#define analyzers_1_cpp

#include "analyzers_1.h"
#include <iostream>
#include <cmath>

namespace ex1
{

vector<mekd::process_description> Make_desc()
{
    mekd::process_description sig;
    mekd::process_description bkg;

    sig.model = mekd::model_HEFTU;
    sig.process = mekd::proc_ttX;
    sig.resonance = mekd::reson_Spin0Pm;
    sig.production = mekd::prod_all;
    sig.decay = mekd::decay_2f;
    sig.final_state = mekd::final_ttbb;

    bkg.model = mekd::model_HEFTU;
    bkg.process = mekd::proc_ttX;
    bkg.resonance = mekd::bkg_tt;
    bkg.production = mekd::prod_all;
    bkg.decay = mekd::decay_2f;
    bkg.final_state = mekd::final_ttbb;

    vector<mekd::process_description> init_desc;
    init_desc.reserve(2);
    init_desc.push_back(sig);
    init_desc.push_back(bkg);

    return init_desc;
}

void My_flags(mekd::flags &flag)
{
    flag.Use_PDF_w_pT0 = true;
    flag.use_prod_ddx = true;
    flag.use_prod_uux = true;
    flag.use_prod_ssx = true;
    flag.use_prod_ccx = true;
    flag.use_prod_bbx = false;
    flag.use_prod_gg = true;
}

/* 
 * analyzer1
 */
analyzer1::analyzer1()
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

analyzer1::~analyzer1()
{
    for (auto i : ps)
        delete i.second;
    for (auto i : event_p)
        delete i;

    delete mes;
}

void analyzer1::analz(const pvec &p, const movec &moth_ID)
{
    smr.smear(p, ps);
    ttbarX ttX = ident.identify(ex1::ttH_SL_bx, p, moth_ID);
    ttbarX ttX_smr = ident.identify(ex1::ttH_SL_bx, ps, moth_ID);
    // ttX.Print_contents();

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

// hereustic particle numbering
inline void analyzer1::eval_MEs_ttbb_bjjbbarlnbbar(const pvec &p)
{
    // inside "in" object
    event_id[0] = 6;
    event_id[1] = -6;
    event_id[2] = 5;
    event_id[3] = -5;
    for (unsigned int i = 0; i < 4; ++i) {
        event_p[0][i] = p[0].second[i] + p[1].second[i] + p[2].second[i];
        event_p[1][i] = p[3].second[i] + p[4].second[i] + p[5].second[i];
        event_p[2][i] = p[6].second[i];
        event_p[3][i] = p[7].second[i];
    }

    mes->eval_MEs(in, MEs);
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

inline void analyzer1::eval_MEs_ttbb(const ttbarX &ttH, vector<double> &MEs)
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