#ifndef analyzers_2p2_cpp
#define analyzers_2p2_cpp

#include "analyzers_2p2.h"
#include <iostream>
#include <cmath>

#define PI acos(-1.0)

namespace ex1
{
/* 
 * analyzer2p2. For single-leptonic case
 */
analyzer2p2::analyzer2p2()
{
    mes = new mekd::MEKD(Make_desc());
    My_flags(mes->flag);

    in.id = &event_id;
    in.p = &event_p;

    event_id.resize(4, 0);
    event_p.resize(4, nullptr);
    for (auto i = event_p.begin(); i != event_p.end(); ++i)
        (*i) = new double[4];
    
    temp.p_t1.resize(3);
    temp.p_t2.resize(3);
    temp.p_X.resize(2);
    for (auto i = temp.p_t1.begin(); i != temp.p_t1.end(); ++i)
        (*i) = pair<int, double *>(0, new double[4]);
    for (auto i = temp.p_t2.begin(); i != temp.p_t2.end(); ++i)
        (*i) = pair<int, double *>(0, new double[4]);
    for (auto i = temp.p_X.begin(); i != temp.p_X.end(); ++i)
        (*i) = pair<int, double *>(0, new double[4]);
    
    // adjust your smearer's resolution
    smr.j_res.pT = 0.1;
    smr.j_res.phi = 0.01;
    smr.j_res.eta = 0.01;

    smr.l_res.pT = 0.01;
    smr.l_res.phi = 0.001;
    smr.l_res.eta = 0.001;
    
    params.m_t = 173;
    params.m_W = 80.4;
    params.SD_m_t = 2.0;
    params.SD_m_W = 2.085;
}

analyzer2p2::~analyzer2p2()
{
    for (auto i : ps)
        delete i.second;
    for (auto i : event_p)
        delete i;
    
    for (auto i : temp.p_t1)
        delete i.second;
    for (auto i : temp.p_t2)
        delete i.second;
    for (auto i : temp.p_X)
        delete i.second;

    delete mes;
}

/*
 * analyzer2p2::analz() block
 * 
 */
inline void CartTXYZ_to_cyl(const double *pTXYZ, double *pcyl)
{
    const double pT = sqrt(pTXYZ[1] * pTXYZ[1] + pTXYZ[2] * pTXYZ[2]);
    const double eta = asinh(pTXYZ[3] / pT);
    double phi = asin(pTXYZ[2] / pT);
    if (pTXYZ[1] < 0)
        phi = PI - phi;
    pcyl[3] = pTXYZ[0];
    pcyl[0] = pT;
    pcyl[1] = phi;
    pcyl[2] = eta;
}

inline void P4_resize(vector<double *> &p, const unsigned int sz)
{
    const auto old_sz = p.size();
    if (old_sz < sz) {
        p.reserve(sz);
        for (auto i = old_sz; i < sz; ++i)
            p.push_back(new double[4]);

        return;
    }
    if (old_sz > sz) {
        for (auto i = (old_sz - 1); i >= sz; --i)
            delete p[i];
        p.resize(sz);

        return;
    }
}

inline void Fill_input_top(const pvec &t, bool &lepton, double b[4],
                           double d1[4], double d2[4])
{
    for (const auto id_p : t) {
        const auto &id = id_p.first;
        if (id == 5 || id == -5) {
            CartTXYZ_to_cyl(id_p.second, b);
            continue;
        }
        if (id == 11 || id == -11 || id == 13 || id == -13) {
            lepton = true;
            CartTXYZ_to_cyl(id_p.second, d1);
            continue;
        }
        if (id >= -4 && id <= 4) {
            if (d1[0] == 0) {    // check
                CartTXYZ_to_cyl(id_p.second, d1);
                continue;
            }
            CartTXYZ_to_cyl(id_p.second, d2);
            continue;
        }
    }
}

inline void Fill_input(const ttbarX &in, recoc::input_x<4> &out)
{
    out.t1_lep = false;
    out.t2_lep = false;
    out.d11[0] = 0;   // check
    out.d21[0] = 0;   // check
    
    out.d12[0] = 0;
    out.d12[1] = 0;
    out.d12[2] = 0;
    out.d12[3] = 0;
    out.d22[0] = 0;
    out.d22[1] = 0;
    out.d22[2] = 0;
    out.d22[3] = 0;
    Fill_input_top(in.p_t1, out.t1_lep, out.b1, out.d11, out.d12);
    Fill_input_top(in.p_t2, out.t2_lep, out.b2, out.d21, out.d22);
    
    // link over stuff
    unsigned int sz;
    if (in.p_X.size() >= 2)
        sz = in.p_X.size() - 2;
    else
        sz = 0;
    P4_resize(out.p_others, sz);
    
    auto it1 = out.p_others.begin();
    for (const auto id_p : in.p_X) {
        const auto &id = id_p.first;
        if (id == 5) {
            CartTXYZ_to_cyl(id_p.second, out.bH1);
            continue;
        }
        if (id == -5) {
            CartTXYZ_to_cyl(id_p.second, out.bH2);
            continue;
        }
        
        // "exotics" go in here
        CartTXYZ_to_cyl(id_p.second, *it1);
        ++it1;
    }
}

inline void Fill_SDs_(const double p[4], const double unc[3], double SD[3])
{
    SD[0] = unc[0] * p[0];
    SD[1] = unc[1] * p[1];
    SD[2] = unc[2] * p[2];
}

void Fill_SDs(recoc::input &in)
{
    auto &p = in.p;
    auto &SD = in.SD;
    const double j_res[3] = {0.1, 0.01, 0.01};  // jet: pT, phi, eta
    const double l_res[3] = {0.01, 0.001, 0.001};   // lepton: pT, phi, eta
    
    Fill_SDs_(p.b1, j_res, SD.b1);
    Fill_SDs_(p.b2, j_res, SD.b2);
    Fill_SDs_(p.bH1, j_res, SD.bH1);
    Fill_SDs_(p.bH2, j_res, SD.bH2);
    if (p.t1_lep)
        Fill_SDs_(p.d11, l_res, SD.d11);
    else {
        Fill_SDs_(p.d11, j_res, SD.d11);
        Fill_SDs_(p.d12, j_res, SD.d12);
    }
    if (p.t2_lep)
        Fill_SDs_(p.d21, l_res, SD.d21);
    else {
        Fill_SDs_(p.d21, j_res, SD.d21);
        Fill_SDs_(p.d22, j_res, SD.d22);
    }
    
    if (p.p_others.size() == 0)
        return;
    
    auto it1 = p.p_others.begin();
    auto it2 = SD.p_others.begin();
    const auto end = p.p_others.end();
    while (it1 != end) {
        Fill_SDs_(*it1, j_res, *it2);
    }
}

inline void Cyl_to_cartTXYZ(const double *pcyl, double *pTXYZ)
{
    const auto &pT = pcyl[0];
    const auto &phi = pcyl[1];
    const auto &eta = pcyl[2];

    pTXYZ[0] = pcyl[3];
    pTXYZ[1] = cos(phi) * pT;
    pTXYZ[2] = sin(phi) * pT;
    pTXYZ[3] = sinh(eta) * pT;
}

inline void Copy_IDs(const ttbarX &ID_src, ttbarX &out)
{
    if (ID_src.p_t1.size() != out.p_t1.size() ||
        ID_src.p_t2.size() != out.p_t2.size() ||
        ID_src.p_X.size() != out.p_X.size()) {
        cerr << "Size mismatch. Copy_IDs failed.\n";
        return;
    }
    auto it = out.p_t1.begin();
    for (auto i : ID_src.p_t1) {
        it->first = i.first;
        ++it;
    }
    it = out.p_t2.begin();
    for (auto i : ID_src.p_t2) {
        it->first = i.first;
        ++it;
    }
    it = out.p_X.begin();
    for (auto i : ID_src.p_X) {
        it->first = i.first;
        ++it;
    }
}

inline void Transform(const recoc::output &in, ttbarX &out)
{
    Cyl_to_cartTXYZ(in.p.b1, out.p_t1[0].second);
    Cyl_to_cartTXYZ(in.p.d11, out.p_t1[1].second);
    Cyl_to_cartTXYZ(in.p.d12, out.p_t1[2].second);
    
    Cyl_to_cartTXYZ(in.p.b2, out.p_t2[0].second);
    Cyl_to_cartTXYZ(in.p.d21, out.p_t2[1].second);
    Cyl_to_cartTXYZ(in.p.d22, out.p_t2[2].second);
    
    Cyl_to_cartTXYZ(in.p.bH1, out.p_X[0].second);
    Cyl_to_cartTXYZ(in.p.bH2, out.p_X[1].second);
}

void analyzer2p2::analz(const pvec &p, const movec &moth_ID)
{
    smr.smear(p, ps);
    ttbarX ttX = ident.identify(ex1::ttH_SL_bx, p, moth_ID);
    // ttbarX ttX_smr = ident.identify(ex1::ttH_SL_bx, ps, moth_ID);
    // ttX.Print_contents();
    
    Fill_input(ttX, in_2_RC.p);
    Fill_SDs(in_2_RC);
    // recoc::Print(in_2_RC);
    recoc::output result = reco_C.reco(in_2_RC, params);
    // recoc::Print(result);
    // recoc::Print_diff(result, in_2_RC);
    Copy_IDs(ttX, temp);
    Transform(result, temp);
    // temp.Print_contents();
    
    eval_MEs_ttbb(ttX, MEs);
    const double lgME_org[2] = {log10(MEs[0]), log10(MEs[1])};
    
    eval_MEs_ttbb(temp, MEs);
    const double lgME_rec[2] = {log10(MEs[0]), log10(MEs[1])};
    
    const double d_ME_sig = lgME_rec[0] - lgME_org[0];
    const double d_ME_bkg = lgME_rec[1] - lgME_org[1];
    const double d_m_sys = temp.M() - ttX.M();
    const double d_m_t1 = temp.M_t1() - ttX.M_t1();
    const double d_m_t2 = temp.M_t2() - ttX.M_t2();
    const double d_m_X = temp.M_X() - ttX.M_X();
    cout << "RESULT: ";
    cout << lgME_org[0] << " " << lgME_org[1] << " "
         << d_ME_sig << " " << d_ME_bkg << " " << d_m_sys << " "
         << d_m_t1 << " " << d_m_t2 << " " << d_m_X
         << endl;
    
    for (auto array : result.p.p_others)
        delete array;
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

inline void analyzer2p2::eval_MEs_ttbb(const ttbarX &ttH, vector<double> &MEs)
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