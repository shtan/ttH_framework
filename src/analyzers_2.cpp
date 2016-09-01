#ifndef analyzers_2_cpp
#define analyzers_2_cpp

#include "analyzers_2.h"
#include <iostream>
#include <cmath>

#define PI acos(-1.0)

namespace ex1
{
/* 
 * analyzer2
 */
analyzer2::analyzer2()
{
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

analyzer2::~analyzer2()
{
    for (auto i : ps)
        delete i.second;
}

/*
 * analyzer2::analz() block
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
                           double d1[4], double d2[4], double &met_px, double &met_py)
{
    for (const auto id_p : t) {
        const auto &id = id_p.first;
        if (id == 5 || id == -5) {
            CartTXYZ_to_cyl(id_p.second, b);
            met_px -= id_p.second[1];
            met_py -= id_p.second[2];
            continue;
        }
        if (id == 11 || id == -11 || id == 13 || id == -13) {
            lepton = true;
            CartTXYZ_to_cyl(id_p.second, d1);
            met_px -= id_p.second[1];
            met_py -= id_p.second[2];
            continue;
        }
        /*if (id == 12 || id == -12 || id == 14 || id == -14) {
            CartTXYZ_to_cyl(id_p.second, d2);
            continue;
        }*/
        if (id >= -4 && id <= 4 && id != 0) {
            if (d1[0] == 0) {    // check
                CartTXYZ_to_cyl(id_p.second, d1);
                met_px -= id_p.second[1];
                met_py -= id_p.second[2];
                continue;
            }
            CartTXYZ_to_cyl(id_p.second, d2);
            met_px -= id_p.second[1];
            met_py -= id_p.second[2];
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

    double met_px = 0;
    double met_py = 0;

    Fill_input_top(in.p_t1, out.t1_lep, out.b1, out.d11, out.d12, met_px, met_py);
    Fill_input_top(in.p_t2, out.t2_lep, out.b2, out.d21, out.d22, met_px, met_py);
    
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
            met_px -= id_p.second[1];
            met_py -= id_p.second[2];
            continue;
        }
        if (id == -5) {
            CartTXYZ_to_cyl(id_p.second, out.bH2);
            met_px -= id_p.second[1];
            met_py -= id_p.second[2];
            continue;
        }
        
        // "exotics" go in here
        CartTXYZ_to_cyl(id_p.second, *it1);
        met_px -= id_p.second[1];
        met_py -= id_p.second[2];
        ++it1;
    }

    out.MET_px = met_px;
    out.MET_py = met_py;

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

void analyzer2::analz(const pvec &p, const movec &moth_ID)
{
    smr.smear(p, ps);
    ttbarX ttX = ident.identify(ex1::ttH_SL_bx, p, moth_ID);
    ttbarX ttX_smr = ident.identify(ex1::ttH_SL_bx, ps, moth_ID);
    // ttX.Print_contents();
    
    Fill_input(ttX, generated.p);
    Fill_input(ttX_smr, in_2_RC.p);

    Fill_SDs(in_2_RC);
cout <<"blah" << endl;
    //recoc::Print(in_2_RC);
    recoc::output result = reco_C.reco(in_2_RC, params);
cout <<"blahh in_2_RC" << endl;
    //recoc::Print(in_2_RC);
cout << "result" << endl;
     recoc::Print(result);
    //recoc::Print_diff(result, in_2_RC);
    cout << "DIFFERENCE BETWEEN SMEARED AND GEN" << endl;
    recoc::Print_diff(in_2_RC, generated);
    cout << "DIFFERENCE BETWEEN FITTED AND GEN" << endl;
    recoc::Print_diff(result, generated);
    
    for (auto array : result.p.p_others)
        delete array;
}
}

#endif
