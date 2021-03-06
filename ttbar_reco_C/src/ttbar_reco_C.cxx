#ifndef ttbar_reco_C_cxx
#define ttbar_reco_C_cxx

#include "ttHReco/ttbar_reco_C.h"
#include <cmath>

using namespace commonstruct;

namespace recoc
{
ttRC::ttRC() {}

/*
 * Diagnostics block
 *
 */
// input/output class/struct (IOC)
template<class IOC>
inline void Print_(const IOC &in)
{
    if (in.p.t1_lep)
        cout << "t1 leptonic:\n";
    else
        cout << "t1 hadronic:\n";

    cout << "b1: " << in.p.b1[0] << " " << in.p.b1[1] << " " << in.p.b1[2]
         << " " << in.p.b1[3] << endl;
    cout << "d11: " << in.p.d11[0] << " " << in.p.d11[1] << " " << in.p.d11[2]
         << " " << in.p.d11[3] << endl;
    cout << "d12: " << in.p.d12[0] << " " << in.p.d12[1] << " " << in.p.d12[2]
         << " " << in.p.d12[3] << endl;

    if (in.p.t2_lep)
        cout << "t2 leptonic:\n";
    else
        cout << "t2 hadronic:\n";

    cout << "b2: " << in.p.b2[0] << " " << in.p.b2[1] << " " << in.p.b2[2]
         << " " << in.p.b2[3] << endl;
    cout << "d21: " << in.p.d21[0] << " " << in.p.d21[1] << " " << in.p.d21[2]
         << " " << in.p.d21[3] << endl;
    cout << "d22: " << in.p.d22[0] << " " << in.p.d22[1] << " " << in.p.d22[2]
         << " " << in.p.d22[3] << endl;

    cout << "Other:\n";
    cout << "bH1: " << in.p.bH1[0] << " " << in.p.bH1[1] << " " << in.p.bH1[2]
         << " " << in.p.bH1[3] << endl;
    cout << "bH2: " << in.p.bH2[0] << " " << in.p.bH2[1] << " " << in.p.bH2[2]
         << " " << in.p.bH2[3] << endl;

    for (auto p : in.p.p_others)
        cout << "p_others: " << p[0] << " " << p[1] << " " << p[2] << " "
             << p[3] << endl;
    
}

void Print(const input &in)
{
    if (in.p.t1_lep)
        cout << "t1 leptonic:\n";
    else
        cout << "t1 hadronic:\n";

    cout << "b1: " << in.p.b1[0] << " " << in.p.b1[1] << " " << in.p.b1[2]
         << " " << in.p.b1[3];
    cout << " +- " << in.SD.b1[0] << " " << in.SD.b1[1] << " " << in.SD.b1[2]
         << " " << endl;
    cout << "d11: " << in.p.d11[0] << " " << in.p.d11[1] << " " << in.p.d11[2]
         << " " << in.p.d11[3];
    cout << " +- " << in.SD.d11[0] << " " << in.SD.d11[1] << " " << in.SD.d11[2]
         << " " << endl;
    cout << "d12: " << in.p.d12[0] << " " << in.p.d12[1] << " " << in.p.d12[2]
         << " " << in.p.d12[3];
    cout << " +- " << in.SD.d12[0] << " " << in.SD.d12[1] << " " << in.SD.d12[2]
         << " " << endl;

    if (in.p.t2_lep)
        cout << "t2 leptonic:\n";
    else
        cout << "t2 hadronic:\n";

    cout << "b2: " << in.p.b2[0] << " " << in.p.b2[1] << " " << in.p.b2[2]
         << " " << in.p.b2[3];
    cout << " +- " << in.SD.b2[0] << " " << in.SD.b2[1] << " " << in.SD.b2[2]
         << " " << endl;
    cout << "d21: " << in.p.d21[0] << " " << in.p.d21[1] << " " << in.p.d21[2]
         << " " << in.p.d21[3];
    cout << " +- " << in.SD.d21[0] << " " << in.SD.d21[1] << " " << in.SD.d21[2]
         << " " << endl;
    cout << "d22: " << in.p.d22[0] << " " << in.p.d22[1] << " " << in.p.d22[2]
         << " " << in.p.d22[3];
    cout << " +- " << in.SD.d22[0] << " " << in.SD.d22[1] << " " << in.SD.d22[2]
         << " " << endl;

    cout << "Other:\n";
    cout << "bH1: " << in.p.bH1[0] << " " << in.p.bH1[1] << " " << in.p.bH1[2]
         << " " << in.p.bH1[3];
    cout << " +- " << in.SD.bH1[0] << " " << in.SD.bH1[1] << " " << in.SD.bH1[2]
         << " " << endl;
    cout << "bH2: " << in.p.bH2[0] << " " << in.p.bH2[1] << " " << in.p.bH2[2]
         << " " << in.p.bH2[3];
    cout << " +- " << in.SD.bH2[0] << " " << in.SD.bH2[1] << " " << in.SD.bH2[2]
         << " " << endl;

    for (auto p : in.p.p_others)
        cout << "p_others: " << p[0] << " " << p[1] << " " << p[2] << " "
             << p[3] << endl;

    cout << in.MET_px << endl;
    cout << in.MET_py << endl;
}

void Print(const output &in)
{
    cout << "w: " << in.weight << endl;
    Print_(in);
}

void Print_diff(const output &final, const input &initial)
{
    output diff;

    diff.p.t1_lep = initial.p.t1_lep;
    diff.p.t2_lep = initial.p.t2_lep;

    for (unsigned int i = 0; i < 4; ++i) {
        diff.p.b1[i] = final.p.b1[i] - initial.p.b1[i];
        diff.p.d11[i] = final.p.d11[i] - initial.p.d11[i];
        diff.p.d12[i] = final.p.d12[i] - initial.p.d12[i];
        
        diff.p.b2[i] = final.p.b2[i] - initial.p.b2[i];
        diff.p.d21[i] = final.p.d21[i] - initial.p.d21[i];
        diff.p.d22[i] = final.p.d22[i] - initial.p.d22[i];
        
        diff.p.bH1[i] = final.p.bH1[i] - initial.p.bH1[i];
        diff.p.bH2[i] = final.p.bH2[i] - initial.p.bH2[i];
    }
    
    const auto sz = final.p.p_others.size();
    if (sz > 0) {
        diff.p.p_others.reserve(sz);
        auto itf = final.p.p_others.begin();
        auto iti = initial.p.p_others.begin();
        
        while (itf != final.p.p_others.end()) {
            double *p = new double[4];
            p[0] = (*itf)[0] - (*iti)[0];
            p[1] = (*itf)[1] - (*iti)[1];
            p[2] = (*itf)[2] - (*iti)[2];
            p[3] = (*itf)[3] - (*iti)[3];
            diff.p.p_others.push_back(p);
        }
    }
    
    Print_(diff);

    cout << "MET pt = " << pow( ( pow(initial.MET_px,2) + pow(initial.MET_py,2) ), 0.5 ) << endl;
    cout << "MET phi = " << atan2( initial.MET_py, initial.MET_px ) << endl;

}

void Print_diff(const input &final, const input &initial)
{
    output diff;

    diff.p.t1_lep = initial.p.t1_lep;
    diff.p.t2_lep = initial.p.t2_lep;

    for (unsigned int i = 0; i < 4; ++i) {
        diff.p.b1[i] = final.p.b1[i] - initial.p.b1[i];
        diff.p.d11[i] = final.p.d11[i] - initial.p.d11[i];
        diff.p.d12[i] = final.p.d12[i] - initial.p.d12[i];
        
        diff.p.b2[i] = final.p.b2[i] - initial.p.b2[i];
        diff.p.d21[i] = final.p.d21[i] - initial.p.d21[i];
        diff.p.d22[i] = final.p.d22[i] - initial.p.d22[i];
        
        diff.p.bH1[i] = final.p.bH1[i] - initial.p.bH1[i];
        diff.p.bH2[i] = final.p.bH2[i] - initial.p.bH2[i];
    }
    
    const auto sz = final.p.p_others.size();
    if (sz > 0) {
        diff.p.p_others.reserve(sz);
        auto itf = final.p.p_others.begin();
        auto iti = initial.p.p_others.begin();
        
        while (itf != final.p.p_others.end()) {
            double *p = new double[4];
            p[0] = (*itf)[0] - (*iti)[0];
            p[1] = (*itf)[1] - (*iti)[1];
            p[2] = (*itf)[2] - (*iti)[2];
            p[3] = (*itf)[3] - (*iti)[3];
            diff.p.p_others.push_back(p);
        }
    }
    
    Print_(diff);

    cout << "MET pt = " << pow( ( pow(initial.MET_px,2) + pow(initial.MET_py,2) ), 0.5 ) << endl;
    cout << "MET phi = " << atan2( initial.MET_py, initial.MET_px ) << endl;

}

/*
 * ttRC::reco() block
 *
 */
inline double Make_px(const double &pT, const double &phi)
{
    return cos(phi) * pT;
}

inline double Make_py(const double &pT, const double &phi)
{
    return sin(phi) * pT;
}

inline double Make_pz(const double &pT, const double &eta)
{
    return sinh(eta) * pT;
}

inline double Make_mass(const double &pT, const double &phi, const double &eta, const double &energy)
{
    LV vec;
    vec.SetPtEtaPhiE(pT, eta, phi, energy);
    return vec.M();
}

inline double Make_pT(const double &px, const double &py)
{
    return pow( (pow(px,2) + pow(py,2)), 0.5 );
}

inline double Make_phi(const double &px, const double &py)
{
    return atan2(py, px);
}

inline void Eval_p_sum_tt_cartesian(const input_x<4> &inp, double outp[3])
{
    double px_tot = 0;
    double py_tot = 0;
    double pz_tot = 0;

    px_tot += Make_px(inp.b1[0], inp.b1[1]);
    py_tot += Make_py(inp.b1[0], inp.b1[1]);
    pz_tot += Make_pz(inp.b1[0], inp.b1[2]);

    px_tot += Make_px(inp.b2[0], inp.b2[1]);
    py_tot += Make_py(inp.b2[0], inp.b2[1]);
    pz_tot += Make_pz(inp.b2[0], inp.b2[2]);

    px_tot += Make_px(inp.d11[0], inp.d11[1]);
    py_tot += Make_py(inp.d11[0], inp.d11[1]);
    pz_tot += Make_pz(inp.d11[0], inp.d11[2]);

    px_tot += Make_px(inp.d12[0], inp.d12[1]);
    py_tot += Make_py(inp.d12[0], inp.d12[1]);
    pz_tot += Make_pz(inp.d12[0], inp.d12[2]);

    if (!inp.t1_lep) {
        px_tot += Make_px(inp.d12[0], inp.d12[1]);
        py_tot += Make_py(inp.d12[0], inp.d12[1]);
        pz_tot += Make_pz(inp.d12[0], inp.d12[2]);
    }
    if (!inp.t2_lep) {
        px_tot += Make_px(inp.d22[0], inp.d22[1]);
        py_tot += Make_py(inp.d22[0], inp.d22[1]);
        pz_tot += Make_pz(inp.d22[0], inp.d22[2]);
    }

    outp[0] = px_tot;
    outp[1] = py_tot;
    outp[2] = pz_tot;
}

inline void Insert_cartesian_XYZT(const double p[4],
                                  vector<LV>::iterator &it_LV)
{
    const double px = Make_px(p[0], p[1]);
    const double py = Make_py(p[0], p[1]);
    const double pz = Make_pz(p[0], p[2]);
    it_LV->SetPxPyPzE(px, py, pz, p[3]);
}


inline void Reset_bigstruct(big_struct &bigstruct)
{
    for (int iTop = 0; iTop < bigstruct.n_tops(); ++iTop){
        delete bigstruct.tops.at(iTop);
        bigstruct.tops.at(iTop) = NULL;
    }
    delete bigstruct.nontops_ptr;
    bigstruct.nontops_ptr = NULL;

    bigstruct.tops.clear();
    bigstruct.nontopsvec.clear();
    bigstruct.innerMinStatus = -1;
    bigstruct.outerMinStatus = -1;
    bigstruct.outerMin_Edm = -1;
}

inline void Fill_bigstruct(const input &in, const parameters_ttRC &pa,
                           big_struct &bigstruct)
{
    // Put in parameters and quantities for first top
    if (in.p.t1_lep) {
        bigstruct.tops.push_back( new top_system(
                    1, in.p.b1[0], in.p.b1[2], in.p.b1[1],
                    Make_mass(in.p.b1[0], in.p.b1[1], in.p.b1[2], in.p.b1[3]),
                    in.SD.b1[0], in.SD.b1[2], in.SD.b1[1],
                    in.p.d11[0], in.p.d11[2], in.p.d11[1],
                    Make_mass(in.p.d11[0], in.p.d11[1], in.p.d11[2], in.p.d11[3]),
                    in.SD.d11[0], in.SD.d11[2], in.SD.d11[1],
                    pa.m_t, pa.SD_m_t, pa.m_W, pa.SD_m_W) );
    } else {
        bigstruct.tops.push_back( new top_system(
                    0, in.p.b1[0], in.p.b1[2], in.p.b1[1],
                    Make_mass(in.p.b1[0], in.p.b1[1], in.p.b1[2], in.p.b1[3]),
                    in.SD.b1[0], in.SD.b1[2], in.SD.b1[1],
                    in.p.d11[0], in.p.d11[2], in.p.d11[1],
                    Make_mass(in.p.d11[0], in.p.d11[1], in.p.d11[2], in.p.d11[3]),
                    in.SD.d11[0], in.SD.d11[2], in.SD.d11[1],
                    pa.m_t, pa.SD_m_t, pa.m_W, pa.SD_m_W,
                    in.p.d12[0], in.p.d12[2], in.p.d12[1],
                    Make_mass(in.p.d12[0], in.p.d12[1], in.p.d12[2], in.p.d12[3]),
                    in.SD.d12[0], in.SD.d12[2], in.SD.d12[1]) );
    }

    // Put in parameters and quantities for second top
    if (in.p.t2_lep) {
        bigstruct.tops.push_back( new top_system(
                    1, in.p.b2[0], in.p.b2[2], in.p.b2[1],
                    Make_mass(in.p.b2[0], in.p.b2[1], in.p.b2[2], in.p.b2[3]),
                    in.SD.b2[0], in.SD.b2[2], in.SD.b2[1],
                    in.p.d21[0], in.p.d21[2], in.p.d21[1],
                    Make_mass(in.p.d21[0], in.p.d21[1], in.p.d21[2], in.p.d21[3]),
                    in.SD.d21[0], in.SD.d21[2], in.SD.d21[1],
                    pa.m_t, pa.SD_m_t, pa.m_W, pa.SD_m_W) );
    } else {
        bigstruct.tops.push_back( new top_system(
                    0, in.p.b2[0], in.p.b2[2], in.p.b2[1],
                    Make_mass(in.p.b2[0], in.p.b2[1], in.p.b2[2], in.p.b2[3]),
                    in.SD.b2[0], in.SD.b2[2], in.SD.b2[1],
                    in.p.d21[0], in.p.d21[2], in.p.d21[1],
                    Make_mass(in.p.d21[0], in.p.d21[1], in.p.d21[2], in.p.d21[3]),
                    in.SD.d21[0], in.SD.d21[2], in.SD.d21[1],
                    pa.m_t, pa.SD_m_t, pa.m_W, pa.SD_m_W,
                    in.p.d22[0], in.p.d22[2], in.p.d22[1],
                    Make_mass(in.p.d22[0], in.p.d22[1], in.p.d22[2], in.p.d22[3]),
                    in.SD.d22[0], in.SD.d22[2], in.SD.d22[1]) );
    }

    //Put in stuff for non-top objects
    vector<double> nontop_pts, nontop_etas, nontop_phis, nontop_ms;
    vector<double> nontop_ptwidths, nontop_etawidths, nontop_phiwidths;

    nontop_pts.push_back(in.p.bH1[0]);
    nontop_etas.push_back(in.p.bH1[2]);
    nontop_phis.push_back(in.p.bH1[1]);
    nontop_ms.push_back( Make_mass(in.p.bH1[0], in.p.bH1[1], in.p.bH1[2], in.p.bH1[3]) );
    nontop_ptwidths.push_back(in.SD.bH1[0]);
    nontop_etawidths.push_back(in.SD.bH1[2]);
    nontop_phiwidths.push_back(in.SD.bH1[1]);

    nontop_pts.push_back(in.p.bH2[0]);
    nontop_etas.push_back(in.p.bH2[2]);
    nontop_phis.push_back(in.p.bH2[1]);
    nontop_ms.push_back( Make_mass(in.p.bH2[0], in.p.bH2[1], in.p.bH2[2], in.p.bH2[3]) );
    nontop_ptwidths.push_back(in.SD.bH2[0]);
    nontop_etawidths.push_back(in.SD.bH2[2]);
    nontop_phiwidths.push_back(in.SD.bH2[1]);

    for (auto p : in.p.p_others) {
        nontop_pts.push_back(p[0]);
        nontop_etas.push_back(p[2]);
        nontop_phis.push_back(p[1]);
        nontop_ms.push_back( Make_mass(p[0], p[1], p[2], p[3]) );
        nontop_ptwidths.push_back(p[0]);
        nontop_etawidths.push_back(p[2]);
        nontop_phiwidths.push_back(p[1]);
    }

    //Put in MET
    bigstruct.MET_px = in.MET_px;
    bigstruct.MET_py = in.MET_py;
    
    bigstruct.nontops_ptr = new nontop_system(
                                      nontop_pts, nontop_etas, nontop_phis, nontop_ms,
                                      nontop_ptwidths, nontop_etawidths, nontop_phiwidths);

}

inline void LV_to_cyl(const LV &in, double out[4])
{
    out[0] = in.Pt();
    out[1] = in.Phi();
    out[2] = in.Eta();
    out[3] = in.E();
}

inline void cyl_to_LV(const double (&in)[4], LV &out)
{
    out.SetPtEtaPhiE( in[0], in[2], in[1], in[3] );
}


void ttRC::daughter_to_parents(const input_x<4> &daughters, parents<4> &prnts)
{
    LV t1, t2, w1, w2, b1, b2, d11, d12, d21, d22, bH1, bH2, h;

    cyl_to_LV(daughters.b1, b1);
    cyl_to_LV(daughters.d11, d11);
    cyl_to_LV(daughters.d12, d12);
    cyl_to_LV(daughters.b2, b2);
    cyl_to_LV(daughters.d21, d21);
    cyl_to_LV(daughters.d22, d22);
    cyl_to_LV(daughters.bH1, bH1);
    cyl_to_LV(daughters.bH2, bH2);

    w1 = d11 + d12;
    t1 = b1 + d11 + d12;
    w2 = d21 + d22;
    t2 = b2 + d21 + d22;
    h = bH1 + bH2;

    cout << "daughter d11 input Pt " << daughters.d11[0] << endl;
    //cout << "daughter d21 input Phi " << daughters.d21[1] << endl;
    //cout << "daughter d21 input Eta " << daughters.d21[2] << endl;
    //cout << "daughter d21 input E " << daughters.d21[3] << endl;
    cout << "daughter d11 Pt " << d11.Pt() << endl;
    cout <<"daughter w1 Pt " << w1.Pt() << endl;
    cout << "daughter d11 Px " << d11.Px() << endl;
    cout << "daughter d12 Px " << d12.Px() << endl;
    cout <<"daughter w1 Px " << w1.Px() << endl;

    LV_to_cyl(w1, prnts.w1);
    LV_to_cyl(t1, prnts.t1);
    LV_to_cyl(w2, prnts.w2);
    LV_to_cyl(t2, prnts.t2);
    LV_to_cyl(h, prnts.h);
}

void ttRC::input_to_output(const input &in, output &out)
{
    out.p = in.p;
    daughter_to_parents(in.p, out.parents_p);
}

void ttRC::met_to_neutrino(const double &met_px, const double &met_py, double neu[4])
{
    LV neutrino;
    double pt = Make_pT( met_px, met_py );
    double phi = Make_phi( met_px, met_py );
    neutrino.SetPtEtaPhiM( pt, 0, phi, 0 );
    LV_to_cyl(neutrino, neu);
}

inline void Plot_chi2s(topEventMinimizer *ev)
{
    ev->plot_chi2s();
}

inline void Get_results(const input &in, big_struct &bigstruct, topEventMinimizer *ev, output &out)
{
    out.p.t1_lep = in.p.t1_lep;
    out.p.t2_lep = in.p.t2_lep;

    out.weight = 1.0/ev->get_best_total_chi2();

    // top 1
    {
        // const auto t1 = ev->getConverter("getTop", 0);
        const auto b = ev->get_b(0);
        // const auto W = ev->getConverter("getW", 0);
        const auto d1 = ev->get_Wd1(0);
        const auto d2 = ev->get_Wd2(0);
        LV_to_cyl(b, out.p.b1);
        LV_to_cyl(d1, out.p.d11);
        LV_to_cyl(d2, out.p.d12);

        const auto t = ev->get_top(0);
        const auto w = ev->get_W(0);
        LV_to_cyl(t, out.parents_p.t1);
        LV_to_cyl(w, out.parents_p.w1);

        out.chi2s.b1 = ev->get_best_b_chi2(0);
        out.chi2s.d11 = ev->get_best_Wd1_chi2(0);
        out.chi2s.d12 = ev->get_best_Wd2_chi2(0);
        out.chi2s.mt1 = ev->get_best_mTop_chi2(0);
        out.chi2s.mw1 = ev->get_best_mW_chi2(0);

    }
    // top 2
    {
        // const auto t1 = ev->getConverter("getTop", 1);
        const auto b = ev->get_b(1);
        // const auto W = ev->getConverter("getW", 1);
        const auto d1 = ev->get_Wd1(1);
        const auto d2 = ev->get_Wd2(1);
        LV_to_cyl(b, out.p.b2);
        LV_to_cyl(d1, out.p.d21);
        LV_to_cyl(d2, out.p.d22);

        const auto t = ev->get_top(1);
        const auto w = ev->get_W(1);
        LV_to_cyl(t, out.parents_p.t2);
        LV_to_cyl(w, out.parents_p.w2);

        out.chi2s.b2 = ev->get_best_b_chi2(1);
        out.chi2s.d21 = ev->get_best_Wd1_chi2(1);
        out.chi2s.d22 = ev->get_best_Wd2_chi2(1);
        out.chi2s.mt2 = ev->get_best_mTop_chi2(1);
        out.chi2s.mw2 = ev->get_best_mW_chi2(1);

    }

    const auto bH1 = ev->get_nontop_object(0);
    const auto bH2 = ev->get_nontop_object(1);
    LV_to_cyl(bH1, out.p.bH1);
    LV_to_cyl(bH2, out.p.bH2);

    const auto h = ev->get_higgs();
    LV_to_cyl(h, out.parents_p.h);

    out.chi2s.nontop_objs = ev->get_best_nontop_chi2();

    out.inner_min_status = ev->get_inner_min_status();
    out.outer_min_status = ev->get_outer_min_status();
    out.inner_edm = ev->get_inner_edm();
    out.outer_edm = ev->get_outer_edm();

}


output ttRC::reco(const input &in, const parameters_ttRC &pa)
{
    big_struct bs;
    Reset_bigstruct(bs);
    Fill_bigstruct(in, pa, bs);

    int debug_verbosity = 2;
    rc = new topEventMinimizer(bs, debug_verbosity);
    rc->create_tops();
    rc->minimizeTotalChiSquare();

    output res;
    Get_results(in, bs, rc, res);

    //Plot_chi2s(rc);

    delete rc;

    return res;
}

}

#endif
