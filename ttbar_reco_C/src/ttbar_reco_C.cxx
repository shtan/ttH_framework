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

    cout << "MET pt = " << pow( ( pow(initial.p.MET_px,2) + pow(initial.p.MET_py,2) ), 0.5 ) << endl;
    cout << "MET phi = " << atan2( initial.p.MET_py, initial.p.MET_px ) << endl;

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

    cout << "MET pt = " << pow( ( pow(initial.p.MET_px,2) + pow(initial.p.MET_py,2) ), 0.5 ) << endl;
    cout << "MET phi = " << atan2( initial.p.MET_py, initial.p.MET_px ) << endl;

}

// void Print_diff(const input &in, const output &res)
// {
//     
// }

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

/*inline void Copy_data_nontop(const input &in, vector<LV> &non_tt,
                             vector<double> &SDs_pT, vector<double> &SDs_phi,
                             vector<double> &SDs_eta)
{
    {
        auto it_LV = non_tt.begin();
        Insert_cartesian_XYZT(in.p.bH1, it_LV);
        ++it_LV;

        Insert_cartesian_XYZT(in.p.bH2, it_LV);
        ++it_LV;

        for (auto p : in.p.p_others) {
            Insert_cartesian_XYZT(p, it_LV);
            ++it_LV;
        }
    }

    {
        auto it_pT = SDs_pT.begin();
        auto it_phi = SDs_phi.begin();
        auto it_eta = SDs_eta.begin();

        *it_pT = in.SD.bH1[0];
        *it_phi = in.SD.bH1[1];
        *it_eta = in.SD.bH1[2];
        ++it_pT;
        ++it_phi;
        ++it_eta;

        *it_pT = in.SD.bH2[0];
        *it_phi = in.SD.bH2[1];
        *it_eta = in.SD.bH2[2];
        ++it_pT;
        ++it_phi;
        ++it_eta;
        for (auto SD : in.SD.p_others) {
            *it_pT = SD[0];
            *it_phi = SD[1];
            *it_eta = SD[2];
            ++it_pT;
            ++it_phi;
            ++it_eta;
        }
    }
}*/

// [4]: pT, phi, eta, E
/*inline void Load_event_data_top(const bool lepton, const double b[4],
                                const double d1[4], const double d2[4],
                                const double SD_b[4], const double SD_d1[4],
                                const double SD_d2[4],
                                const parameters_ttRC &pa,
                                topEventMinimizer *ev)
{
    const auto &m_t = pa.m_t;
    const auto &m_W = pa.m_W;
    const auto &SD_m_t = pa.SD_m_t;
    const auto &SD_m_W = pa.SD_m_W;

    const auto b1_px = Make_px(b[0], b[1]);
    const auto b1_py = Make_py(b[0], b[1]);
    const auto b1_pz = Make_pz(b[0], b[2]);
    const auto &b1_E = b[3];
    const auto &b1_SD_pT = SD_b[0];
    const auto &b1_SD_phi = SD_b[1];
    const auto &b1_SD_eta = SD_b[2];

    const auto d1_px = Make_px(d1[0], d1[1]);
    const auto d1_py = Make_py(d1[0], d1[1]);
    const auto d1_pz = Make_pz(d1[0], d1[2]);
    const auto &d1_E = d1[3];
    const auto &d1_SD_pT = SD_d1[0];
    const auto &d1_SD_phi = SD_d1[1];
    const auto &d1_SD_eta = SD_d1[2];
    if (lepton) {
        ev->addLeptonicTop(b1_px, b1_py, b1_pz, b1_E, b1_SD_pT, b1_SD_eta,
                           b1_SD_phi, d1_px, d1_py, d1_pz, d1_E, d1_SD_pT,
                           d1_SD_eta, d1_SD_phi, m_t, SD_m_t, m_W, SD_m_W);
    } else {
        const auto d2_px = Make_px(d2[0], d2[1]);
        const auto d2_py = Make_py(d2[0], d2[1]);
        const auto d2_pz = Make_pz(d2[0], d2[2]);
        const auto &d2_E = d2[3];
        const auto &d2_SD_pT = SD_d2[0];
        const auto &d2_SD_phi = SD_d2[1];
        const auto &d2_SD_eta = SD_d2[2];

        ev->addHadronicTop(b1_px, b1_py, b1_pz, b1_E, b1_SD_pT, b1_SD_eta,
                           b1_SD_phi, d1_px, d1_py, d1_pz, d1_E, d1_SD_pT,
                           d1_SD_eta, d1_SD_phi, d2_px, d2_py, d2_pz, d2_E,
                           d2_SD_pT, d2_SD_eta, d2_SD_phi, m_t, SD_m_t, m_W,
                           SD_m_W);
    }
}*/

/*inline void Load_event_data_tt(const input &in, const parameters_ttRC &pa,
                               topEventMinimizer *ev)
{
    Load_event_data_top(in.p.t1_lep, in.p.b1, in.p.d11, in.p.d12, in.SD.b1,
                        in.SD.d11, in.SD.d12, pa, ev);
    Load_event_data_top(in.p.t2_lep, in.p.b2, in.p.d21, in.p.d22, in.SD.b2,
                        in.SD.d21, in.SD.d22, pa, ev);
}*/

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

/*    //Put in neutrinos as MET
    double d12_px = Make_px( in.p.d12[0], in.p.d12[1] );
    double d12_py = Make_py( in.p.d12[0], in.p.d12[1] );
    double d22_px = Make_px( in.p.d22[0], in.p.d22[1] );
    double d22_py = Make_py( in.p.d22[0], in.p.d22[1] );
    bigstruct.MET_px = 0;
    bigstruct.MET_py = 0;
    if (in.p.t1_lep) {
        bigstruct.MET_px += d12_px;
        bigstruct.MET_py += d12_py;
    }
    if (in.p.t2_lep) {
        bigstruct.MET_px += d22_px;
        bigstruct.MET_py += d22_py;
    }*/

    //Put in MET
    bigstruct.MET_px = in.p.MET_px;
    bigstruct.MET_py = in.p.MET_py;
    
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
    }

    const auto bH1 = ev->get_nontop_object(0);
    const auto bH2 = ev->get_nontop_object(1);
    LV_to_cyl(bH1, out.p.bH1);
    LV_to_cyl(bH2, out.p.bH2);

}

/*inline void Get_results(const input &in, topEventMinimizer *ev, output &out)
{
    out.p.t1_lep = in.p.t1_lep;
    out.p.t2_lep = in.p.t2_lep;
    // Minimizer stuff
    //     ev->getNonTopChiSquare();
    //     ev->getHadronicChiSquare();
    //     ev->getTopMassChiSquare();
    //     ev->getTopChiSquare();
    //
    //     ev->getOneTopMassChiSquare(0);
    //     ev->getOneBChiSquare(0);
    //     ev->getOneWMassChiSquare(0);
    //     ev->getOneWDaughter1ChiSquare(0);
    //
    //     ev->getOneTopMassChiSquare(1);
    //     ev->getOneBChiSquare(1);
    //     ev->getOneWMassChiSquare(1);
    //     ev->getOneWDaughter1ChiSquare(1);

    out.weight = 1.0 / ev->getChiSquare(); // total

    // top 1
    {
        // const auto t1 = ev->getConverter("getTop", 0);
        const auto b = ev->getConverter("getBJet", 0);
        // const auto W = ev->getConverter("getW", 0);
        const auto d1 = ev->getConverter("getWDaughter1", 0);
        const auto d2 = ev->getConverter("getWDaughter2", 0);
        LV_to_cyl(b, out.p.b1);
        LV_to_cyl(d1, out.p.d11);
        LV_to_cyl(d2, out.p.d12);
    }
    // top 2
    {
        // const auto t1 = ev->getConverter("getTop", 1);
        const auto b = ev->getConverter("getBJet", 1);
        // const auto W = ev->getConverter("getW", 1);
        const auto d1 = ev->getConverter("getWDaughter1", 1);
        const auto d2 = ev->getConverter("getWDaughter2", 1);
        LV_to_cyl(b, out.p.b2);
        LV_to_cyl(d1, out.p.d21);
        LV_to_cyl(d2, out.p.d22);
    }

    const auto bH1 = ev->getConverter("getNonTopObject4", 0);
    const auto bH2 = ev->getConverter("getNonTopObject4", 1);
    LV_to_cyl(bH1, out.p.bH1);
    LV_to_cyl(bH2, out.p.bH2);

    //     const auto sz = in.p.p_others.size();
    //     if (sz > 0) {
    //         out.p.p_others.reserve(sz);
    //         for (unsigned int i = 0; i < sz; ++i) {
    //             const auto lv = ev->getConverter("getNonTopObject4", 2 + i);
    //             double *p = new double[4];
    //             LV_to_cyl(lv, p);
    //             out.p.p_others.push_back(p);
    //         }
    //     }
}*/

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
    delete rc;

    return res;
}

/*output ttRC::reco(const input &in, const parameters_ttRC &pa)
{
    const auto non_tt_sz = in.p.p_others.size() + 2; // + 2 bH
    vector<LV> non_tt(non_tt_sz);
    vector<double> non_tt_SD_pT(non_tt_sz);
    vector<double> non_tt_SD_phi(non_tt_sz);
    vector<double> non_tt_SD_eta(non_tt_sz);

    Copy_data_nontop(in, non_tt, non_tt_SD_pT, non_tt_SD_phi, non_tt_SD_eta);

    double p_sum_tt[3];
    Eval_p_sum_tt_cartesian(in.p, p_sum_tt);

    rc = new topEventMinimizer(
        non_tt, non_tt_SD_pT, non_tt_SD_eta, non_tt_SD_phi, pa.m_t, pa.SD_m_t,
        pa.m_W, pa.SD_m_W, p_sum_tt[0], p_sum_tt[1], p_sum_tt[2]);

    Load_event_data_tt(in, pa, rc);

    rc->initializeDeltas();
    rc->minimizeTotalChiSquare();

    output res;
    Get_results(in, rc, res);
    delete rc;

    return res;
}*/

// void topReconstructionFromLHE::Loop_fill_results_SM(topEventMinimizer &ev,
//                                                     handleEvent &evh)
// {
//     // Minimizer stuff
//     evh.chiSquareds["nonTop"] = ev.getNonTopChiSquare();
//     evh.chiSquareds["hadronic"] = ev.getHadronicChiSquare();
//     evh.chiSquareds["topMass"] = ev.getTopMassChiSquare();
//     evh.chiSquareds["topSystem"] = ev.getTopChiSquare();
//     evh.chiSquareds["total"] = ev.getChiSquare();
//     evh.chiSquareds["qbarFromW"] = ev.getHadronicChiSquare();
//
//     evh.chiSquareds["top_topMass"] = ev.getOneTopMassChiSquare(0);
//     evh.chiSquareds["antiTop_topMass"] = ev.getOneTopMassChiSquare(1);
//     evh.chiSquareds["bottom"] = ev.getOneBChiSquare(0);
//     evh.chiSquareds["antiBottom"] = ev.getOneBChiSquare(1);
//     evh.chiSquareds["Wplus_Wmass"] = ev.getOneWMassChiSquare(0);
//     evh.chiSquareds["Wminus_Wmass"] = ev.getOneWMassChiSquare(1);
//
//     if (evh.leptonFlag == 0) {
//         evh.chiSquareds["leptonicTopMass"] = ev.getOneTopMassChiSquare(0);
//         evh.chiSquareds["hadronicTopMass"] = ev.getOneTopMassChiSquare(1);
//         evh.chiSquareds["leptonicBottom"] = ev.getOneBChiSquare(0);
//         evh.chiSquareds["hadronicBottom"] = ev.getOneBChiSquare(1);
//         evh.chiSquareds["leptonicWMass"] = ev.getOneWMassChiSquare(0);
//         evh.chiSquareds["hadronicWMass"] = ev.getOneWMassChiSquare(1);
//         evh.chiSquareds["qFromW"] = ev.getOneWDaughter1ChiSquare(1);
//         evh.chiSquareds["lepton"] = ev.getOneWDaughter1ChiSquare(0);
//         evh.chiSquareds["lepton_or_antilepton"] =
//             ev.getOneWDaughter1ChiSquare(0);
//
//     } else if (evh.leptonFlag == 1) {
//         evh.chiSquareds["leptonicTopMass"] = ev.getOneTopMassChiSquare(1);
//         evh.chiSquareds["hadronicTopMass"] = ev.getOneTopMassChiSquare(0);
//         evh.chiSquareds["leptonicBottom"] = ev.getOneBChiSquare(1);
//         evh.chiSquareds["hadronicBottom"] = ev.getOneBChiSquare(0);
//         evh.chiSquareds["leptonicWMass"] = ev.getOneWMassChiSquare(1);
//         evh.chiSquareds["hadronicWMass"] = ev.getOneWMassChiSquare(0);
//         evh.chiSquareds["qFromW"] = ev.getOneWDaughter1ChiSquare(0);
//         evh.chiSquareds["antiLepton"] = ev.getOneWDaughter1ChiSquare(1);
//         evh.chiSquareds["lepton_or_antilepton"] =
//             ev.getOneWDaughter1ChiSquare(1);
//     }
//
//     // for (vector<string>::const_iterator t = chinames.begin(); t <
//     // chinames.end(); t++){
//     //    cout<<*t<<" " << *(evh.chiSquareds[*t])<<endl;
//     //}
//
//     *evh.bestParticles["top"] = ev.getConverter("getTop", 0);
//     *evh.bestParticles["bottom"] = ev.getConverter("getBJet", 0);
//     *evh.bestParticles["Wplus"] = ev.getConverter("getW", 0);
//     // FIXME
//     if (evh.leptonFlag == false) {
//         *evh.bestParticles["antiLepton"] = ev.getConverter("getWDaughter1",
//         0);
//         *evh.bestParticles["neutrino"] = ev.getConverter("getWDaughter2", 0);
//     } else {
//         *evh.bestParticles["qFromW"] = ev.getConverter("getWDaughter1", 0);
//         *evh.bestParticles["qbarFromW"] = ev.getConverter("getWDaughter2",
//         0);
//     }
//     *evh.bestParticles["bFromH"] = ev.getConverter("getNonTopObject4", 0);
//
//     *evh.bestParticles["antiTop"] = ev.getConverter("getTop", 1);
//     *evh.bestParticles["antiBottom"] = ev.getConverter("getBJet", 1);
//     *evh.bestParticles["Wminus"] = ev.getConverter("getW", 1);
//     // FIXME
//     if (evh.leptonFlag == false) {
//         *evh.bestParticles["qFromW"] = ev.getConverter("getWDaughter1", 1);
//         *evh.bestParticles["qbarFromW"] = ev.getConverter("getWDaughter2",
//         1);
//     } else {
//         *evh.bestParticles["lepton"] = ev.getConverter("getWDaughter1", 1);
//         *evh.bestParticles["antiNeutrino"] =
//             ev.getConverter("getWDaughter2", 1);
//     }
//     *evh.bestParticles["bbarFromH"] = ev.getConverter("getNonTopObject4", 1);
// }
}

#endif
