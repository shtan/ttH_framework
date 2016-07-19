inline void Boost_(const double *boost, double *vector)
{
    const double Boost_trigger_gamma = 1e-10; // min. value of 1-gamma to
                                              // trigger boost

    double ovec[4];
    for (int count = 0; count < 4; ++count)
        ovec[count] = vector[count];

    const double beta2 =
        boost[0] * boost[0] + boost[1] * boost[1] + boost[2] * boost[2];
    const double gamma = 1 / sqrt(1 - beta2);

    if ((gamma - 1) > Boost_trigger_gamma) {
        vector[0] = gamma * (ovec[0] - boost[0] * ovec[1] - boost[1] * ovec[2] -
                             boost[2] * ovec[3]);
        vector[1] = -gamma * boost[0] * ovec[0] +
                    (1 + (gamma - 1) * boost[0] * boost[0] / beta2) * ovec[1] +
                    (gamma - 1) * boost[0] / beta2 *
                        (boost[1] * ovec[2] + boost[2] * ovec[3]);
        vector[2] = -gamma * boost[1] * ovec[0] +
                    (1 + (gamma - 1) * boost[1] * boost[1] / beta2) * ovec[2] +
                    (gamma - 1) * boost[1] / beta2 *
                        (boost[0] * ovec[1] + boost[2] * ovec[3]);
        vector[3] = -gamma * boost[2] * ovec[0] +
                    (1 + (gamma - 1) * boost[2] * boost[2] / beta2) * ovec[3] +
                    (gamma - 1) * boost[2] / beta2 *
                        (boost[0] * ovec[1] + boost[1] * ovec[2]);
    }
}

inline void unboost_sys_v0(handleEvent &evh, double *boost)
{
    const unsigned int bp_sz = 8;
    const string boost_partc[bp_sz] = {"bottom",     "antiBottom", "lepton",
                                       "antiLepton", "qFromW",     "qbarFromW",
                                       "bFromH",     "bbarFromH"};
    double p[bp_sz][4];
    double tot_p[4] = {0, 0, 0, 0};
    for (unsigned int i = 0; i < bp_sz; ++i) {
        p[i][0] = evh.smearedParticles[boost_partc[i]]->E();
        p[i][1] = evh.smearedParticles[boost_partc[i]]->Px();
        p[i][2] = evh.smearedParticles[boost_partc[i]]->Py();
        p[i][3] = evh.smearedParticles[boost_partc[i]]->Pz();
        tot_p[0] += p[i][0];
        //         tot_p[1] += p[i][1];
        //         tot_p[2] += p[i][2];
        tot_p[3] += p[i][3];
    }

    boost[0] = 0;
    boost[1] = 0;
    boost[2] = tot_p[3] / tot_p[0];
    for (unsigned int i = 0; i < bp_sz; ++i) {
        Boost_(boost, p[i]);
        evh.smearedParticles[boost_partc[i]]->SetE(p[i][0]);
        evh.smearedParticles[boost_partc[i]]->SetPx(p[i][1]);
        evh.smearedParticles[boost_partc[i]]->SetPy(p[i][2]);
        evh.smearedParticles[boost_partc[i]]->SetPz(p[i][3]);
    }
}

inline void boost_v0(handleEvent &evh, double *boost)
{
    const unsigned int bp_sz = 8;
    const string boost_partc[bp_sz] = {"bottom",     "antiBottom", "lepton",
                                       "antiLepton", "qFromW",     "qbarFromW",
                                       "bFromH",     "bbarFromH"};
    double p[bp_sz][4];
    double tot[4] = {0, 0, 0, 0};
    for (unsigned int i = 0; i < bp_sz; ++i) {
        p[i][0] = evh.smearedParticles[boost_partc[i]]->E();
        p[i][1] = evh.smearedParticles[boost_partc[i]]->Px();
        p[i][2] = evh.smearedParticles[boost_partc[i]]->Py();
        p[i][3] = evh.smearedParticles[boost_partc[i]]->Pz();
        tot[0] += p[i][0];
        tot[1] += p[i][1];
        tot[2] += p[i][2];
        tot[3] += p[i][3];
        cout << tot[0] << " " << tot[1] << " " << tot[2] << endl;
    }

    tot[0] = 0;
    tot[1] = 0;
    tot[2] = 0;
    tot[3] = 0;
    for (unsigned int i = 0; i < bp_sz; ++i) {
        Boost_(boost, p[i]);
        evh.smearedParticles[boost_partc[i]]->SetE(p[i][0]);
        evh.smearedParticles[boost_partc[i]]->SetPx(p[i][1]);
        evh.smearedParticles[boost_partc[i]]->SetPy(p[i][2]);
        evh.smearedParticles[boost_partc[i]]->SetPz(p[i][3]);
        tot[0] += p[i][0];
        tot[1] += p[i][1];
        tot[2] += p[i][2];
        tot[3] += p[i][3];
        cout << tot[0] << " " << tot[1] << " " << tot[2] << endl;
    }
}

inline void boost_v0_all_best(handleEvent &evh, double *boost)
{
    double p[4];
    for (auto m : evh.bestParticles) {
        auto part = m.second;
        p[0] = part->E();
        p[1] = part->Px();
        p[2] = part->Py();
        p[3] = part->Pz();

        Boost_(boost, p);
        part->SetE(p[0]);
        part->SetPx(p[1]);
        part->SetPy(p[2]);
        part->SetPz(p[3]);
    }
}