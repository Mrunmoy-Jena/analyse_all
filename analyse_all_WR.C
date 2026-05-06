#include "analyse_all_WR.h"
#include "TSystem.h"
#include<algorithm>
using namespace std;

void run_analysis(TString inputFile, TFile* outputFile)
{
    TChain * tree = new TChain("evt");
    
    tree->Add(inputFile);
    
    std::cout << "\n\n-- Analyse TTree with " << tree->GetEntries() << " entries \n";

    // Create histograms
    //auto array_th2_tofd = CreateHistograms<TH2F>("bar", 45, 1000,-145,-137,1000,31,35);
    //auto array_th1_tofd = CreateHistograms<TH1F>("toffset_bar", 45, 1000,170,178);
    //auto h2_tofd_pid = CreateHistograms<TH2F>("tofd_pid", 1,  45,0,45,1000,3,12);
    auto h2_charge_foot_vs_tofd = CreateHistograms<TH2F>("h2_charge_foot_vs_tofd", 1, 1000,0,15,1000,0,15);
    auto h2_frs_pid = CreateHistograms<TH2F>("frs_pid", 1, 1000, 2.4, 3.2, 1000, 3, 12);
    auto h1_neuland = CreateHistograms<TH1F>("h1_NeuLAND_ana", 3, 2000,-1000,1000);
    
   /* auto h1_dwr_main_foot = CreateHistograms<TH1F>("main - foot", 16 ,1000,-500,500);
    auto h1_dwr_main_califa = CreateHistograms<TH1F>("main - califa", 4, 1000,-2000,2000);
    auto h1_dwr_main_califa_2hits = CreateHistograms<TH1F>("main - califa (2 hits in CALIFA)", 4, 1000,-2000,2000);
    auto h1_dwr_main_califa_25F_p2p_24O = CreateHistograms<TH1F>("main - califa 25F(p2p)24O", 4, 1000,-2000,2000);
    auto h1_dwr_main_alpide = CreateHistograms<TH1F>("main - alpide", 2, 1000,-500,500);
    auto h1_dwr_main_los = CreateHistograms<TH1F>("main - los", 1, 1000,-500,500);
    auto h1_dwr_main_s2 = CreateHistograms<TH1F>("main - s2", 1, 1000,-500,500);
    */

        auto h1_califa_gamma_corr_25F_p2p_24O_mult_exact = CreateHistograms<TH1F>("h1_califa_gamma_corr_25F_p2p_24O_mult_exact", 10, 200, 0, 10);
        auto h1_califa_gamma_corr_25F_p2p_24O_mult_leq = CreateHistograms<TH1F>("h1_califa_gamma_corr_25F_p2p_24O_mult_leq", 10, 200, 0, 10);
    auto h2_neuland = CreateHistograms<TH2F>("h2_NeuLAND_ana", 3, 1000,-1000,1000,1000,-1000,1000);
    auto h1_neuland_Z = CreateHistograms<TH1F>("h1_NeuLAND_ana_Z", 3, 1000,1500,1800);
    auto h2_neuland_Z = CreateHistograms<TH2F>("h2_NeuLAND_ana_Z", 3, 1000,1500,1800,1000,-500,500);
    
    auto h2_AoZ_vs_Z = CreateHistograms<TH2F>("h2_AoZ_vs_Z", 1, 1000, 1.5,3.8,1000,1,12);
    
    auto h2_califa_theta_corr = CreateHistograms<TH2F>("h2_califa_theta_corr",1, 100,0,100,100,0,100);
    auto h1_califa_opang = CreateHistograms<TH1F>("h1_califa_opang", 1 , 200,0,200);
    auto h1_califa_deltaphi = CreateHistograms<TH1F>("h1_califa_deltaphi", 1 , 360,0,360);
    auto h1_califa_deltaphifolded = CreateHistograms<TH1F>("h1_califa_deltaphifolded", 1 ,250,0,250);
    auto h1_califa_opang_25F_p2p_24O = CreateHistograms<TH1F>("h1_califa_opang_25F_p2p_24O", 1 , 200,0,200);
    auto h1_califa_opang_25F_p2p_23O = CreateHistograms<TH1F>("h1_califa_opang_25F_p2p_23O", 1 , 200,0,200);
    auto h1_califa_opang_25F_ppn_24F = CreateHistograms<TH1F>("h1_califa_opang_25F_ppn_24F", 1 , 200,0,200);
    auto h1_califa_opang_25F_p2p_22O = CreateHistograms<TH1F>("h1_califa_opang_25F_p2p_22O", 1 , 200,0,200);
    auto h2_califa_theta_corr_25F_p2p_24O = CreateHistograms<TH2F>("h2_califa_theta_corr_25F_p2p_24O",1, 100,0,100,100,0,100);
    auto h2_califa_theta_corr_25F_p2p_23O = CreateHistograms<TH2F>("h2_califa_theta_corr_25F_p2p_23O",1, 100,0,100,100,0,100);
    auto h2_califa_theta_corr_25F_ppn_24F = CreateHistograms<TH2F>("h2_califa_theta_corr_25F_ppn_24F",1, 100,0,100,100,0,100);
    auto h2_califa_E_vs_theta = CreateHistograms<TH2F>("h2_califa_E_vs_theta",1,100,0,100,1000,0,1000);
    auto h2_califa_Ecorr = CreateHistograms<TH2F>("h2_califa_Ecorr",1, 1000,0,1000,1000,0,1000);
    auto h2_califa_Ecorr_25F_p2p_24O = CreateHistograms<TH2F>("h2_califa_Ecorr_25F_p2p_24O",1, 1000,0,1000,1000,0,1000);
    auto h2_califa_Ecorr_25F_p2p_23O = CreateHistograms<TH2F>("h2_califa_Ecorr_25F_p2p_23O",1, 1000,0,1000,1000,0,1000);
    auto h2_califa_Ecorr_25F_ppn_24F = CreateHistograms<TH2F>("h2_califa_Ecorr_25F_ppn_24F",1, 1000,0,1000,1000,0,1000);

    // Set axis titles for coincidence plots
    h2_califa_theta_corr[0]->GetXaxis()->SetTitle("Theta1 [deg]");
    h2_califa_theta_corr[0]->GetYaxis()->SetTitle("Theta2 [deg]");
    h2_califa_theta_corr_25F_p2p_24O[0]->GetXaxis()->SetTitle("Theta1 [deg]");
    h2_califa_theta_corr_25F_p2p_24O[0]->GetYaxis()->SetTitle("Theta2 [deg]");
    h2_califa_theta_corr_25F_p2p_23O[0]->GetXaxis()->SetTitle("Theta1 [deg]");
    h2_califa_theta_corr_25F_p2p_23O[0]->GetYaxis()->SetTitle("Theta2 [deg]");
    h2_califa_theta_corr_25F_ppn_24F[0]->GetXaxis()->SetTitle("Theta1 [deg]");
    h2_califa_theta_corr_25F_ppn_24F[0]->GetYaxis()->SetTitle("Theta2 [deg]");
    h2_califa_E_vs_theta[0]->GetXaxis()->SetTitle("Theta [deg]");
    h2_califa_E_vs_theta[0]->GetYaxis()->SetTitle("Energy [MeV]");
    h2_califa_Ecorr[0]->GetXaxis()->SetTitle("Energy1 [MeV]");
    h2_califa_Ecorr[0]->GetYaxis()->SetTitle("Energy2 [MeV]");
    h2_califa_Ecorr_25F_p2p_24O[0]->GetXaxis()->SetTitle("Energy1 [MeV]");
    h2_califa_Ecorr_25F_p2p_24O[0]->GetYaxis()->SetTitle("Energy2 [MeV]");
    h2_califa_Ecorr_25F_p2p_23O[0]->GetXaxis()->SetTitle("Energy1 [MeV]");
    h2_califa_Ecorr_25F_p2p_23O[0]->GetYaxis()->SetTitle("Energy2 [MeV]");
    h2_califa_Ecorr_25F_ppn_24F[0]->GetXaxis()->SetTitle("Energy1 [MeV]");
    h2_califa_Ecorr_25F_ppn_24F[0]->GetYaxis()->SetTitle("Energy2 [MeV]");
    h1_califa_opang[0]->GetXaxis()->SetTitle("Opening Angle [deg]");
    h1_califa_opang[0]->GetYaxis()->SetTitle("Counts");
    h1_califa_opang_25F_p2p_24O[0]->GetXaxis()->SetTitle("Opening Angle [deg]");
    h1_califa_opang_25F_p2p_24O[0]->GetYaxis()->SetTitle("Counts");
    h1_califa_opang_25F_p2p_23O[0]->GetXaxis()->SetTitle("Opening Angle [deg]");
    h1_califa_opang_25F_p2p_23O[0]->GetYaxis()->SetTitle("Counts");
    h1_califa_opang_25F_ppn_24F[0]->GetXaxis()->SetTitle("Opening Angle [deg]");
    h1_califa_opang_25F_ppn_24F[0]->GetYaxis()->SetTitle("Counts");
    h1_califa_opang_25F_p2p_22O[0]->GetXaxis()->SetTitle("Opening Angle [deg]");
    h1_califa_opang_25F_p2p_22O[0]->GetYaxis()->SetTitle("Counts");
    h1_califa_deltaphi[0]->GetXaxis()->SetTitle("Delta Phi [deg]");
    h1_califa_deltaphi[0]->GetYaxis()->SetTitle("Counts");
    h1_califa_deltaphifolded[0]->GetXaxis()->SetTitle("Folded Delta Phi [deg]");
    h1_califa_deltaphifolded[0]->GetYaxis()->SetTitle("Counts");
   
    auto h2_califa_gamma_all = CreateHistograms<TH2F>("h2_califa_gamma_all",1,2000,0,100,200,0,10);
    auto h2_califa_gamma_25F_p2p_24O = CreateHistograms<TH2F>("h2_califa_gamma_25F_p2p_24O",1,2000,0,100,200,0,10);
    auto h2_califa_gamma_25F_p2p_23O = CreateHistograms<TH2F>("h2_califa_gamma_25F_p2p_23O",1,2000,0,100,200,0,10);
    auto h2_califa_gamma_25F_ppn_24F = CreateHistograms<TH2F>("h2_califa_gamma_25F_ppn_24F",1,2000,0,100,200,0,10);
    auto h1_califa_gamma_all = CreateHistograms<TH1F>("h1_califa_gamma_all",1,200,0,10);
    auto h1_califa_gamma_25F_p2p_24O = CreateHistograms<TH1F>("h1_califa_gamma_25F_p2p_24O",1,200,0,10);
    auto h1_califa_gamma_25F_p2p_23O = CreateHistograms<TH1F>("h1_califa_gamma_25F_p2p_23O",1,200,0,10);
    auto h1_califa_gamma_25F_ppn_24F = CreateHistograms<TH1F>("h1_califa_gamma_25F_ppn_24F",1,200,0,10);

    // Set axis titles for gamma histograms
    h2_califa_gamma_all[0]->GetXaxis()->SetTitle("Theta [deg]");
    h2_califa_gamma_all[0]->GetYaxis()->SetTitle("Energy [MeV]");
    h2_califa_gamma_25F_p2p_24O[0]->GetXaxis()->SetTitle("Theta [deg]");
    h2_califa_gamma_25F_p2p_24O[0]->GetYaxis()->SetTitle("Energy [MeV]");
    h2_califa_gamma_25F_p2p_23O[0]->GetXaxis()->SetTitle("Theta [deg]");
    h2_califa_gamma_25F_p2p_23O[0]->GetYaxis()->SetTitle("Energy [MeV]");
    h2_califa_gamma_25F_ppn_24F[0]->GetXaxis()->SetTitle("Theta [deg]");
    h2_califa_gamma_25F_ppn_24F[0]->GetYaxis()->SetTitle("Energy [MeV]");
    h1_califa_gamma_all[0]->GetXaxis()->SetTitle("Energy [MeV]");
    h1_califa_gamma_all[0]->GetYaxis()->SetTitle("Counts");
    h1_califa_gamma_25F_p2p_24O[0]->GetXaxis()->SetTitle("Energy [MeV]");
    h1_califa_gamma_25F_p2p_24O[0]->GetYaxis()->SetTitle("Counts");
    h1_califa_gamma_25F_p2p_23O[0]->GetXaxis()->SetTitle("Energy [MeV]");
    h1_califa_gamma_25F_p2p_23O[0]->GetYaxis()->SetTitle("Counts");
    h1_califa_gamma_25F_ppn_24F[0]->GetXaxis()->SetTitle("Energy [MeV]");
    h1_califa_gamma_25F_ppn_24F[0]->GetYaxis()->SetTitle("Counts");

    auto h2_califa_gamma_corr_all = CreateHistograms<TH2F>("h2_califa_gamma_corr_all",1,2000,0,100,200,0,10);
    auto h2_califa_gamma_corr_25F_p2p_24O = CreateHistograms<TH2F>("h2_califa_gamma_corr_25F_p2p_24O",1,2000,0,100,200,0,10);
    auto h2_califa_gamma_corr_25F_p2p_23O = CreateHistograms<TH2F>("h2_califa_gamma_corr_25F_p2p_23O",1,2000,0,100,200,0,10);
    auto h1_gamma_521_gated = CreateHistograms<TH1F>("h1_gamma_521_gated",1,200,0,10);
    auto h1_gamma_3118_gated = CreateHistograms<TH1F>("h1_gamma_3118_gated",1,500,0,5);
    auto h1_gamma_Esum_mult2 = CreateHistograms<TH1F>("h1_gamma_Esum_mult2",1,200,0,10);

    

    auto h2_califa_gamma_corr_25F_ppn_24F = CreateHistograms<TH2F>("h2_califa_gamma_corr_25F_ppn_24F",1,2000,0,100,200,0,10);
    auto h1_califa_gamma_corr_all = CreateHistograms<TH1F>("h1_califa_gamma_corr_all",1,200,0,10);
    auto h1_califa_gamma_corr_25F_p2p_24O = CreateHistograms<TH1F>("h1_califa_gamma_corr_25F_p2p_24O",1,200,0,10);
    auto h1_califa_gamma_corr_25F_p2p_23O = CreateHistograms<TH1F>("h1_califa_gamma_corr_25F_p2p_23O",1,200,0,10);
    auto h2_califa_gamma_corr_25F_p2p_22O = CreateHistograms<TH2F>("h2_califa_gamma_corr_25F_p2p_22O",1,2000,0,100,200,0,10);
    auto h1_califa_gamma_corr_25F_p2p_22O = CreateHistograms<TH1F>("h1_califa_gamma_corr_25F_p2p_22O",1,200,0,10);
    auto h2_califa_gamma_corr_25F_p2p_22O_mleq3 = CreateHistograms<TH2F>("h2_califa_gamma_corr_25F_p2p_22O_,leq3",1,2000,0,100,200,0,10);
    
    // Set axis titles for corrected gamma histograms
    h2_califa_gamma_corr_all[0]->GetXaxis()->SetTitle("Theta [deg]");
    h2_califa_gamma_corr_all[0]->GetYaxis()->SetTitle("Energy [MeV]");
    h2_califa_gamma_corr_25F_p2p_24O[0]->GetXaxis()->SetTitle("Theta [deg]");
    h2_califa_gamma_corr_25F_p2p_24O[0]->GetYaxis()->SetTitle("Energy [MeV]");
    h2_califa_gamma_corr_25F_p2p_23O[0]->GetXaxis()->SetTitle("Theta [deg]");
    h2_califa_gamma_corr_25F_p2p_23O[0]->GetYaxis()->SetTitle("Energy [MeV]");
    h2_califa_gamma_corr_25F_ppn_24F[0]->GetXaxis()->SetTitle("Theta [deg]");
    h2_califa_gamma_corr_25F_ppn_24F[0]->GetYaxis()->SetTitle("Energy [MeV]");
    h1_califa_gamma_corr_all[0]->GetXaxis()->SetTitle("Energy [MeV]");
    h1_califa_gamma_corr_all[0]->GetYaxis()->SetTitle("Counts");
    h1_califa_gamma_corr_25F_p2p_24O[0]->GetXaxis()->SetTitle("Energy [MeV]");
    h1_califa_gamma_corr_25F_p2p_24O[0]->GetYaxis()->SetTitle("Counts");
    h1_califa_gamma_corr_25F_p2p_23O[0]->GetXaxis()->SetTitle("Energy [MeV]");
    h1_califa_gamma_corr_25F_p2p_23O[0]->GetYaxis()->SetTitle("Counts");
    h1_gamma_521_gated[0]->GetXaxis()->SetTitle("Energy [MeV]");
    h1_gamma_521_gated[0]->GetYaxis()->SetTitle("Counts");
    h1_gamma_3118_gated[0]->GetXaxis()->SetTitle("Energy [MeV]");
    h1_gamma_3118_gated[0]->GetYaxis()->SetTitle("Counts");
    h1_gamma_Esum_mult2[0]->GetXaxis()->SetTitle("Energy [MeV]");
    h1_gamma_Esum_mult2[0]->GetYaxis()->SetTitle("Counts");
    
    auto h1_califa_gamma_corr_25F_ppn_24F_allmult = CreateHistograms<TH1F>("h1_califa_gamma_corr_25F_ppn_24F_allmult",1,200,0,10);
    auto h1_califa_gamma_corr_25F_ppn_24F_mult1 = CreateHistograms<TH1F>("h1_califa_gamma_corr_25F_ppn_24F_mult1",1,200,0,10);
    auto h1_califa_gamma_corr_25F_ppn_24F_mult2 = CreateHistograms<TH1F>("h1_califa_gamma_corr_25F_ppn_24F_mult2",1,200,0,10);
    auto h1_califa_gamma_corr_25F_ppn_24F_mult3 = CreateHistograms<TH1F>("h1_califa_gamma_corr_25F_ppn_24F_mult3",1,200,0,10);
    auto h1_califa_gamma_corr_25F_ppn_24F_mult4 = CreateHistograms<TH1F>("h1_califa_gamma_corr_25F_ppn_24F_mult4",1,200,0,10);
    auto h1_califa_gamma_corr_25F_ppn_24F_mult5 = CreateHistograms<TH1F>("h1_califa_gamma_corr_25F_ppn_24F_mult5",1,200,0,10);

    // fragment-beta corrected histograms use the existing *_gamma_corr names

    // All-multiplicity (all gamma multiplicities) 22O spectrum - match 24F binning
    auto h1_califa_gamma_corr_25F_p2p_22O_allmult = CreateHistograms<TH1F>("h1_califa_gamma_corr_25F_p2p_22O_allmult",1,200,0,10);

    // Set axis titles for 22O gamma histograms (fragment-beta corrected)
    h2_califa_gamma_corr_25F_p2p_22O[0]->GetXaxis()->SetTitle("Theta [deg]");
    h2_califa_gamma_corr_25F_p2p_22O[0]->GetYaxis()->SetTitle("Energy [MeV]");
    h1_califa_gamma_corr_25F_p2p_22O[0]->GetXaxis()->SetTitle("Energy [MeV]");
    h1_califa_gamma_corr_25F_p2p_22O[0]->GetYaxis()->SetTitle("Counts");
    // Axis titles for 22O all-multiplicity histogram
    h1_califa_gamma_corr_25F_p2p_22O_allmult[0]->GetXaxis()->SetTitle("Energy [MeV]");
    h1_califa_gamma_corr_25F_p2p_22O_allmult[0]->GetYaxis()->SetTitle("Counts");

    auto h1_califa_gamma_mult_25F_p2p_22O = CreateHistograms<TH1F>("h1_califa_gamma_mult_25F_p2p_22O", 1, 16, -0.5, 15.5);
    auto h1_califa_gamma_corr_25F_p2p_22O_mult_exact = CreateHistograms<TH1F>("h1_califa_gamma_corr_25F_p2p_22O_mult_exact", 10, 200, 0, 10);
    auto h1_califa_gamma_corr_25F_p2p_22O_mult_leq = CreateHistograms<TH1F>("h1_califa_gamma_corr_25F_p2p_22O_mult_leq", 10, 200, 0, 10);
    // 1D histograms counting gamma energies but only for events that produce at least one pair (i.e. multiplicity>=2)
    auto h1_califa_gamma_corr_25F_p2p_22O_mult_leq_pairs_only = CreateHistograms<TH1F>("h1_califa_gamma_corr_25F_p2p_22O_mult_leq_pairs_only", 10, 200, 0, 10);
    auto h1_frag_beta = CreateHistograms<TH1F>("h1_frag_beta", 1, 1000, 0.78, 0.84);
    auto h2_proton_cluster_energy_corr_25F_p2p_22O = CreateHistograms<TH2F>("h2_proton_cluster_energy_corr_25F_p2p_22O", 1, 400, 0, 400, 400, 0, 400);
    auto h2_gamma_theta_vs_proton1_theta_25F_p2p_22O = CreateHistograms<TH2F>("h2_gamma_theta_vs_proton1_theta_25F_p2p_22O", 1, 180, 0, 180, 180, 0, 180);
    auto h2_gamma_theta_vs_proton2_theta_25F_p2p_22O = CreateHistograms<TH2F>("h2_gamma_theta_vs_proton2_theta_25F_p2p_22O", 1, 180, 0, 180, 180, 0, 180);
    auto h1_gamma_gate_3199_corr_remaining_corr_25F_p2p_22O = CreateHistograms<TH1F>("h1_gamma_gate_3199_corr_remaining_corr_25F_p2p_22O", 1, 200, 0, 10);
    auto h1_gamma_gate_3199_corr_remaining_raw_25F_p2p_22O = CreateHistograms<TH1F>("h1_gamma_gate_3199_corr_remaining_raw_25F_p2p_22O", 1, 200, 0, 10);
    auto h1_gamma_gate_1383_corr_remaining_corr_25F_p2p_22O = CreateHistograms<TH1F>("h1_gamma_gate_1383_corr_remaining_corr_25F_p2p_22O", 1, 200, 0, 10);
    auto h1_gamma_gate_1383_corr_remaining_raw_25F_p2p_22O = CreateHistograms<TH1F>("h1_gamma_gate_1383_corr_remaining_raw_25F_p2p_22O", 1, 200, 0, 10);

    // Gated spectra per multiplicity condition (Mleq1..Mleq10) for 3.199 MeV gate
    auto h1_gamma_gate_3199_corr_remaining_corr_25F_p2p_22O_mult_leq = CreateHistograms<TH1F>("h1_gamma_gate_3199_corr_remaining_corr_25F_p2p_22O_mult_leq", 10, 200, 0, 10);
    auto h1_gamma_gate_3199_corr_remaining_raw_25F_p2p_22O_mult_leq = CreateHistograms<TH1F>("h1_gamma_gate_3199_corr_remaining_raw_25F_p2p_22O_mult_leq", 10, 200, 0, 10);
    // Gated spectra per multiplicity condition (Mleq1..Mleq10) for 1.383 MeV gate
    auto h1_gamma_gate_1383_corr_remaining_corr_25F_p2p_22O_mult_leq = CreateHistograms<TH1F>("h1_gamma_gate_1383_corr_remaining_corr_25F_p2p_22O_mult_leq", 10, 200, 0, 10);
    auto h1_gamma_gate_1383_corr_remaining_raw_25F_p2p_22O_mult_leq = CreateHistograms<TH1F>("h1_gamma_gate_1383_corr_remaining_raw_25F_p2p_22O_mult_leq", 10, 200, 0, 10);
    h1_gamma_gate_3199_corr_remaining_corr_25F_p2p_22O[0]->GetXaxis()->SetTitle("Doppler-corrected energy [MeV]");
    h1_gamma_gate_3199_corr_remaining_corr_25F_p2p_22O[0]->GetYaxis()->SetTitle("Counts");
    h1_gamma_gate_3199_corr_remaining_raw_25F_p2p_22O[0]->GetXaxis()->SetTitle("Uncorrected energy [MeV]");
    h1_gamma_gate_3199_corr_remaining_raw_25F_p2p_22O[0]->GetYaxis()->SetTitle("Counts");
    // Axis titles for multiplicity-gated versions
    for (size_t ii = 0; ii < h1_gamma_gate_3199_corr_remaining_corr_25F_p2p_22O_mult_leq.size(); ++ii) {
        h1_gamma_gate_3199_corr_remaining_corr_25F_p2p_22O_mult_leq[ii]->GetXaxis()->SetTitle("Doppler-corrected energy [MeV]");
        h1_gamma_gate_3199_corr_remaining_corr_25F_p2p_22O_mult_leq[ii]->GetYaxis()->SetTitle("Counts");
        h1_gamma_gate_3199_corr_remaining_raw_25F_p2p_22O_mult_leq[ii]->GetXaxis()->SetTitle("Uncorrected energy [MeV]");
        h1_gamma_gate_3199_corr_remaining_raw_25F_p2p_22O_mult_leq[ii]->GetYaxis()->SetTitle("Counts");
        h1_gamma_gate_1383_corr_remaining_corr_25F_p2p_22O_mult_leq[ii]->GetXaxis()->SetTitle("Doppler-corrected energy [MeV]");
        h1_gamma_gate_1383_corr_remaining_corr_25F_p2p_22O_mult_leq[ii]->GetYaxis()->SetTitle("Counts");
        h1_gamma_gate_1383_corr_remaining_raw_25F_p2p_22O_mult_leq[ii]->GetXaxis()->SetTitle("Uncorrected energy [MeV]");
        h1_gamma_gate_1383_corr_remaining_raw_25F_p2p_22O_mult_leq[ii]->GetYaxis()->SetTitle("Counts");
    }
    h1_gamma_gate_1383_corr_remaining_corr_25F_p2p_22O[0]->GetXaxis()->SetTitle("Doppler-corrected energy [MeV]");
    h1_gamma_gate_1383_corr_remaining_corr_25F_p2p_22O[0]->GetYaxis()->SetTitle("Counts");
    h1_gamma_gate_1383_corr_remaining_raw_25F_p2p_22O[0]->GetXaxis()->SetTitle("Uncorrected energy [MeV]");
    h1_gamma_gate_1383_corr_remaining_raw_25F_p2p_22O[0]->GetYaxis()->SetTitle("Counts");
    h1_frag_beta[0]->GetXaxis()->SetTitle("Fragment beta");
    h1_frag_beta[0]->GetYaxis()->SetTitle("Counts");
    auto h2_gamma_Epair_Meq2_25F_p2p_22O = CreateHistograms<TH2F>("h2_gamma_Epair_Meq2_25F_p2p_22O", 1, 400, 0, 4, 400, 0, 4);
    auto h2_gamma_Tpair_Meq2_25F_p2p_22O = CreateHistograms<TH2F>("h2_gamma_Tpair_Meq2_25F_p2p_22O", 1, 400, -200, 200, 400, -200, 200);
    auto h2_gamma_thetapair_Meq2_25F_p2p_22O = CreateHistograms<TH2F>("h2_gamma_thetapair_Meq2_25F_p2p_22O", 1, 180, 0, 180, 180, 0, 180);
    auto h2_gamma_phipair_Meq2_25F_p2p_22O = CreateHistograms<TH2F>("h2_gamma_phipair_Meq2_25F_p2p_22O", 1, 360, -180, 180, 360, -180, 180);
    // axis labels for Meq2 coincidence plots
    h2_gamma_Epair_Meq2_25F_p2p_22O[0]->GetXaxis()->SetTitle("Lower energy [MeV]");
    h2_gamma_Epair_Meq2_25F_p2p_22O[0]->GetYaxis()->SetTitle("Higher energy [MeV]");
    h2_gamma_Tpair_Meq2_25F_p2p_22O[0]->GetXaxis()->SetTitle("Lower time [ns]");
    h2_gamma_Tpair_Meq2_25F_p2p_22O[0]->GetYaxis()->SetTitle("Higher time [ns]");
    h2_gamma_thetapair_Meq2_25F_p2p_22O[0]->GetXaxis()->SetTitle("Lower theta [deg]");
    h2_gamma_thetapair_Meq2_25F_p2p_22O[0]->GetYaxis()->SetTitle("Higher theta [deg]");
    h2_gamma_phipair_Meq2_25F_p2p_22O[0]->GetXaxis()->SetTitle("Lower phi [deg]");
    h2_gamma_phipair_Meq2_25F_p2p_22O[0]->GetYaxis()->SetTitle("Higher phi [deg]");
    auto h2_gamma_Epair_Mleq2_25F_p2p_22O = CreateHistograms<TH2F>("h2_gamma_Epair_Mleq2_25F_p2p_22O", 1, 200, 0, 10, 200, 0, 10);
    auto h2_gamma_Tpair_Mleq2_25F_p2p_22O = CreateHistograms<TH2F>("h2_gamma_Tpair_Mleq2_25F_p2p_22O", 1, 400, -200, 200, 400, -200, 200);
    auto h2_gamma_Epair_Mleq3_25F_p2p_22O = CreateHistograms<TH2F>("h2_gamma_Epair_Mleq3_25F_p2p_22O", 3, 400, 0, 4, 400, 0, 4);
    auto h2_gamma_Tpair_Mleq3_25F_p2p_22O = CreateHistograms<TH2F>("h2_gamma_Tpair_Mleq3_25F_p2p_22O", 3, 400, -200, 200, 400, -200, 200);
    auto h2_gamma_thetapair_Mleq3_25F_p2p_22O = CreateHistograms<TH2F>("h2_gamma_thetapair_Mleq3_25F_p2p_22O", 3, 180, 0, 180, 180, 0, 180);
    auto h2_gamma_phipair_Mleq3_25F_p2p_22O = CreateHistograms<TH2F>("h2_gamma_phipair_Mleq3_25F_p2p_22O", 3, 360, -180, 180, 360, -180, 180);
    h2_gamma_Epair_Mleq2_25F_p2p_22O[0]->GetXaxis()->SetTitle("Lower energy [MeV]");
    h2_gamma_Epair_Mleq2_25F_p2p_22O[0]->GetYaxis()->SetTitle("Higher energy [MeV]");
    h2_gamma_Tpair_Mleq2_25F_p2p_22O[0]->GetXaxis()->SetTitle("Lower time [ns]");
    h2_gamma_Tpair_Mleq2_25F_p2p_22O[0]->GetYaxis()->SetTitle("Higher time [ns]");
    for (size_t ii=0; ii<3; ++ii) { h2_gamma_Epair_Mleq3_25F_p2p_22O[ii]->GetXaxis()->SetTitle("Lower energy [MeV]"); h2_gamma_Epair_Mleq3_25F_p2p_22O[ii]->GetYaxis()->SetTitle("Higher energy [MeV]"); }
    auto h2_gamma_Epair_Mleq4_25F_p2p_22O = CreateHistograms<TH2F>("h2_gamma_Epair_Mleq4_25F_p2p_22O", 6, 200, 0, 10, 200, 0, 10);
    auto h2_gamma_Tpair_Mleq4_25F_p2p_22O = CreateHistograms<TH2F>("h2_gamma_Tpair_Mleq4_25F_p2p_22O", 6, 400, -200, 200, 400, -200, 200);
    auto h2_gamma_thetapair_Mleq4_25F_p2p_22O = CreateHistograms<TH2F>("h2_gamma_thetapair_Mleq4_25F_p2p_22O", 6, 180, 0, 180, 180, 0, 180);
    auto h2_gamma_phipair_Mleq4_25F_p2p_22O = CreateHistograms<TH2F>("h2_gamma_phipair_Mleq4_25F_p2p_22O", 6, 360, -180, 180, 360, -180, 180);
    auto h2_gamma_Epair_Mleq5_25F_p2p_22O = CreateHistograms<TH2F>("h2_gamma_Epair_Mleq5_25F_p2p_22O", 10, 200, 0, 10, 200, 0, 10);
    for (size_t ii=0; ii<6; ++ii) { h2_gamma_Epair_Mleq4_25F_p2p_22O[ii]->GetXaxis()->SetTitle("Lower energy [MeV]"); h2_gamma_Epair_Mleq4_25F_p2p_22O[ii]->GetYaxis()->SetTitle("Higher energy [MeV]"); }
    for (size_t ii=0; ii<10; ++ii) { h2_gamma_Epair_Mleq5_25F_p2p_22O[ii]->GetXaxis()->SetTitle("Lower energy [MeV]"); h2_gamma_Epair_Mleq5_25F_p2p_22O[ii]->GetYaxis()->SetTitle("Higher energy [MeV]"); }
    
    // Rename multiplicity-binned histograms so ROOT names reflect actual multiplicity numbers (1..N)
    for (size_t i=0; i < h1_califa_gamma_corr_25F_p2p_24O_mult_exact.size(); ++i) {
        h1_califa_gamma_corr_25F_p2p_24O_mult_exact[i]->SetName(Form("h1_califa_gamma_corr_25F_p2p_24O_mult_exact_%d", int(i+1)));
    }
    for (size_t i=0; i < h1_califa_gamma_corr_25F_p2p_24O_mult_leq.size(); ++i) {
        h1_califa_gamma_corr_25F_p2p_24O_mult_leq[i]->SetName(Form("h1_califa_gamma_corr_25F_p2p_24O_mult_leq_%d", int(i+1)));
    }
    for (size_t i=0; i < h1_califa_gamma_corr_25F_p2p_22O_mult_exact.size(); ++i) {
        h1_califa_gamma_corr_25F_p2p_22O_mult_exact[i]->SetName(Form("h1_califa_gamma_corr_25F_p2p_22O_mult_exact_%d", int(i+1)));
    }
    for (size_t i=0; i < h1_califa_gamma_corr_25F_p2p_22O_mult_leq.size(); ++i) {
        h1_califa_gamma_corr_25F_p2p_22O_mult_leq[i]->SetName(Form("h1_califa_gamma_corr_25F_p2p_22O_mult_leq_%d", int(i+1)));
    }
    for (size_t i=0; i < h1_califa_gamma_corr_25F_p2p_22O_mult_leq_pairs_only.size(); ++i) {
        h1_califa_gamma_corr_25F_p2p_22O_mult_leq_pairs_only[i]->SetName(Form("h1_califa_gamma_corr_25F_p2p_22O_mult_leq_pairs_only_%d", int(i+1)));
    }

    // Rename multiplicity-gated arrays for gate histograms to explicit multiplicity numbers
    for (size_t i=0; i < h1_gamma_gate_3199_corr_remaining_corr_25F_p2p_22O_mult_leq.size(); ++i) {
        h1_gamma_gate_3199_corr_remaining_corr_25F_p2p_22O_mult_leq[i]->SetName(Form("h1_gamma_gate_3199_corr_remaining_corr_25F_p2p_22O_mult_leq_%d", int(i+1)));
        h1_gamma_gate_3199_corr_remaining_raw_25F_p2p_22O_mult_leq[i]->SetName(Form("h1_gamma_gate_3199_corr_remaining_raw_25F_p2p_22O_mult_leq_%d", int(i+1)));
        h1_gamma_gate_1383_corr_remaining_corr_25F_p2p_22O_mult_leq[i]->SetName(Form("h1_gamma_gate_1383_corr_remaining_corr_25F_p2p_22O_mult_leq_%d", int(i+1)));
        h1_gamma_gate_1383_corr_remaining_raw_25F_p2p_22O_mult_leq[i]->SetName(Form("h1_gamma_gate_1383_corr_remaining_raw_25F_p2p_22O_mult_leq_%d", int(i+1)));
    }

    // Diagnostic histograms: number of pairs filled per event for multiplicity conditions
    auto h1_paircount_Mleq2 = CreateHistograms<TH1F>("h1_paircount_Mleq2", 1, 10, 0, 10);
    auto h1_paircount_Mleq3 = CreateHistograms<TH1F>("h1_paircount_Mleq3", 1, 20, 0, 20);
    auto h1_paircount_Mleq4 = CreateHistograms<TH1F>("h1_paircount_Mleq4", 1, 30, 0, 30);
    auto h1_paircount_Mleq5 = CreateHistograms<TH1F>("h1_paircount_Mleq5", 1, 40, 0, 40);
    auto h2_gamma_Tpair_Mleq5_25F_p2p_22O = CreateHistograms<TH2F>("h2_gamma_Tpair_Mleq5_25F_p2p_22O", 10, 400, -200, 200, 400, -200, 200);
    auto h2_gamma_thetapair_Mleq5_25F_p2p_22O = CreateHistograms<TH2F>("h2_gamma_thetapair_Mleq5_25F_p2p_22O", 10, 180, 0, 180, 180, 0, 180);
    auto h2_gamma_phipair_Mleq5_25F_p2p_22O = CreateHistograms<TH2F>("h2_gamma_phipair_Mleq5_25F_p2p_22O", 10, 360, -180, 180, 360, -180, 180);
    auto h2_vertex_XY = CreateHistograms<TH2F>("h2_vertex_XY", 1, 400,-100,100,400,-100,100);
    auto h2_vertex_XZ = CreateHistograms<TH2F>("h2_vertex_XZ", 1, 1000,-500,500,1000,-500,500);
    auto h2_vertex_YZ = CreateHistograms<TH2F>("h2_vertex_YZ", 1, 1000,-500,500,1000,-500,500);
    auto h1_DCA = CreateHistograms<TH1F>("h1_DCA", 1, 1000,0,10);
    auto h1_vertex_Z = CreateHistograms<TH1F>("h1_vertex_Z", 1, 1000,-500,500);

    // Beam angular divergence: projected angles of the fragment momentum in x-z and y-z planes.
    // Fill for all 25F-gated events so that projecting onto each axis gives the
    // sigma_x and sigma_y inputs for SetBeamAngularDivergence() in the CALIFA simulation.
    auto h2_frag_beam_divergence = CreateHistograms<TH2F>("h2_frag_beam_divergence", 1,
                                                          500, -50., 50.,   // theta_x [mrad]
                                                          500, -50., 50.);  // theta_y [mrad]
    h2_frag_beam_divergence[0]->GetXaxis()->SetTitle("#theta_{x} = atan(p_{x}/p_{z}) [mrad]");
    h2_frag_beam_divergence[0]->GetYaxis()->SetTitle("#theta_{y} = atan(p_{y}/p_{z}) [mrad]");
    
    // Set branch addresses
    R3BEventHeader* header = nullptr;
    tree->SetBranchAddress("EventHeader.", &header);

    /*TClonesArray* WR_s2_array = nullptr;
    tree->SetBranchAddress("WRS2Data", &WR_s2_array);

    TClonesArray* WR_los_array = nullptr;
    tree->SetBranchAddress("WRLosData", &WR_los_array);

    TClonesArray* WR_foot_array = nullptr;
    tree->SetBranchAddress("WRFootData", &WR_foot_array);

    TClonesArray* WR_alpide_array = nullptr;
    tree->SetBranchAddress("WRAlpideData", &WR_alpide_array);

    TClonesArray* WR_califa_array = nullptr;
    tree->SetBranchAddress("WRCalifaData", &WR_califa_array);
    */
    TClonesArray* frs_array = nullptr;
    tree->SetBranchAddress("FrsData", &frs_array);

    TClonesArray* tofd_array = nullptr;
    tree->SetBranchAddress("TofdHit", &tofd_array);

    TClonesArray* neuland_hit_array = nullptr;
    tree->SetBranchAddress("NeulandHits", &neuland_hit_array);

    TClonesArray* foot_array = nullptr;
    tree->SetBranchAddress("FootHitData", &foot_array);

    TClonesArray* fragment_tracks = nullptr;
    tree->SetBranchAddress("FragmentMDFTrack", &fragment_tracks);

    TClonesArray* incoming_tracks = nullptr;
    tree->SetBranchAddress("IncomingTrackFoot", &incoming_tracks);

    TClonesArray* califa_data = nullptr;
    tree->SetBranchAddress("CalifaClusterData", &califa_data);

    Long64_t nEntries = tree->GetEntries();

    std::map<double, TVector3> map_neuland_positions;

    auto fillSortedGammaPair = [&](TH2* hE, TH2* hT, TH2* hTheta, TH2* hPhi,
                                   double e1, double e2,
                                   double t1, double t2,
                                   double th1, double th2,
                                   double ph1, double ph2)
    {
        // Ensure lower energy on X and higher energy on Y, along with their associated quantities
        double ex = e1, ey = e2;
        double tx = t1, ty = t2;
        double thx = th1, thy = th2;
        double phx = ph1, phy = ph2;
        if(ex > ey) { std::swap(ex, ey); std::swap(tx, ty); std::swap(thx, thy); std::swap(phx, phy); }
        if(hE) hE->Fill(ex, ey);
        if(hT) hT->Fill(tx, ty);
        if(hTheta) hTheta->Fill(thx, thy);
        if(hPhi) hPhi->Fill(phx, phy);
    };

    auto fillGammaPairFamily = [&](const std::vector<double>& eCorr, const std::vector<double>& times,
                                   const std::vector<double>& thetas, const std::vector<double>& phis,
                                   size_t maxMult,
                                   std::vector<TH2F*>& hE, std::vector<TH2F*>& hT,
                                   std::vector<TH2F*>& hTheta, std::vector<TH2F*>& hPhi)
    {
        if(eCorr.size() < 2 || eCorr.size() > maxMult) return;
        size_t idx = 0;
        for(size_t i = 0; i < eCorr.size(); ++i)
        {
            for(size_t j = i + 1; j < eCorr.size(); ++j)
            {
                if(idx >= hE.size() || idx >= hT.size() || idx >= hTheta.size() || idx >= hPhi.size()) return;
                fillSortedGammaPair(hE[idx], hT[idx], hTheta[idx], hPhi[idx],
                                    eCorr[i], eCorr[j], times[i], times[j], thetas[i], thetas[j], phis[i], phis[j]);
                ++idx;
            }
        }
    };

    for(Long64_t i = 0; i < nEntries; ++i)
    {
        // Clear arrays
        /*if(WR_s2_array) WR_s2_array->Clear();
        if(WR_los_array) WR_los_array->Clear();
        if(WR_foot_array) WR_foot_array->Clear();
        if(WR_alpide_array) WR_alpide_array->Clear();
        if(WR_califa_array) WR_califa_array->Clear();
        */
	if(tofd_array) tofd_array->Clear();
        if(foot_array) foot_array->Clear();
        if(frs_array) frs_array->Clear();
        if(fragment_tracks) fragment_tracks->Clear();
        if(incoming_tracks) incoming_tracks->Clear();
        if(califa_data) califa_data->Clear();

        tree->GetEntry(i);
        
        if(i % 100000 == 0) 
            cout << "\r Number of events processed: " << i << flush;

        // Early cuts to skip expensive processing
        if(header->GetTrigger() != 1) continue;
        if(header->GetTpat() > 32) continue;
        if(incoming_tracks->GetEntriesFast() == 0) continue;
        if(fragment_tracks->GetEntriesFast() == 0) continue;
	
        
	// Get FRS data
        //uint64_t ts_main = header->GetTimeStamp();
        //auto ts_foot = extractFromTCA<R3BWRData>(WR_foot_array, &R3BWRData::GetTimeStamp);
        //auto ts_alpide = extractFromTCA<R3BWRData>(WR_alpide_array, &R3BWRData::GetTimeStamp);
        //auto ts_califa = extractFromTCA<R3BWRData>(WR_califa_array, &R3BWRData::GetTimeStamp);
        //auto ts_s2 = extractFromTCA<R3BWRData>(WR_s2_array, &R3BWRData::GetTimeStamp);
        //auto ts_los = extractFromTCA<R3BWRData>(WR_los_array, &R3BWRData::GetTimeStamp);

        // Fill WR histograms with bounds checking
      /*  for(size_t j = 0; j < 16 && j < ts_foot.size(); ++j){
            h1_dwr_main_foot[j]->Fill((long)ts_foot[j] - long(ts_main));
        }
        for(size_t j = 0; j < 4 && j < ts_califa.size(); ++j){
            h1_dwr_main_califa[j]->Fill((long)ts_main - (long)ts_califa[j]);
        }
        //for(size_t j = 0; j < 2 && j < ts_alpide.size(); ++j){
        //    h1_dwr_main_alpide[j]->Fill((long)ts_main - (long)ts_alpide[j]);
        //}
        if(ts_s2.size() > 0)
            h1_dwr_main_s2[0]->Fill((long)ts_main - (long)ts_s2[0]);
        if(ts_los.size() > 0)
            h1_dwr_main_los[0]->Fill((long)ts_main - (long)ts_los[0]);
	*/

        // Get FRS data
        double frs_Z = 0; 
        double frs_AoZ = 0;
        double frs_Beta = 0;
        double frs_Gamma = 0;
        
        for(int j = 0; j < frs_array->GetEntriesFast(); ++j)
        {
            auto frs_item = dynamic_cast<R3BFrsData*>(frs_array->At(j));
          if(frs_item->GetStaId() != 2 && frs_item->GetStaId() !=1) continue;
            frs_Z = frs_item->GetZ();
            frs_AoZ = frs_item->GetAq();
            frs_Beta = frs_item->GetBeta();
            frs_Gamma = 1./sqrt(1. - pow(frs_Beta, 2));
        } 
        
        if(frs_Z == 0 || frs_AoZ == 0) continue;
        h2_frs_pid[0]->Fill(frs_AoZ, frs_Z);

        auto fragment = dynamic_cast<R3BTrackingParticle*>(fragment_tracks->At(0));
        auto intrack = dynamic_cast<R3BTrackingParticle*>(incoming_tracks->At(0));
        if(!fragment || !intrack) continue;

        const auto out_track_mom = fragment->GetStartMomentum();
        const auto out_track_pos = fragment->GetStartPosition();
        const auto frag_Z = fragment->GetCharge();
        const auto frag_AoQ = fragment->GetMass();
        const auto PoQ = out_track_mom.Mag();
        double frag_beta = fragment->GetBeta();
        if(frag_beta <= 0.0 || frag_beta >= 1.0)
        {
            double pFrag = PoQ * std::max(1.0, frag_Z);
            double mFrag = std::max(0.0, frag_AoQ * std::max(1.0, frag_Z));
            if(pFrag > 0.0 && mFrag > 0.0)
            {
                double eFrag = std::sqrt(pFrag * pFrag + mFrag * mFrag);
                frag_beta = pFrag / eFrag;
            }
        }
        double frag_gamma = (frag_beta > 0.0 && frag_beta < 1.0) ? 1./std::sqrt(1. - frag_beta * frag_beta) : 1.0;
        const auto in_track_pos = intrack->GetStartPosition();
        const auto in_track_mom = intrack->GetStartMomentum();
        const auto in_track_mom_unit = in_track_mom.Unit();
        const auto out_track_mom_unit = out_track_mom.Unit();

        // Get CALIFA data
        vector<R3BCalifaClusterData*> califa_hits;
        vector<R3BCalifaClusterData*> califa_gammas;
        vector<double> califa_theta_g;
	vector<double> califa_phi_g;
        vector<double> califa_time_g;
        // record fragment beta for diagnostics
        h1_frag_beta[0]->Fill(frag_beta);
        vector<double> califa_e_g;
        vector<double> califa_e_g_corr;
	vector<double> deltaphi_proton_gamma;  //phi difference between a gamma cluster and the proton

         // Cut for incoming 25F
        bool is25F = (frs_AoZ > 2.765 && frs_AoZ < 2.79) && (frs_Z > 8 && frs_Z < 10.3); //correct cut !

        // Beam angular divergence: for 25F events fill the projected fragment angles.
        // out_track_mom is the MDF-reconstructed fragment momentum (TVector3).
        // atan2(px, pz) and atan2(py, pz) give the projected angles in the x-z and y-z planes.
        if(is25F && out_track_mom.Z() != 0.)
        {
            const double theta_x_mrad = 1000. * std::atan2(out_track_mom.X(), out_track_mom.Z());
            const double theta_y_mrad = 1000. * std::atan2(out_track_mom.Y(), out_track_mom.Z());
            h2_frag_beam_divergence[0]->Fill(theta_x_mrad, theta_y_mrad);
        }

        for(int j = 0; j < califa_data->GetEntriesFast(); ++j)
        {
            auto califa_item = dynamic_cast<R3BCalifaClusterData*>(califa_data->At(j));
            if(is25F){
            
	    if(califa_item->GetEnergy()>1000*20) //incoming 25F and cluster > 20 MeV, prefilter for p/n hits
            {
                califa_hits.push_back(califa_item);
            }
            else if(califa_item->/*GetClusterType()==1*/GetEnergy()<1000*10) //gamma clusters
            {
                califa_gammas.push_back(califa_item);
                califa_e_g.push_back(califa_item->GetEnergy()/1000.);
                califa_theta_g.push_back(califa_item->GetTheta()*TMath::RadToDeg());
		califa_phi_g.push_back(califa_item->GetPhi()*TMath::RadToDeg());
                // corrected energies will be computed later using fragment beta
            }
	    }
        }

        bool isatleasttwohits = (califa_hits.size()>=2); //relaxed conditions
        
	if(!isatleasttwohits) continue; 
        
        //function to sort the cluster vector in descending order of energies (first one would typically be a proton cluster)
	std::sort(califa_hits.begin(), califa_hits.end(),
          [](const R3BCalifaClusterData* a,
             const R3BCalifaClusterData* b)
          {
              return a->GetEnergy() > b->GetEnergy();
          });
	
	//same for gamma clusters
	std::sort(califa_gammas.begin(), califa_gammas.end(),
          [](const R3BCalifaClusterData* a,
             const R3BCalifaClusterData* b)
          {
              return a->GetEnergy() > b->GetEnergy();
          });

	bool twohighenergyhits = (califa_hits[0]->GetEnergy()/1000. > 35. && califa_hits[1]->GetEnergy()/1000. > 35.); // colleague-style selected proton thresholds
	bool onehighenergyhit = (califa_hits[0]->GetEnergy()/1000. > 35.); // colleague-style selected proton threshold

        califa_e_g_corr.clear();
        califa_time_g.resize(califa_gammas.size());
        califa_theta_g.resize(califa_gammas.size());
        califa_phi_g.resize(califa_gammas.size());
        califa_e_g.resize(califa_gammas.size());
        for(size_t ig = 0; ig < califa_gammas.size(); ++ig)
        {
            auto* g = califa_gammas[ig];
            califa_e_g[ig] = g->GetEnergy()/1000.;
            califa_theta_g[ig] = g->GetTheta()*TMath::RadToDeg();
            califa_phi_g[ig] = g->GetPhi()*TMath::RadToDeg();
            califa_time_g[ig] = g->GetTime();
            // Compute the angle between the gamma direction and the actual fragment
            // momentum direction (out_track_mom_unit) from the MDF track.
            // This replaces g->GetTheta() which was the crystal angle from the z-axis,
            // ignoring the real beam tilt event-by-event.
            TVector3 gamma_dir;
            gamma_dir.SetMagThetaPhi(1.0, g->GetTheta(), g->GetPhi());
            const double theta_frag_gamma = gamma_dir.Angle(out_track_mom_unit); // [rad]
            califa_e_g_corr.push_back(DopplerCorrectGammaEnergy(califa_e_g[ig], theta_frag_gamma, frag_beta));
        }
	
	int it_1 = (gRandom->Uniform(-1,1) > 0) ? 1 : 0;
        int it_2 = (it_1 > 0) ? 0 : 1;
        
	double theta1 = califa_hits[it_1]->GetTheta();
	double phi1   = califa_hits[it_1]->GetPhi();
	double theta2 = califa_hits[it_2]->GetTheta();
	double phi2   = califa_hits[it_2]->GetPhi();

	TVector3 vec1, vec2;
	vec1.SetMagThetaPhi(1.0, theta1, phi1);
	vec2.SetMagThetaPhi(1.0, theta2, phi2);

	double opang = vec1.Angle(vec2) * TMath::RadToDeg();
        double deltaphi =  TMath::Abs(phi1 - phi2) * TMath::RadToDeg();
        double deltaphifolded = deltaphi;
	double DeltaPhiVeto = 30.;
	if (deltaphi > 180.) deltaphifolded = 360. - deltaphi;

	//fill general (no cuts yet) histograms
	h2_califa_theta_corr[0]->Fill(theta1*TMath::RadToDeg(), theta2*TMath::RadToDeg());
	h1_califa_opang[0]->Fill(opang);  
	h2_califa_Ecorr[0]->Fill(califa_hits[it_1]->GetEnergy()/1000., califa_hits[it_2]->GetEnergy()/1000.);
        h1_califa_deltaphi[0]->Fill(deltaphi);
	h1_califa_deltaphifolded[0]->Fill(deltaphifolded);

	       
       	bool isQFS = false;
        if (opang > 77.34 && opang < 94.53 && TMath::Abs(180. - deltaphifolded) < DeltaPhiVeto) isQFS = true;
        
	//Loop over all clusters in a QFS event and fill in the Esum_vs_theta (of proton) correlation plot
	if(isQFS){
	double sumenergy = 0;	
	for(int i = 0; i< califa_hits.size(); ++i)
	{
		sumenergy+=califa_hits[i]->GetEnergy()/1000.;
		
	}
        h2_califa_E_vs_theta[0]->Fill(califa_hits[0]->GetTheta()*TMath::RadToDeg(),sumenergy);
	}
        

        const int foot5_idx = FindFootHitIndexByDetId(foot_array, 5);
        double foot_charge = -999.;
        if(foot5_idx >= 0)
        {
            auto* foot_item = dynamic_cast<R3BFootHitData*>(foot_array->At(foot5_idx));
            if(foot_item) foot_charge = foot_item->GetZCharge();
        }
        
	        auto target_vertex = get_vertex(in_track_pos, out_track_pos, in_track_mom_unit, out_track_mom_unit);
        
        if(target_vertex.dca < 0.2)
        {
            if(target_vertex.vertexPos.Z() > (-25) && target_vertex.vertexPos.Z() < 120)
            {
                h2_vertex_XY[0]->Fill(target_vertex.vertexPos.X(), target_vertex.vertexPos.Y());
            }
            h2_vertex_XZ[0]->Fill(target_vertex.vertexPos.Z(), target_vertex.vertexPos.X());
            h2_vertex_YZ[0]->Fill(target_vertex.vertexPos.Z(), target_vertex.vertexPos.Y());
            h1_vertex_Z[0]->Fill(target_vertex.vertexPos.Z());
        }
       
        //CALIFA with specific PID
	bool is_p2p_24O =  (frag_Z > 7.6 && frag_Z < 8.5) && (frag_AoQ > 2.9 && frag_AoQ < 3.02);
    bool is_p2p_23O =  (frag_Z > 7.6 && frag_Z < 8.5) && (frag_AoQ > 2.79 && frag_AoQ < 2.88);
	bool is_p2p_22O =  (frag_Z > 7.6 && frag_Z < 8.5) && (frag_AoQ > 2.67 && frag_AoQ < 2.76);
	bool is_ppn_24F =  (frag_Z > 8.5 && frag_Z < 9.5) && (frag_AoQ > 2.57 && frag_AoQ < 2.67);
	bool is_ppn_23F =  (frag_Z > 8.5 && frag_Z < 9.5) && (frag_AoQ > 2.45 && frag_AoQ < 2.56);
        
              
	    
		for(size_t a = 0; a < califa_gammas.size(); ++a){   //a loop over all califa gamma clusters
            h2_califa_gamma_all[0]->Fill(califa_theta_g.at(a), califa_e_g.at(a));
	        h1_califa_gamma_all[0]->Fill(califa_e_g.at(a));
		    h2_califa_gamma_corr_all[0]->Fill(califa_theta_g.at(a), califa_e_g_corr.at(a));
            h1_califa_gamma_corr_all[0]->Fill(califa_e_g_corr.at(a));
	        double absdeltaphi_p_g = std::abs(califa_phi_g.at(a)-(califa_hits[0]->GetPhi()*TMath::RadToDeg()));
            if(absdeltaphi_p_g > 180) absdeltaphi_p_g = 360 - absdeltaphi_p_g; //fold back
		        deltaphi_proton_gamma.push_back(absdeltaphi_p_g);
		
		}
	    
	size_t nGamma = std::min({
			califa_gammas.size(),
			califa_theta_g.size(),
			califa_phi_g.size(),
			califa_e_g.size(),
			califa_e_g_corr.size(),
			deltaphi_proton_gamma.size()});  //find min size of vectors and then resize them to min, so as to avoid out of bound problems
	
	auto shrink = [&](auto& v) {v.resize(nGamma);}; //lambda function
/*	shrink(califa_gammas);
	shrink(califa_theta_g);
	shrink(califa_phi_g);
	shrink(califa_e_g);
	shrink(califa_e_g_corr);
	shrink(deltaphi_proton_gamma);
*/		
	if(is_p2p_24O)
        {
        h2_califa_theta_corr_25F_p2p_24O[0]->Fill(theta1*TMath::RadToDeg(), theta2*TMath::RadToDeg());
        if (TMath::Abs(180. - deltaphifolded) < DeltaPhiVeto && twohighenergyhits) 
            {
            h1_califa_opang_25F_p2p_24O[0]->Fill(opang);
            }
            h2_califa_Ecorr_25F_p2p_24O[0]->Fill(califa_hits[it_1]->GetEnergy()/1000, califa_hits[it_2]->GetEnergy()/1000);
            if(isQFS && twohighenergyhits){
		        for(size_t a = 0; a < califa_gammas.size(); ++a){
                    h2_califa_gamma_25F_p2p_24O[0]->Fill(califa_theta_g.at(a), califa_e_g.at(a));
	                h1_califa_gamma_25F_p2p_24O[0]->Fill(califa_e_g.at(a));
		            h2_califa_gamma_corr_25F_p2p_24O[0]->Fill(califa_theta_g.at(a), califa_e_g_corr.at(a));
                    h1_califa_gamma_corr_25F_p2p_24O[0]->Fill(califa_e_g_corr.at(a));
			
    		    }

                // fill multiplicity-specific 24O histograms (exact and leq)
                const size_t gammaMult24 = califa_gammas.size();
                for(size_t m = 1; m <= 10; ++m)
                {
                    if(gammaMult24 <= m)
                    {
                        for(size_t a = 0; a < gammaMult24; ++a) h1_califa_gamma_corr_25F_p2p_24O_mult_leq[m-1]->Fill(califa_e_g_corr[a]);
                    }
                    if(gammaMult24 == m)
                    {
                        for(size_t a = 0; a < gammaMult24; ++a) h1_califa_gamma_corr_25F_p2p_24O_mult_exact[m-1]->Fill(califa_e_g_corr[a]);
                    }
                }
	    }
	    
        }

        if(is_p2p_23O)
        {
            h2_califa_theta_corr_25F_p2p_23O[0]->Fill(theta1*TMath::RadToDeg(), theta2*TMath::RadToDeg());
            // fill opening-angle when coplanarity (delta-phi) condition is satisfied
            if (TMath::Abs(180. - deltaphifolded) < DeltaPhiVeto && twohighenergyhits) h1_califa_opang_25F_p2p_23O[0]->Fill(opang);
            h2_califa_Ecorr_25F_p2p_23O[0]->Fill(califa_hits[it_1]->GetEnergy()/1000, califa_hits[it_2]->GetEnergy()/1000);
            if(isQFS && twohighenergyhits){
                for(size_t a = 0; a < califa_gammas.size(); ++a){
                    h2_califa_gamma_25F_p2p_23O[0]->Fill(califa_theta_g.at(a), califa_e_g.at(a));
                    h1_califa_gamma_25F_p2p_23O[0]->Fill(califa_e_g.at(a));
                    h2_califa_gamma_corr_25F_p2p_23O[0]->Fill(califa_theta_g.at(a), califa_e_g_corr.at(a));
                    h1_califa_gamma_corr_25F_p2p_23O[0]->Fill(califa_e_g_corr.at(a));
                }
            }

        }

        if(is_p2p_22O)
        {
            h2_proton_cluster_energy_corr_25F_p2p_22O[0]->Fill(califa_hits[1]->GetEnergy()/1000., califa_hits[0]->GetEnergy()/1000.);
            // fill opening-angle when delta-phi (coplanarity) is satisfied; keep rest of QFS logic unchanged
            if (TMath::Abs(180. - deltaphifolded) < DeltaPhiVeto && twohighenergyhits) h1_califa_opang_25F_p2p_22O[0]->Fill(opang);
            if(isQFS && twohighenergyhits)
            {
                const size_t gammaMult = califa_gammas.size();
                h1_califa_gamma_mult_25F_p2p_22O[0]->Fill(static_cast<double>(gammaMult));
                
		if(gammaMult <=3){    //A 2D plot of energy vs theta only for m<=3
			for(size_t a = 0; a < gammaMult; ++a)
                           {
                            h2_califa_gamma_corr_25F_p2p_22O_mleq3[0]->Fill(califa_theta_g[a], califa_e_g_corr[a]);
                           }
		}

		for(size_t a = 0; a < gammaMult; ++a)
                {
                    h2_califa_gamma_corr_25F_p2p_22O[0]->Fill(califa_theta_g[a], califa_e_g_corr[a]);
                    h1_califa_gamma_corr_25F_p2p_22O[0]->Fill(califa_e_g_corr[a]);
                    // Fill all-multiplicity 22O histogram (same binning as 24F allmult)
                    h1_califa_gamma_corr_25F_p2p_22O_allmult[0]->Fill(califa_e_g_corr[a]);
                    h2_gamma_theta_vs_proton1_theta_25F_p2p_22O[0]->Fill(califa_hits[0]->GetTheta()*TMath::RadToDeg(), califa_theta_g[a]);
                    h2_gamma_theta_vs_proton2_theta_25F_p2p_22O[0]->Fill(califa_hits[1]->GetTheta()*TMath::RadToDeg(), califa_theta_g[a]);
                }

                for(size_t m = 1; m <= 10; ++m)
                {
                    if(gammaMult <= m)
                    {
                        for(size_t a = 0; a < gammaMult; ++a) h1_califa_gamma_corr_25F_p2p_22O_mult_leq[m-1]->Fill(califa_e_g_corr[a]);
                        // also fill the "pairs-only" histogram for this multiplicity if event actually produces pairs
                        if(gammaMult >= 2) {
                            for(size_t a = 0; a < gammaMult; ++a) h1_califa_gamma_corr_25F_p2p_22O_mult_leq_pairs_only[m-1]->Fill(califa_e_g_corr[a]);
                        }
                    }
                    if(gammaMult == m)
                    {
                        for(size_t a = 0; a < gammaMult; ++a) h1_califa_gamma_corr_25F_p2p_22O_mult_exact[m-1]->Fill(califa_e_g_corr[a]);
                    }
                }

                constexpr double gate3199 = 3.199;
                constexpr double gate1383 = 1.383;
                constexpr double gateWidth = 0.12;
                for(size_t iGamma = 0; iGamma < gammaMult; ++iGamma)
                {
                    if(std::abs(califa_e_g_corr[iGamma] - gate3199) < gateWidth)
                    {
                        for(size_t jGamma = 0; jGamma < gammaMult; ++jGamma)
                        {
                            if(jGamma == iGamma) continue;
                                    h1_gamma_gate_3199_corr_remaining_corr_25F_p2p_22O[0]->Fill(califa_e_g_corr[jGamma]);
                                    h1_gamma_gate_3199_corr_remaining_raw_25F_p2p_22O[0]->Fill(califa_e_g[jGamma]);
                                    // also fill multiplicity-conditioned gated histograms (Mleq)
                                    for(size_t m = 1; m <= h1_gamma_gate_3199_corr_remaining_corr_25F_p2p_22O_mult_leq.size(); ++m) {
                                        if(gammaMult <= m) {
                                            h1_gamma_gate_3199_corr_remaining_corr_25F_p2p_22O_mult_leq[m-1]->Fill(califa_e_g_corr[jGamma]);
                                            h1_gamma_gate_3199_corr_remaining_raw_25F_p2p_22O_mult_leq[m-1]->Fill(califa_e_g[jGamma]);
                                        }
                                    }
                        }
                    }
                    if(std::abs(califa_e_g_corr[iGamma] - gate1383) < gateWidth)
                    {
                        for(size_t jGamma = 0; jGamma < gammaMult; ++jGamma)
                        {
                            if(jGamma == iGamma) continue;
                                    h1_gamma_gate_1383_corr_remaining_corr_25F_p2p_22O[0]->Fill(califa_e_g_corr[jGamma]);
                                    h1_gamma_gate_1383_corr_remaining_raw_25F_p2p_22O[0]->Fill(califa_e_g[jGamma]);
                                    // also fill multiplicity-conditioned gated histograms (Mleq)
                                    for(size_t m = 1; m <= h1_gamma_gate_1383_corr_remaining_corr_25F_p2p_22O_mult_leq.size(); ++m) {
                                        if(gammaMult <= m) {
                                            h1_gamma_gate_1383_corr_remaining_corr_25F_p2p_22O_mult_leq[m-1]->Fill(califa_e_g_corr[jGamma]);
                                            h1_gamma_gate_1383_corr_remaining_raw_25F_p2p_22O_mult_leq[m-1]->Fill(califa_e_g[jGamma]);
                                        }
                                    }
                        }
                    }
                }

                if(gammaMult == 2)
                {
                    fillSortedGammaPair(h2_gamma_Epair_Meq2_25F_p2p_22O[0], h2_gamma_Tpair_Meq2_25F_p2p_22O[0],
                                        h2_gamma_thetapair_Meq2_25F_p2p_22O[0], h2_gamma_phipair_Meq2_25F_p2p_22O[0],
                                        califa_e_g_corr[0], califa_e_g_corr[1], califa_time_g[0], califa_time_g[1],
                                        califa_theta_g[0], califa_theta_g[1], califa_phi_g[0], califa_phi_g[1]);
                }
                {
                    int pairCount = 0;
                    if(gammaMult <= 2) pairCount = static_cast<int>(gammaMult * (gammaMult - 1) / 2);
                    h1_paircount_Mleq2[0]->Fill(pairCount);
                }
                fillGammaPairFamily(califa_e_g_corr, califa_time_g, califa_theta_g, califa_phi_g, 2,
                                    h2_gamma_Epair_Mleq2_25F_p2p_22O, h2_gamma_Tpair_Mleq2_25F_p2p_22O,
                                    h2_gamma_thetapair_Meq2_25F_p2p_22O, h2_gamma_phipair_Meq2_25F_p2p_22O);
                {
                    int pairCount = 0;
                    if(gammaMult <= 3) pairCount = static_cast<int>(gammaMult * (gammaMult - 1) / 2);
                    h1_paircount_Mleq3[0]->Fill(pairCount);
                }
                fillGammaPairFamily(califa_e_g_corr, califa_time_g, califa_theta_g, califa_phi_g, 3,
                                    h2_gamma_Epair_Mleq3_25F_p2p_22O, h2_gamma_Tpair_Mleq3_25F_p2p_22O,
                                    h2_gamma_thetapair_Mleq3_25F_p2p_22O, h2_gamma_phipair_Mleq3_25F_p2p_22O);
                {
                    int pairCount = 0;
                    if(gammaMult <= 4) pairCount = static_cast<int>(gammaMult * (gammaMult - 1) / 2);
                    h1_paircount_Mleq4[0]->Fill(pairCount);
                }
                fillGammaPairFamily(califa_e_g_corr, califa_time_g, califa_theta_g, califa_phi_g, 4,
                                    h2_gamma_Epair_Mleq4_25F_p2p_22O, h2_gamma_Tpair_Mleq4_25F_p2p_22O,
                                    h2_gamma_thetapair_Mleq4_25F_p2p_22O, h2_gamma_phipair_Mleq4_25F_p2p_22O);
                {
                    int pairCount = 0;
                    if(gammaMult <= 5) pairCount = static_cast<int>(gammaMult * (gammaMult - 1) / 2);
                    h1_paircount_Mleq5[0]->Fill(pairCount);
                }
                fillGammaPairFamily(califa_e_g_corr, califa_time_g, califa_theta_g, califa_phi_g, 5,
                                    h2_gamma_Epair_Mleq5_25F_p2p_22O, h2_gamma_Tpair_Mleq5_25F_p2p_22O,
                                    h2_gamma_thetapair_Mleq5_25F_p2p_22O, h2_gamma_phipair_Mleq5_25F_p2p_22O);
            }
        }

        if(is_ppn_24F)
        {
            //h2_proton_cluster_energy_corr_25F_ppn_24F[0]->Fill(califa_hits[1]->GetEnergy()/1000., califa_hits[0]->GetEnergy()/1000.);
            if(isQFS && twohighenergyhits)
            {
                
                const size_t gammaMult = califa_gammas.size();
		for(size_t a = 0; a < gammaMult; ++a){
		h1_califa_gamma_25F_ppn_24F[0]->Fill(califa_e_g[a]);
                h2_califa_gamma_corr_25F_ppn_24F[0]->Fill(califa_theta_g[a], califa_e_g_corr[a]);
                h1_califa_gamma_corr_25F_ppn_24F_allmult[0]->Fill(califa_e_g_corr[a]);
	        if(califa_gammas.size()==1)h1_califa_gamma_corr_25F_ppn_24F_mult1[0]->Fill(califa_e_g_corr[a]); //gamma cluster multiplicity 1
                if(califa_gammas.size()==2)h1_califa_gamma_corr_25F_ppn_24F_mult2[0]->Fill(califa_e_g_corr[a]); //gamma cluster multiplicity 2
                if(califa_gammas.size()==3)h1_califa_gamma_corr_25F_ppn_24F_mult3[0]->Fill(califa_e_g_corr[a]); //gamma cluster multiplicity 3
                if(califa_gammas.size()==4)h1_califa_gamma_corr_25F_ppn_24F_mult4[0]->Fill(califa_e_g_corr[a]); //gamma cluster multiplicity 4
                if(califa_gammas.size()==5)h1_califa_gamma_corr_25F_ppn_24F_mult5[0]->Fill(califa_e_g_corr[a]); //gamma cluster multiplicity 5
		}
           	 }
		
        }
 
        // Excitation energy calculation
        double dalt = 931.494;
        double Sp_25F = 14.4552;
        double mp = 938.272;
        
            double sp = (frag_AoQ*frag_Z*dalt) + mp - (frs_Z*frs_AoZ*dalt);
            double Ex = sp - Sp_25F;
           
    	    h1_DCA[0]->Fill(target_vertex.dca);
        //h2_charge_foot_vs_tofd[0]->Fill(foot_charge, frag_Z);
        h2_AoZ_vs_Z[0]->Fill(frag_AoQ, frag_Z);
        //h2_tofd_pid[0]->Fill(bar_tofd, charge_tofd);

        double Mass = 0; 
        double Charge = 0;
        if(Init_A_and_Z(Mass, Charge, frag_AoQ, frag_Z, foot_charge)) 
        {
            double tof_mdf = fragment->GetLength() / SPEED_OF_LIGHT * sqrt(pow((Mass/Charge/PoQ), 2) + 1);
            /*if(array_th2_tofd[bar_tofd]->GetEntries() < 1e5)
            {
                array_th2_tofd[bar_tofd]->Fill(tof_tofd, tof_mdf);
                array_th1_tofd[bar_tofd]->Fill(tof_mdf - tof_tofd);
            }*/
        }

        // Neutron analysis
        auto neuland_positions = extractFromTCA<R3BNeulandHit>(neuland_hit_array, &R3BNeulandHit::GetPosition);
        if(neuland_positions.size() < 1) continue;
        auto neuland_times = extractFromTCA<R3BNeulandHit>(neuland_hit_array, &R3BNeulandHit::GetT);

        map_neuland_positions.clear(); 
        for(size_t j = 0; j < neuland_times.size(); ++j)
            map_neuland_positions[neuland_times[j]] = neuland_positions[j];

        auto neuland_first_time = map_neuland_positions.begin()->first;
        auto neuland_first_position = map_neuland_positions.begin()->second;

        h2_neuland_Z[0]->Fill(neuland_first_position.Z(), neuland_first_time);

        for(size_t j = 0; j < neuland_positions.size(); ++j)
        {
            h2_neuland_Z[1]->Fill(neuland_positions[j].Z(), neuland_positions[j].Y());
            h2_neuland_Z[2]->Fill(neuland_positions[j].Z(), neuland_positions[j].X());
            h2_neuland[1]->Fill(neuland_positions[j].X(), neuland_positions[j].Y());
            h1_neuland_Z[0]->Fill(neuland_positions[j].Z());
            h1_neuland[0]->Fill(neuland_times[j]);
            h1_neuland[1]->Fill(neuland_positions[j].X());
            h1_neuland[2]->Fill(neuland_positions[j].Y());
        }
    }

    // Write histograms to file
    std::vector<TH1*> all_hists;
    all_hists.clear();
    /*appendVector(all_hists, h1_dwr_main_foot);
    appendVector(all_hists, h1_dwr_main_califa);
    appendVector(all_hists, h1_dwr_main_califa_2hits);
    appendVector(all_hists, h1_dwr_main_califa_25F_p2p_24O);    
    appendVector(all_hists, h1_dwr_main_alpide);
    appendVector(all_hists, h1_dwr_main_los);
    appendVector(all_hists, h1_dwr_main_s2);
    */
    appendVector(all_hists, h1_neuland);
    appendVector(all_hists, h2_neuland);
    appendVector(all_hists, h1_neuland_Z);
    appendVector(all_hists, h2_neuland_Z);
    appendVector(all_hists, h2_frs_pid);
    //appendVector(all_hists, h2_tofd_pid);
    appendVector(all_hists, h2_AoZ_vs_Z);
    appendVector(all_hists, h2_califa_theta_corr);
    appendVector(all_hists, h2_califa_theta_corr_25F_p2p_24O);
    appendVector(all_hists, h2_califa_theta_corr_25F_p2p_23O);
    appendVector(all_hists, h2_califa_theta_corr_25F_ppn_24F);
    appendVector(all_hists, h2_califa_E_vs_theta);
 
    appendVector(all_hists, h2_califa_Ecorr);
    appendVector(all_hists, h2_califa_Ecorr_25F_p2p_24O);
    appendVector(all_hists, h2_califa_Ecorr_25F_p2p_23O);
    appendVector(all_hists, h2_califa_Ecorr_25F_ppn_24F);
    appendVector(all_hists, h2_vertex_XY);
    appendVector(all_hists, h2_vertex_YZ);
    appendVector(all_hists, h2_vertex_XZ);
    //appendVector(all_hists, h2_charge_foot_vs_tofd);
    appendVector(all_hists, h1_califa_opang);
    appendVector(all_hists, h1_califa_deltaphi);
    appendVector(all_hists, h1_califa_deltaphifolded);
    appendVector(all_hists, h1_califa_opang_25F_p2p_24O);
    appendVector(all_hists, h1_califa_opang_25F_p2p_23O);
    appendVector(all_hists, h1_califa_opang_25F_p2p_22O);
    appendVector(all_hists, h1_califa_opang_25F_ppn_24F);
    
    appendVector(all_hists, h2_califa_gamma_all);

    appendVector(all_hists, h2_califa_gamma_25F_ppn_24F);
    appendVector(all_hists, h2_califa_gamma_25F_p2p_24O);
    appendVector(all_hists, h2_califa_gamma_25F_p2p_23O);
    appendVector(all_hists, h1_califa_gamma_all);
    appendVector(all_hists, h1_califa_gamma_25F_ppn_24F);
    appendVector(all_hists, h1_califa_gamma_25F_p2p_24O);
    appendVector(all_hists, h1_califa_gamma_25F_p2p_23O);
    appendVector(all_hists, h2_califa_gamma_corr_all);
    appendVector(all_hists, h2_califa_gamma_corr_25F_ppn_24F);
    appendVector(all_hists, h2_califa_gamma_corr_25F_p2p_24O);
    appendVector(all_hists, h2_califa_gamma_corr_25F_p2p_23O);
    appendVector(all_hists, h2_califa_gamma_corr_25F_p2p_22O);
    appendVector(all_hists, h2_califa_gamma_corr_25F_p2p_22O_mleq3); 
    appendVector(all_hists, h1_califa_gamma_corr_all);
    appendVector(all_hists, h1_gamma_521_gated);
    appendVector(all_hists, h1_gamma_3118_gated);
    appendVector(all_hists, h1_gamma_Esum_mult2);

    appendVector(all_hists, h1_califa_gamma_corr_25F_ppn_24F_allmult);
    appendVector(all_hists, h1_califa_gamma_corr_25F_ppn_24F_mult1);
    appendVector(all_hists, h1_califa_gamma_corr_25F_ppn_24F_mult2);
    appendVector(all_hists, h1_califa_gamma_corr_25F_ppn_24F_mult3);
    appendVector(all_hists, h1_califa_gamma_corr_25F_ppn_24F_mult4);
    appendVector(all_hists, h1_califa_gamma_corr_25F_ppn_24F_mult5);
    
    appendVector(all_hists, h1_califa_gamma_corr_25F_p2p_24O);
    appendVector(all_hists, h1_califa_gamma_corr_25F_p2p_23O);
    appendVector(all_hists, h1_califa_gamma_corr_25F_p2p_22O);
    appendVector(all_hists, h1_califa_gamma_corr_25F_p2p_22O_allmult);
    // also save 24O multiplicity-specific histograms used for subtraction
    appendVector(all_hists, h1_califa_gamma_corr_25F_p2p_24O_mult_exact);
    appendVector(all_hists, h1_califa_gamma_corr_25F_p2p_24O_mult_leq);
    appendVector(all_hists, h1_califa_gamma_mult_25F_p2p_22O);
    appendVector(all_hists, h1_califa_gamma_corr_25F_p2p_22O_mult_exact);
    appendVector(all_hists, h1_califa_gamma_corr_25F_p2p_22O_mult_leq);
    appendVector(all_hists, h1_califa_gamma_corr_25F_p2p_22O_mult_leq_pairs_only);
    appendVector(all_hists, h2_proton_cluster_energy_corr_25F_p2p_22O);
    appendVector(all_hists, h2_gamma_theta_vs_proton1_theta_25F_p2p_22O);
    appendVector(all_hists, h2_gamma_theta_vs_proton2_theta_25F_p2p_22O);
    appendVector(all_hists, h1_gamma_gate_3199_corr_remaining_corr_25F_p2p_22O);
    appendVector(all_hists, h1_gamma_gate_3199_corr_remaining_raw_25F_p2p_22O);
    appendVector(all_hists, h1_gamma_gate_3199_corr_remaining_corr_25F_p2p_22O_mult_leq);
    appendVector(all_hists, h1_gamma_gate_3199_corr_remaining_raw_25F_p2p_22O_mult_leq);
    appendVector(all_hists, h1_gamma_gate_1383_corr_remaining_corr_25F_p2p_22O);
    appendVector(all_hists, h1_gamma_gate_1383_corr_remaining_raw_25F_p2p_22O);
    appendVector(all_hists, h1_gamma_gate_1383_corr_remaining_corr_25F_p2p_22O_mult_leq);
    appendVector(all_hists, h1_gamma_gate_1383_corr_remaining_raw_25F_p2p_22O_mult_leq);
    appendVector(all_hists, h1_frag_beta);
    appendVector(all_hists, h2_gamma_Epair_Meq2_25F_p2p_22O);
    appendVector(all_hists, h2_gamma_Tpair_Meq2_25F_p2p_22O);
    appendVector(all_hists, h2_gamma_thetapair_Meq2_25F_p2p_22O);
    appendVector(all_hists, h2_gamma_phipair_Meq2_25F_p2p_22O);
    appendVector(all_hists, h2_gamma_Epair_Mleq2_25F_p2p_22O);
    appendVector(all_hists, h2_gamma_Tpair_Mleq2_25F_p2p_22O);
    appendVector(all_hists, h2_gamma_Epair_Mleq3_25F_p2p_22O);
    appendVector(all_hists, h2_gamma_Tpair_Mleq3_25F_p2p_22O);
    appendVector(all_hists, h2_gamma_thetapair_Mleq3_25F_p2p_22O);
    appendVector(all_hists, h2_gamma_phipair_Mleq3_25F_p2p_22O);
    appendVector(all_hists, h2_gamma_Epair_Mleq4_25F_p2p_22O);
    appendVector(all_hists, h2_gamma_Tpair_Mleq4_25F_p2p_22O);
    appendVector(all_hists, h2_gamma_thetapair_Mleq4_25F_p2p_22O);
    appendVector(all_hists, h2_gamma_phipair_Mleq4_25F_p2p_22O);
    appendVector(all_hists, h2_gamma_Epair_Mleq5_25F_p2p_22O);
    appendVector(all_hists, h2_gamma_Tpair_Mleq5_25F_p2p_22O);
    appendVector(all_hists, h2_gamma_thetapair_Mleq5_25F_p2p_22O);
    appendVector(all_hists, h2_gamma_phipair_Mleq5_25F_p2p_22O);
    appendVector(all_hists, h1_DCA);
    appendVector(all_hists, h1_vertex_Z);
    appendVector(all_hists, h2_frag_beam_divergence);
    
    for(auto* h : all_hists)
        if(h) h->Write();
    
    outputFile->Close();
}

int main(int argc, char** argv)
{
    if (argc != 2) {
        std::cerr << "Usage: analyse_all_WR <inputrootfile>" << std::endl;
        return 1;
    }
    
    TString outrootfilename = argv[1];
    TString filename = argv[1];

    gStyle->SetOptStat("ieoumr");

	outrootfilename.ReplaceAll(".root", "_analysed.root");  // Replace .lmd by .root
	TString outFile = "/lustre/r3b/mjena/G249/rootfiles/analyzed/ts_4mus/with_randomization/" + outrootfilename;
	TFile* fout = new TFile(outFile, "RECREATE");
	if (fout->IsZombie() || !fout->IsOpen()) {
   	 std::cerr << "ERROR: Cannot create output file: " << outFile << std::endl;
   	 return 1;
	}
	const TString inputdir = "/lustre/r3b/mjena/G249/rootfiles/ts_4mus/with_randomization/";

	//Redefine root file name to include path to .root directory
	TString inputfilename = inputdir + "/" + filename;
	inputfilename.ReplaceAll("//", "/");

	run_analysis(inputfilename,fout);
	//std::cout << "Histograms saved to " << outFile << "\n";
	return 0;
    
}

