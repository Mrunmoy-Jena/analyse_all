// CALIFA gamma-response simulation using R3BCALIFATestGenerator
// with INCL-style vertex smearing conventions:
//   SetXYZ(x,y,z)       -> central vertex [cm]
//   SetDxDyDz(sx,sy,sz) -> x,y Gaussian sigmas [cm], z uniform in [z-sz, z+sz]
//
// Requires the updated R3BCALIFATestGenerator.{h,cxx} that adds SetDxDyDz().
//
// Example usage:
//   .L califaGammaResponseSim.C
//   califaGammaResponseSimSingle(1383., 200000, 0.8, 0.0, 0.1, 0.1, 2.5, 0.0, kTRUE);
//   .q
//   root -l -q 'califaGammaResponseSim.C'
//   buildGammaResponseHist("./sim_gamma_1383keV.root", "hResp_1383keV", 1383., 0.8, "./gamma_response_functions.root", kTRUE);

#include <TClonesArray.h>
#include <TFile.h>
#include <TH1D.h>
#include <TMath.h>
#include <TStopwatch.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>
#include <TRandom3.h>
#include <TVector3.h>
#include <TROOT.h>

namespace
{
    void ConfigureEnvironment()
    {
        TString dir = gSystem->Getenv("VMCWORKDIR");

        TString r3b_geomdir = dir + "/geometry/";
        gSystem->Setenv("GEOMPATH", r3b_geomdir.Data());
        r3b_geomdir.ReplaceAll("//", "/");

        TString r3b_confdir = dir + "/gconfig/";
        gSystem->Setenv("CONFIG_DIR", r3b_confdir.Data());
        r3b_confdir.ReplaceAll("//", "/");
    }

    Double_t DopplerCorrectEnergy(const Double_t eLabMeV, const Double_t beta, const Double_t thetaLabRad)
    {
        const Double_t gamma = 1.0 / TMath::Sqrt(1.0 - beta * beta);
        return gamma * (1.0 - beta * TMath::Cos(thetaLabRad)) * eLabMeV;
    }
}

void califaGammaResponseSimSingle(Double_t eGamma_keV = 1383.,
                                  Int_t nEvents = 200000,
                                  Double_t beta = 0.8,
                                  Double_t betaDispersion = 0.0,
                                  Double_t sigmaX_cm = 0.1,
                                  Double_t sigmaY_cm = 0.1,
                                  Double_t halfThicknessZ_cm = 2.5,
                                  Double_t zCenter_cm = 0.0,
                                  Bool_t includeTargetArea = kTRUE,
                                  Int_t randomSeed = 0,
                                  TString outFile = "",
                                  TString parFile = "")
{
    ConfigureEnvironment();
    gROOT->SetBatch(kTRUE);

    if (outFile.IsNull())
        outFile = Form("./sim_gamma_%dkeV.root", (int)TMath::Nint(eGamma_keV));
    if (parFile.IsNull())
        parFile = Form("./par_gamma_%dkeV.root", (int)TMath::Nint(eGamma_keV));

    if (beta <= 0.0 || beta >= 1.0)
    {
        Error("califaGammaResponseSimSingle", "beta must satisfy 0 < beta < 1. Current value = %.6f", beta);
        return;
    }

    TStopwatch timer;
    timer.Start();

    FairRunSim* run = new FairRunSim();
    run->SetName("TGeant4");
    run->SetOutputFile(outFile.Data());
    FairRuntimeDb* rtdb = run->GetRuntimeDb();
    
    FairParAsciiFileIo *parIo1 = new FairParAsciiFileIo();
    parIo1->open("califaDigi.par", "in");  //important
    rtdb->setFirstInput(parIo1);
    rtdb->print();

    rtdb->getContainer("califaCrystalPars4Sim");

    UInt_t runId = 1;
    rtdb->initContainers(runId);
    
    run->SetMaterials("media_r3b.geo");

    FairModule* cave = new R3BCave("CAVE");
    cave->SetGeometryFileName("r3b_cave.geo");
    run->AddModule(cave);

    if (includeTargetArea)
    {
        auto* foots = new R3BTra("target_area_v2025_5cm.geo.root");
        run->AddModule(foots);
    }

    // Keep the same CALIFA placement logic as your current gamma macro.
    auto* calsim = new R3BCalifa("califa_v2025.6.geo.root");
    calsim->SelectGeometryVersion(2025);
    run->AddModule(calsim);

    // --------------------------------------------------------------------
    // Generator: use updated R3BCALIFATestGenerator with INCL-style smearing
    // --------------------------------------------------------------------
    auto* gammaGen = new R3BCALIFATestGenerator(22, 1); // PDG 22 = gamma, multiplicity 1

    const Double_t eGammaGeV = eGamma_keV * 1.e-6; // keV -> GeV
    gammaGen->SetPRange(eGammaGeV, eGammaGeV);     // for gamma: p = E in GeV
    gammaGen->SetThetaRange(0., 180.);
    gammaGen->SetCosTheta();
    gammaGen->SetPhiRange(0., 360.);

    // INCL-style vertex convention:
    // x,y Gaussian around central point; z uniform in [zCenter-sz, zCenter+sz]
    gammaGen->SetXYZ(0., 0., zCenter_cm);
    gammaGen->SetDxDyDz(sigmaX_cm, sigmaY_cm, halfThicknessZ_cm);

    // Lorentz boost of emitting fragment.
    gammaGen->SetFragmentVelocity(beta, betaDispersion);
    gammaGen->SetLorentzBoost(beta);

    auto* primGen = new FairPrimaryGenerator();
    primGen->AddGenerator(gammaGen);
    run->SetGenerator(primGen);
    run->SetStoreTraj(kTRUE);

    FairLogger::GetLogger()->SetLogVerbosityLevel("LOW");
    FairLogger::GetLogger()->SetLogScreenLevel("INFO");

    auto* califaDig = new R3BCalifaDigitizer();
    califaDig->SetNonUniformity(1.);
    califaDig->SetRealConfig(true);  //should apply effects of res. smearing crystal by crystal (better than a global, fixed smearing)
    //califaDig->SetExpGammaEnergyRes(6.);
    //califaDig->SetComponentRes(6.);
    //califaDig->SetDetectionThreshold(0.000010); // GeV = 10 keV
    run->AddTask(califaDig);

    auto* cal2Clus = new R3BCalifaCrystalCal2Cluster();
    cal2Clus->SetCrystalThreshold(0.1);       // 100 keV
    cal2Clus->SetProtonClusterThreshold(18.); // irrelevant for gammas but keep same logic
    cal2Clus->SetGammaClusterThreshold(0.1);  // 100 keV
    cal2Clus->SelectGeometryVersion(2025);
    cal2Clus->SetRoundWindow(0.4);
    cal2Clus->SetRandomization(kTRUE);
    TString randomization_file = "/home/e12exp/ssd/R3BParams_g249/califa/gammas_minus_4.root";
    randomization_file.ReplaceAll("//", "/");
    cal2Clus->SetRandomizationFile(randomization_file);
    run->AddTask(cal2Clus);

    run->Init();
    TVirtualMC::GetMC()->SetRandom(new TRandom3(randomSeed));
    TVirtualMC::GetMC()->SetMaxNStep(-15000);

    Bool_t kParameterMerged = kTRUE;
    auto* parOut = new FairParRootFileIo(kParameterMerged);
    parOut->open(parFile.Data());
    rtdb->setOutput(parOut);
    rtdb->saveOutput();

    if (nEvents > 0)
        run->Run(nEvents);

    timer.Stop();
    std::cout << "\nFinished single-line gamma simulation.\n"
              << "  E_gamma(rest)          = " << eGamma_keV << " keV\n"
              << "  beta(mean)             = " << beta << "\n"
              << "  beta(dispersion)       = " << betaDispersion << "\n"
              << "  sigmaX                 = " << sigmaX_cm << " cm\n"
              << "  sigmaY                 = " << sigmaY_cm << " cm\n"
              << "  z uniform half-width   = " << halfThicknessZ_cm << " cm\n"
              << "  z center               = " << zCenter_cm << " cm\n"
              << "  include target area    = " << (includeTargetArea ? "yes" : "no") << "\n"
              << "  output                 = " << outFile << "\n"
              << "  parameter              = " << parFile << "\n"
              << "  real time              = " << timer.RealTime() << " s\n"
              << "  cpu time               = " << timer.CpuTime() << " s\n\n";
}

void buildGammaResponseHist(TString simFile,
                            TString histName,
                            Double_t eGamma_keV,
			    Double_t betaForCorrection,
                            TString outRoot = "./gamma_response_functions.root",
                            Bool_t recreate = kFALSE,
                            TString treeName = "evt",
                            TString branchName = "CalifaClusterData")
{
    auto* inFile = TFile::Open(simFile, "READ");
    if (!inFile || inFile->IsZombie())
    {
        Error("buildGammaResponseHist", "Cannot open %s", simFile.Data());
        return;
    }

    TTree* tree = dynamic_cast<TTree*>(inFile->Get(treeName));
    if (!tree)
    {
        Error("buildGammaResponseHist", "Tree '%s' not found in %s", treeName.Data(), simFile.Data());
        inFile->Close();
        return;
    }

    auto* clusters = new TClonesArray("R3BCalifaClusterData");
    tree->SetBranchAddress(branchName.Data(), &clusters);

    const Double_t eGammaMeV = eGamma_keV * 1.e-3;
    auto* hResp = new TH1D(histName.Data(),
                           Form("CALIFA response, E_{0}=%.0f keV", eGamma_keV),
                           1600,
                           0,
                           eGammaMeV * 1.8);
    hResp->GetXaxis()->SetTitle("E_{#gamma}^{Doppler corr.} [MeV]");
    hResp->GetYaxis()->SetTitle("Counts");
    hResp->Sumw2();

    const Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i)
    {
        tree->GetEntry(i);
        const Int_t nClusters = clusters->GetEntriesFast();
        if (nClusters <= 0)
            continue;

        Double_t eRestAddback = 0.;
        for (Int_t j = 0; j < nClusters; ++j)
        {
            auto* cl = dynamic_cast<R3BCalifaClusterData*>(clusters->At(j));
            if (!cl)
                continue;

            const Double_t eLab = cl->GetEnergy();
            const Double_t theta = cl->GetTheta();
            const Double_t eRest = DopplerCorrectEnergy(eLab, betaForCorrection, theta);

            if (eRest > 0.)
                eRestAddback += eRest;
        }

        if (eRestAddback > 0.)
            hResp->Fill(eRestAddback);
    }

    auto* outFile = TFile::Open(outRoot, recreate ? "RECREATE" : "UPDATE");
    hResp->Write(histName.Data(), TObject::kOverwrite);
    outFile->Close();
    inFile->Close();

    std::cout << "Wrote histogram '" << histName << "' to " << outRoot << "\n";
}
