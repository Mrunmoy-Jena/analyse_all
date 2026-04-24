#include <TFile.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TString.h>
#include <iostream>
#include <TF1.h>
#include <TCanvas.h>
#include <TPaveText.h>
#include <TMath.h>
#include <TAxis.h>
#include <TROOT.h>

// ------------------------------------------------------------
// Global response histograms
// ------------------------------------------------------------
static TH1D* gResp1 = nullptr; // 1218 keV response
static TH1D* gResp2 = nullptr; // 1383 keV response
static TH1D* gResp3 = nullptr; // 1710 keV response
static TH1D* gResp4 = nullptr; // 2354 keV response
static TH1D* gResp5 = nullptr; // 2600 keV response
static TH1D* gResp6 = nullptr; // 3199 keV response
static TH1D* gResp7 = nullptr; // 3310 keV response
static TH1D* gResp8 = nullptr; // 3710 keV response

// ------------------------------------------------------------
// Fixed background shape parameters from 24O fit
// Background model:
//   A_bkg * ( exp(p0 + p1*x) + exp(p2 + p3*x) )
// ------------------------------------------------------------
static const Double_t kBgP0 = 6.675;
static const Double_t kBgP1 = -2.45;
static const Double_t kBgP2 = 5.204;
static const Double_t kBgP3 = -0.4046;

// ------------------------------------------------------------
// Evaluate histogram template with LINEAR INTERPOLATION between
// bin centres.
// ------------------------------------------------------------
Double_t EvalTemplate(TH1D* h, Double_t x)
{
    if (!h)
        return 0.0;

    int bin = h->FindBin(x);

    if (bin <= 0 || bin > h->GetNbinsX())
        return 0.0;

    Double_t xc = h->GetXaxis()->GetBinCenter(bin);

    int binN = (x >= xc) ? bin + 1 : bin - 1;

    if (binN <= 0 || binN > h->GetNbinsX())
        return h->GetBinContent(bin);

    Double_t xN = h->GetXaxis()->GetBinCenter(binN);
    Double_t y0 = h->GetBinContent(bin);
    Double_t y1 = h->GetBinContent(binN);

    Double_t t = (x - xc) / (xN - xc);
    return y0 + t * (y1 - y0);
}

// ------------------------------------------------------------
// Make a TF1 from a histogram template with a fixed scale,
// for drawing only.
// ------------------------------------------------------------
TF1* MakeTemplateTF1(TH1D* h,
                     const TString& name,
                     Double_t scale)
{
    if (!h)
        return nullptr;

    auto f = new TF1(
        name,
        [h, scale](double* xx, double*) -> double
        {
            return scale * EvalTemplate(h, xx[0]);
        },
        h->GetXaxis()->GetXmin(),
        h->GetXaxis()->GetXmax(),
        0);

    f->SetNpx(2000);
    return f;
}

// ------------------------------------------------------------
// Fixed-shape background evaluator
// ------------------------------------------------------------
Double_t EvalFixedBackgroundShape(Double_t x)
{
    return TMath::Exp(kBgP0 + kBgP1 * x) + TMath::Exp(kBgP2 + kBgP3 * x);
}

// ------------------------------------------------------------
// Model: response templates (no shift) + fixed-shape background
//
// par[0] = A1       scale of response 1
// par[1] = A2       scale of response 2
// par[2] = A3       scale of response 3
// par[3] = A4       scale of response 4
// par[4] = A5       scale of response 5
// par[5] = A6       scale of response 6
// par[6] = A7       scale of response 7
// par[7] = A8       scale of response 8
// par[8] = A_bkg    overall scale of fixed 24O background shape
// ------------------------------------------------------------
Double_t fitModelResp(Double_t* x, Double_t* par)
{
    const Double_t xx = x[0];

    const Double_t r1  = par[0] * EvalTemplate(gResp1, xx);
    const Double_t r2  = par[1] * EvalTemplate(gResp2, xx);
    const Double_t r3  = par[2] * EvalTemplate(gResp3, xx);
    const Double_t r4  = par[3] * EvalTemplate(gResp4, xx);
    const Double_t r5  = par[4] * EvalTemplate(gResp5, xx);
    const Double_t r6  = par[5] * EvalTemplate(gResp6, xx);
    const Double_t r7  = par[6] * EvalTemplate(gResp7, xx);
    const Double_t r8  = par[7] * EvalTemplate(gResp8, xx);
    const Double_t bkg = par[8] * EvalFixedBackgroundShape(xx);

    return r1 + r2 + r3 + r4 + r5 + r6 + r7 + r8 + bkg;
}

// ------------------------------------------------------------
// Initialize fit parameters on a fresh TF1 (9 parameters)
// ------------------------------------------------------------
void InitFitParameters(TF1& f, TH1F* h)
{
    const Double_t maxv = h->GetMaximum();

    f.SetParName(0, "A1");
    f.SetParName(1, "A2");
    f.SetParName(2, "A3");
    f.SetParName(3, "A4");
    f.SetParName(4, "A5");
    f.SetParName(5, "A6");
    f.SetParName(6, "A7");
    f.SetParName(7, "A8");
    f.SetParName(8, "A_bkg");

    f.SetParameter(0, maxv * 0.10);
    f.SetParameter(1, maxv * 0.50);
    f.SetParameter(2, maxv * 0.10);
    f.SetParameter(3, maxv * 0.10);
    f.SetParameter(4, maxv * 0.10);
    f.SetParameter(5, maxv * 0.20);
    f.SetParameter(6, maxv * 0.10);
    f.SetParameter(7, maxv * 0.10);
    f.SetParameter(8, 1.0); // global scale on fixed 24O background shape

    //f.SetParLimits(0, 0.0, 1.0e9);
    //f.SetParLimits(1, 0.0, 1.0e9);
    //f.SetParLimits(2, 0.0, 1.0e9);
    //f.SetParLimits(3, 0.0, 1.0e9);
    //f.SetParLimits(4, 0.0, 1.0e9);
    //f.SetParLimits(5, 0.0, 1.0e9);
    //f.SetParLimits(6, 0.0, 1.0e9);
    //f.SetParLimits(7, 0.0, 1.0e9);
    //f.SetParLimits(8, 0.0, 1.0e3);
}

// ------------------------------------------------------------
// Draw and save one fitted histogram
// ------------------------------------------------------------
void DrawAndWriteFit(TFile* fout,
                     TH1F* h,
                     TF1& fitFunc,
                     const TString& canvasName,
                     const TString& fitName,
                     const TString& labelText,
                     const Double_t xminFit,
                     const Double_t xmaxFit)
{
    if (!fout || !h)
        return;

    fout->cd();

    const Double_t A1    = fitFunc.GetParameter(0);
    const Double_t A2    = fitFunc.GetParameter(1);
    const Double_t A3    = fitFunc.GetParameter(2);
    const Double_t A4    = fitFunc.GetParameter(3);
    const Double_t A5    = fitFunc.GetParameter(4);
    const Double_t A6    = fitFunc.GetParameter(5);
    const Double_t A7    = fitFunc.GetParameter(6);
    const Double_t A8    = fitFunc.GetParameter(7);
    const Double_t A_bkg = fitFunc.GetParameter(8);

    const Double_t sigma1 = gResp1 ? gResp1->GetRMS() : 0.0;
    const Double_t sigma2 = gResp2 ? gResp2->GetRMS() : 0.0;
    const Double_t sigma3 = gResp3 ? gResp3->GetRMS() : 0.0;
    const Double_t sigma4 = gResp4 ? gResp4->GetRMS() : 0.0;
    const Double_t sigma5 = gResp5 ? gResp5->GetRMS() : 0.0;
    const Double_t sigma6 = gResp6 ? gResp6->GetRMS() : 0.0;
    const Double_t sigma7 = gResp7 ? gResp7->GetRMS() : 0.0;
    const Double_t sigma8 = gResp8 ? gResp8->GetRMS() : 0.0;
    const Double_t peak1  = gResp1 ? gResp1->GetXaxis()->GetBinCenter(gResp1->GetMaximumBin()) : 0.0;
    const Double_t peak2  = gResp2 ? gResp2->GetXaxis()->GetBinCenter(gResp2->GetMaximumBin()) : 0.0;
    const Double_t peak3  = gResp3 ? gResp3->GetXaxis()->GetBinCenter(gResp3->GetMaximumBin()) : 0.0;
    const Double_t peak4  = gResp4 ? gResp4->GetXaxis()->GetBinCenter(gResp4->GetMaximumBin()) : 0.0;
    const Double_t peak5  = gResp5 ? gResp5->GetXaxis()->GetBinCenter(gResp5->GetMaximumBin()) : 0.0;
    const Double_t peak6  = gResp6 ? gResp6->GetXaxis()->GetBinCenter(gResp6->GetMaximumBin()) : 0.0;
    const Double_t peak7  = gResp7 ? gResp7->GetXaxis()->GetBinCenter(gResp7->GetMaximumBin()) : 0.0;
    const Double_t peak8  = gResp8 ? gResp8->GetXaxis()->GetBinCenter(gResp8->GetMaximumBin()) : 0.0;

    TCanvas* c = new TCanvas(canvasName, canvasName, 900, 650);

    h->SetDirectory(nullptr);
    h->GetListOfFunctions()->Clear();
    h->SetMarkerStyle(24);
    h->SetMarkerSize(0.65);
    h->SetLineColor(kBlack);
    h->Draw("E");

    // Background: fixed-shape, only scaled by A_bkg
    TF1 f_bkg(Form("%s_bkg", canvasName.Data()),
              [A_bkg](double* xx, double*) -> double
              {
                  return A_bkg * EvalFixedBackgroundShape(xx[0]);
              },
              xminFit, xmaxFit, 0);
    f_bkg.SetLineColor(kBlue);
    f_bkg.SetLineStyle(2);
    f_bkg.SetLineWidth(1);
    f_bkg.SetNpx(2000);
    f_bkg.Draw("same");

    TF1* f_r1 = MakeTemplateTF1(gResp1, Form("%s_r1", canvasName.Data()), A1);
    if (f_r1) { f_r1->SetLineColor(kRed);        f_r1->SetLineWidth(1); f_r1->Draw("same"); }

    TF1* f_r2 = MakeTemplateTF1(gResp2, Form("%s_r2", canvasName.Data()), A2);
    if (f_r2) { f_r2->SetLineColor(kGreen + 2);  f_r2->SetLineWidth(1); f_r2->Draw("same"); }

    TF1* f_r3 = MakeTemplateTF1(gResp3, Form("%s_r3", canvasName.Data()), A3);
    if (f_r3) { f_r3->SetLineColor(kMagenta + 1);f_r3->SetLineWidth(1); f_r3->Draw("same"); }

    TF1* f_r4 = MakeTemplateTF1(gResp4, Form("%s_r4", canvasName.Data()), A4);
    if (f_r4) { f_r4->SetLineColor(kOrange + 7); f_r4->SetLineWidth(1); f_r4->Draw("same"); }

    TF1* f_r5 = MakeTemplateTF1(gResp5, Form("%s_r5", canvasName.Data()), A5);
    if (f_r5) { f_r5->SetLineColor(kCyan + 2);   f_r5->SetLineWidth(1); f_r5->Draw("same"); }

    TF1* f_r6 = MakeTemplateTF1(gResp6, Form("%s_r6", canvasName.Data()), A6);
    if (f_r6) { f_r6->SetLineColor(kViolet + 1); f_r6->SetLineWidth(1); f_r6->Draw("same"); }

    TF1* f_r7 = MakeTemplateTF1(gResp7, Form("%s_r7", canvasName.Data()), A7);
    if (f_r7) { f_r7->SetLineColor(kSpring + 5); f_r7->SetLineWidth(1); f_r7->Draw("same"); }

    TF1* f_r8 = MakeTemplateTF1(gResp8, Form("%s_r8", canvasName.Data()), A8);
    if (f_r8) { f_r8->SetLineColor(kPink + 7);   f_r8->SetLineWidth(1); f_r8->Draw("same"); }

    fitFunc.SetNpx(2000);
    fitFunc.SetLineColor(kBlack);
    fitFunc.SetLineWidth(1);
    fitFunc.Draw("same");

    TPaveText* pt = new TPaveText(0.50, 0.35, 0.95, 0.90, "NDC");
    pt->SetFillColor(0);
    pt->SetBorderSize(1);
    pt->AddText(labelText);
    pt->AddText(Form("#chi^{2}/ndf = %.2f",
                     fitFunc.GetNDF() > 0 ? fitFunc.GetChisquare() / fitFunc.GetNDF() : -1.0));
    pt->AddText(Form("A1 = %.1f", A1));
    pt->AddText(Form("A2 = %.1f", A2));
    pt->AddText(Form("A3 = %.1f", A3));
    pt->AddText(Form("A4 = %.1f", A4));
    pt->AddText(Form("A5 = %.1f", A5));
    pt->AddText(Form("A6 = %.1f", A6));
    pt->AddText(Form("A7 = %.1f", A7));
    pt->AddText(Form("A8 = %.1f", A8));
    pt->AddText(Form("A_bkg = %.4f", A_bkg));
    pt->AddText(Form("Resp1 peak = %.4f, sigma = %.3f", peak1, sigma1));
    pt->AddText(Form("Resp2 peak = %.4f, sigma = %.3f", peak2, sigma2));
    pt->AddText(Form("Resp3 peak = %.4f, sigma = %.3f", peak3, sigma3));
    pt->AddText(Form("Resp4 peak = %.4f, sigma = %.3f", peak4, sigma4));
    pt->AddText(Form("Resp5 peak = %.4f, sigma = %.3f", peak5, sigma5));
    pt->AddText(Form("Resp6 peak = %.4f, sigma = %.3f", peak6, sigma6));
    pt->AddText(Form("Resp7 peak = %.4f, sigma = %.3f", peak7, sigma7));
    pt->AddText(Form("Resp8 peak = %.4f, sigma = %.3f", peak8, sigma8));
    pt->AddText(Form("Fixed bg: exp(%.3f + %.4f x) + exp(%.3f + %.4f x)",
                     kBgP0, kBgP1, kBgP2, kBgP3));
    pt->Draw("same");

    c->Write();

    delete f_r1;
    delete f_r2;
    delete f_r3;
    delete f_r4;
    delete f_r5;
    delete f_r6;
    delete f_r7;
    delete f_r8;

    fitFunc.Write(fitName);
    f_bkg.Write(Form("%s_bkg", fitName.Data()));
}

// ------------------------------------------------------------
// Main macro
// ------------------------------------------------------------
void fit_gamma_with_response(
    TString inputFile = "/lustre/r3b/mjena/G249/rootfiles/analyzed/ts_4mus/wo_randomization/combined_analysed.root",
    TString respFile  = "gamma_response_functions.root",
    Double_t xminFit  = 0.5,
    Double_t xmaxFit  = 10.0,
    Int_t nSmooth     = 1)
{
    gROOT->SetBatch(kTRUE);

    TFile* f = TFile::Open(inputFile, "UPDATE");
    if (!f || f->IsZombie())
    {
        std::cerr << "Cannot open input file: " << inputFile << std::endl;
        return;
    }

    TFile* fResp = TFile::Open(respFile, "READ");
    if (!fResp || fResp->IsZombie())
    {
        std::cerr << "Cannot open response file: " << respFile << std::endl;
        f->Close();
        return;
    }

    gResp1 = (TH1D*)fResp->Get("hResp_1218keV");
    gResp2 = (TH1D*)fResp->Get("hResp_1383keV");
    gResp3 = (TH1D*)fResp->Get("hResp_1710keV");
    gResp4 = (TH1D*)fResp->Get("hResp_2354keV");
    gResp5 = (TH1D*)fResp->Get("hResp_2600keV");
    gResp6 = (TH1D*)fResp->Get("hResp_3199keV");
    gResp7 = (TH1D*)fResp->Get("hResp_3310keV");
    gResp8 = (TH1D*)fResp->Get("hResp_3710keV");

    if (!gResp1 || !gResp2 || !gResp3 || !gResp4 || !gResp5 || !gResp6 || !gResp7 || !gResp8)
    {
        std::cerr << "Could not find one or more response histograms in " << respFile << std::endl;
        fResp->Close();
        f->Close();
        return;
    }

    gResp1 = (TH1D*)gResp1->Clone("gResp1_local");
    gResp2 = (TH1D*)gResp2->Clone("gResp2_local");
    gResp3 = (TH1D*)gResp3->Clone("gResp3_local");
    gResp4 = (TH1D*)gResp4->Clone("gResp4_local");
    gResp5 = (TH1D*)gResp5->Clone("gResp5_local");
    gResp6 = (TH1D*)gResp6->Clone("gResp6_local");
    gResp7 = (TH1D*)gResp7->Clone("gResp7_local");
    gResp8 = (TH1D*)gResp8->Clone("gResp8_local");
    gResp1->SetDirectory(nullptr);
    gResp2->SetDirectory(nullptr);
    gResp3->SetDirectory(nullptr);
    gResp4->SetDirectory(nullptr);
    gResp5->SetDirectory(nullptr);
    gResp6->SetDirectory(nullptr);
    gResp7->SetDirectory(nullptr);
    gResp8->SetDirectory(nullptr);
    fResp->Close();

    const Double_t int1 = gResp1->Integral();
    const Double_t int2 = gResp2->Integral();
    const Double_t int3 = gResp3->Integral();
    const Double_t int4 = gResp4->Integral();
    const Double_t int5 = gResp5->Integral();
    const Double_t int6 = gResp6->Integral();
    const Double_t int7 = gResp7->Integral();
    const Double_t int8 = gResp8->Integral();

    if (int1 > 0.0) gResp1->Scale(1.0 / int1); else std::cerr << "Warning: response 1 has zero integral." << std::endl;
    if (int2 > 0.0) gResp2->Scale(1.0 / int2); else std::cerr << "Warning: response 2 has zero integral." << std::endl;
    if (int3 > 0.0) gResp3->Scale(1.0 / int3); else std::cerr << "Warning: response 3 has zero integral." << std::endl;
    if (int4 > 0.0) gResp4->Scale(1.0 / int4); else std::cerr << "Warning: response 4 has zero integral." << std::endl;
    if (int5 > 0.0) gResp5->Scale(1.0 / int5); else std::cerr << "Warning: response 5 has zero integral." << std::endl;
    if (int6 > 0.0) gResp6->Scale(1.0 / int6); else std::cerr << "Warning: response 6 has zero integral." << std::endl;
    if (int7 > 0.0) gResp7->Scale(1.0 / int7); else std::cerr << "Warning: response 7 has zero integral." << std::endl;
    if (int8 > 0.0) gResp8->Scale(1.0 / int8); else std::cerr << "Warning: response 8 has zero integral." << std::endl;

    if (nSmooth > 0)
    {
        for (Int_t i = 0; i < nSmooth; ++i)
        {
            gResp1->Smooth(1);
            gResp2->Smooth(1);
            gResp3->Smooth(1);
            gResp4->Smooth(1);
            gResp5->Smooth(1);
            gResp6->Smooth(1);
            gResp7->Smooth(1);
            gResp8->Smooth(1);
        }
    }

    f->cd();

    TH1F* h22O_all = (TH1F*)f->Get("h1_califa_gamma_corr_25F_p2p_22O_allmult_0");
    if (!h22O_all)
    {
        std::cerr << "Cannot find histogram: h1_califa_gamma_corr_25F_p2p_22O_allmult_0" << std::endl;
    }
    else
    {
        TF1 f_all("f_all", fitModelResp, xminFit, xmaxFit, 9);
        InitFitParameters(f_all, h22O_all);

        h22O_all->Fit(&f_all, "RQN", "", xminFit, xmaxFit);

        DrawAndWriteFit(f,
                        h22O_all,
                        f_all,
                        "c_fit_22O_all",
                        "fit_22O_all",
                        "All multiplicities",
                        xminFit,
                        xmaxFit);
    }

    for (int m = 1; m <= 8; ++m)
    {
        const char* types[] = { "mult_exact", "mult_leq" };

        for (const char* t : types)
        {
            TString histName = Form("h1_califa_gamma_corr_25F_p2p_22O_%s_%d", t, m);
            TH1F* h = (TH1F*)f->Get(histName);

            if (!h)
            {
                std::cout << "Warning: Cannot find histogram: " << histName << std::endl;
                continue;
            }

            TF1 f_mult(Form("f_%s_%d", t, m), fitModelResp, xminFit, xmaxFit, 9);
            InitFitParameters(f_mult, h);

            h->Fit(&f_mult, "RQN", "", xminFit, xmaxFit);

            DrawAndWriteFit(f,
                            h,
                            f_mult,
                            Form("c_fit_22O_%s_%d", t, m),
                            Form("fit_22O_%s_%d", t, m),
                            Form("%s, m = %d", t, m),
                            xminFit,
                            xmaxFit);
        }
    }

    f->Close();

    delete gResp1;
    delete gResp2;
    delete gResp3;
    delete gResp4;
    delete gResp5;
    delete gResp6;
    delete gResp7;
    delete gResp8;
    gResp1 = nullptr;
    gResp2 = nullptr;
    gResp3 = nullptr;
    gResp4 = nullptr;
    gResp5 = nullptr;
    gResp6 = nullptr;
    gResp7 = nullptr;
    gResp8 = nullptr;

    std::cout << "Done. Fits appended to: " << inputFile << std::endl;
}
