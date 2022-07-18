
/*
 
 plot of sigma(gamma-p) for semiforward and mid-rapidity paper, Alice data,
 other experiments and models
 
 compile:
 
 g++ plotForPaper.cxx `root-config --libs --cflags` -lEG -lMinuit -o plotForPaper
 
 run:
 
 ./plotForPaper
 
 */

// c++ headers
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>

// root headers
#include <Rtypes.h>
#include <TMath.h>
#include <TH2D.h>
#include <TFile.h>
#include <TF2.h>
#include <TF1.h>
#include <TRandom.h>
#include <TGraphErrors.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TMinuit.h>
#include <TStyle.h>
#include <TGaxis.h>
#include <TVirtualFitter.h>
#include <TDatabasePDG.h>
#include <TLatex.h>
#include <TLegend.h>

using namespace std;

// local headers
#include "dataAndModels.h"

//_____________________________________________________________________________
int main(void) {
    
    //configure layout of the plot
    gStyle->SetPadTickY(1);
    gStyle->SetTickLength(0.02,"Y");
    
    //TCanvas* c1 = new TCanvas("c1","c1",1000,700);
    TCanvas* c1 = new TCanvas("c1","c1",800,700);
    
    TH1F* frame1 = gPad->DrawFrame(gxmin,gymin,gxmax,gymax);
    frame1->GetXaxis()->SetMoreLogLabels();
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.11);//0.025
    gPad->SetBottomMargin(0.13);
    gPad->SetLogx();
    //gPad->SetLogy();
    
    gStyle->SetOptStat("");
    gStyle->SetPalette(1);
    gStyle->SetLineWidth(2);      //axis line
    gStyle->SetFrameLineWidth(2); //frame line
    TGaxis::SetMaxDigits(3);
    
    frame1->SetTitle(";W_{#gammap} (GeV);#frac{#sigma(#gamma+p #rightarrow J/#psi+X)}{#sigma(#gamma+p #rightarrow J/#psi+p)} ");
    Float_t siz = 0.045;
    frame1->SetTitleSize(siz);       frame1->SetLabelSize(siz);
    frame1->SetTitleSize(siz, "Y");  frame1->SetLabelSize(siz, "Y");
    frame1->GetXaxis()->SetTitleOffset(1.4);
    frame1->GetYaxis()->SetTitleOffset(1.5);
    
    //upper horizontal axis with Bjorken-x
    Double_t mjpsi = TDatabasePDG::Instance()->GetParticle(443)->Mass();
    Double_t bxmin = TMath::Power((mjpsi/gxmax),2.);
    Double_t bxmax = TMath::Power((mjpsi/gxmin),2.);
    TF1 *fbx = new TF1("fbx","TMath::Power(([0]/x),2.)", bxmin, bxmax);
    fbx->SetParameter(0, mjpsi);
    
    TGaxis *axis = new TGaxis(gxmax, gymax, gxmin, gymax, "fbx", 510, "+G");
    axis->SetTextFont(42);
    axis->SetLabelFont(42);
    axis->SetTitleSize(siz); axis->SetLabelSize(siz);
    axis->SetLabelOffset(-0.035);
    axis->SetTitleOffset();
    axis->Draw("same");
    
    TLatex *bxtit = new TLatex();
    bxtit->SetTextFont(42);
    bxtit->SetTextSize(siz);
    bxtit->SetTextAlign(31);
    bxtit->DrawLatex(gxmax, gymax+(gymax-gymin)*0.1, "Bjorken-#it{x}");
    
    // Draw the logo
    //  0: Just "ALICE" (for final data), to be added only if ALICE does not appear otherwise (e.g. in the legend)
    //  >0: ALICE Preliminary
    //DrawLogo(1, 0.58, 0.81);
    
    //end of layout configuration
    
    //alice fit
    alice *alic = alice::instance();
    /*
     TGraph *gShade = alic->fit();
     gShade->Draw("fsame");
     */
    
    //theoretical models
    
    //JIMWLK evolution
    TGraphErrors* gJIMWLK  = read_jimwlk("cgc_incoh_coh_ratio_evolution.txt");
    //SetGraphStyle(gJIMWLK  ,     1,kDot             ,  1,kRed   , 7,2);
    //SetGraphStyle(gJIMWLK, 1, kLine, 1, kBlack, 1, 3);
    //gJIMWLK->Draw("CF same");
    gJIMWLK->Draw("E3 same");
    TGraphErrors* gJIMWLK2  = read_jimwlk("cgc_incoh_coh_ratio_evolution.txt");
    SetGraphStyle(gJIMWLK2, 1, kLine, 1, kBlack, 1, 2);
    gJIMWLK2->Draw("lxsame");
    
    /*
    TGraphErrors* gJIMWLK3  = read_jimwlk2();
    SetGraphStyle(gJIMWLK3, 1, kLine, 1, kRed, 1, 2);
    gJIMWLK3->Draw("lsame");
     */
    
    // CCT
    TGraphErrors* gCCT = read_cct();
    SetGraphStyle(gCCT, kBlue, kOpenCircle, 0, kBlue, 1, 2);
    gCCT->Draw("lsame");
    
    
    //experiments
    
    //H1
    TGraphErrors* gH1_ee = read_h1();
    SetGraphStyle(gH1_ee, kBlue, kFullCircle, 1.4, kBlue, 1, 1);
    //gH1_ee->Draw("pzsame");
    gH1_ee->Draw("psame");
    
    //alice data
    TGraphErrors *alice8TDissStat       = alic->getForwardEightTeVDiss(true, false);    // first is stat, second is syst
    SetGraphStyle(alice8TDissStat, kMagenta, kFullDiamond, 0., kMagenta, 1, 4);
    alice8TDissStat->Draw("pzsame");
    TGraphErrors *alice8TDissTotal       = alic->getForwardEightTeVDiss(true, true);   // first is stat, second is syst
    SetGraphStyle(alice8TDissTotal, kMagenta, kFullDiamond, 2., kMagenta, 1, 2);
    alice8TDissTotal->Draw("psame");

    //legend
    Double_t xl = 0.55, dxl = 0.2;
    Double_t yl = 0.51, dyl = 0.2;
    TLegend* leg1 = new TLegend(xl, yl, xl+dxl, yl+dyl);
    leg1->SetFillStyle(0);
    leg1->SetBorderSize(0);
    leg1->SetTextSize(siz-0.005);
    leg1->AddEntry(alice8TDissTotal,"ALICE (#sqrt{#it{s}_{NN}} = 8.16 TeV)","p");
    leg1->AddEntry(gH1_ee,"H1","p");
    //leg1->AddEntry(gJIMWLK,"JIMWLK","l");
    leg1->AddEntry(gJIMWLK,"MS","l");
    leg1->AddEntry(gCCT,"CCT","l");
    leg1->Draw("same");
 
    
    //put vertical scale atop the models
    gPad->RedrawAxis("Y");
    
    c1->SaveAs("fig-ratio-final.eps");
    
    return 0;
    
}//main
