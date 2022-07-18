
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
int main() {  // if 1: 2.5 < y < 3.25; if 2: 3.25 < y < 4.0;
    int rapBin = 1;
    if (rapBin != 1 && rapBin != 2) return -1;
    //cout << "rapBin = " << rapBin << endl;
    //configure layout of the plot
    gStyle->SetPadTickY(1);
    gStyle->SetTickLength(0.02,"Y");
    
    TCanvas* c1 = new TCanvas("c1","c1",1000,700);
    
    TH1F* frame1 = gPad->DrawFrame(gxmin,gymin,gxmax,gymax);
    frame1->GetXaxis()->SetMoreLogLabels();
    gPad->SetLeftMargin(0.11);
    gPad->SetRightMargin(0.03);
    gPad->SetTopMargin(0.11);//0.025
    gPad->SetBottomMargin(0.13);
    //gPad->SetLogx();
    gPad->SetLogy();
    
    gStyle->SetOptStat("");
    gStyle->SetPalette(1);
    gStyle->SetLineWidth(2);      //axis line
    gStyle->SetFrameLineWidth(2); //frame line
    TGaxis::SetMaxDigits(3);
    
    frame1->SetTitle(";M_{#mu#mu} (GeV/#it{c}^{2});#frac{d #sigma}{dM_{#mu#mu}}(#gamma#gamma #rightarrow #mu#mu) (#mub/(GeV/#it{c}^{2}))");
    Float_t siz = 0.045;
    frame1->SetTitleSize(siz);       frame1->SetLabelSize(siz);
    frame1->SetTitleSize(siz, "Y");  frame1->SetLabelSize(siz, "Y");
    frame1->GetXaxis()->SetTitleOffset(1.4);
    frame1->GetYaxis()->SetTitleOffset(1.0);
    
    /*
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
     */
    
    TLatex *bxtit = new TLatex();
    bxtit->SetTextFont(42);
    bxtit->SetTextSize(siz);
    bxtit->SetTextAlign(31);
    //bxtit->DrawLatex(gxmax, gymax+450, "Bjorken-#it{x}");
    
    TLatex * textRap = nullptr;
    if (rapBin == 1) textRap = new TLatex(0.69,0.78,"#bf{2.50 < y < 3.25}");
    else textRap = new TLatex(0.69,0.78,"#bf{3.25 < y < 4.00}");
    textRap->SetNDC();
    textRap->Draw();
    TLatex * textPt = new TLatex(0.69,0.70,"#bf{p_{T} < 3 GeV/#it{c}}");
    textPt->SetNDC();
    textPt->Draw();
    
    // Draw the logo
    //  0: Just "ALICE" (for final data), to be added only if ALICE does not appear otherwise (e.g. in the legend)
    //  >0: ALICE Preliminary
    //DrawLogo(1, 0.58, 0.81);
    
    //end of layout configuration
    
    //alice fit
    alice *alic = alice::instance();
    
    
    //alice data
    TGraphErrors *aliceFwd1      = alic->getForwardEightTeV1();
    TGraphErrors *aliceFwd2      = alic->getForwardEightTeV2();
    if (rapBin == 1) aliceFwd2->Draw("psame");
    else aliceFwd1->Draw("psame");
    
    // STARlight
    TGraphErrors *sl1      = alic->getForwardEightTeVSL1();
    TGraphErrors *sl2      = alic->getForwardEightTeVSL2();
    if (rapBin == 1) sl2->Draw("psame");
    else sl1->Draw("psame");

    
    //legend
    Double_t xl = 0.17, dxl = 0.4;
    Double_t yl = 0.18, dyl = 0.28;
    dyl = 0.17;
    TLegend* leg2 = new TLegend(xl, yl, xl+dxl, yl+dyl);
    leg2->SetFillStyle(0);
    leg2->SetBorderSize(0);
    //leg2->SetTextSize(siz-0.005);
    leg2->SetTextSize(siz*0.9);
    if (rapBin == 1) leg2->AddEntry(aliceFwd2,"ALICE (#sqrt{#it{s}_{NN}} = 8.16 TeV)","lp");
    else leg2->AddEntry(aliceFwd1,"ALICE (#sqrt{#it{s}_{NN}} = 8.16 TeV)","lp");
    if (rapBin == 1) leg2->AddEntry(sl2,"STARlight","lp");
    else leg2->AddEntry(sl1,"STARlight","lp");
    leg2->Draw("same");
    
    
    //put vertical scale atop the models
    gPad->RedrawAxis("Y");
    
    if (rapBin == 1) c1->SaveAs("fig-gg-mass-final-1.eps");
    else c1->SaveAs("fig-gg-mass-final-2.eps");
    
    return 0;
    
}//main
