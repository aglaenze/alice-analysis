
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
#include <TGraphAsymmErrors.h>
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
    
    bool drawPreviousMeas = true;
    bool drawModels = false;
    
    TCanvas* c1 = new TCanvas("c1","c1",1000,700);
    
    TH1F* frame1 = gPad->DrawFrame(gxmin,gymin,gxmax,gymax);
    frame1->GetXaxis()->SetMoreLogLabels();
    gPad->SetLeftMargin(0.09);
    gPad->SetRightMargin(0.08);
    gPad->SetTopMargin(0.11);//0.025
    gPad->SetBottomMargin(0.13);
    gPad->SetLogx();
    gPad->SetLogy();
    
    gStyle->SetOptStat("");
    gStyle->SetPalette(1);
    gStyle->SetLineWidth(2);      //axis line
    gStyle->SetFrameLineWidth(2); //frame line
    TGaxis::SetMaxDigits(3);
    
    frame1->SetTitle(";W_{#gammap} (GeV);#sigma(#gamma+p #rightarrow J/#psi+p) (nb)");
    Float_t siz = 0.045;
    frame1->SetTitleSize(siz);       frame1->SetLabelSize(siz);
    frame1->SetTitleSize(siz, "Y");  frame1->SetLabelSize(siz, "Y");
    frame1->GetXaxis()->SetTitleOffset(1.4);
    frame1->GetYaxis()->SetTitleOffset(1.0);
    
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
    //bxtit->DrawLatex(gxmax, gymax+450, "Bjorken-#it{x}");
    double yPosBjorken = pow(10., TMath::Log10(gymax) + 0.1* (TMath::Log10(gymax)-TMath::Log10(gymin)));
    bxtit->DrawLatex(gxmax, yPosBjorken, "Bjorken-#it{x}");
    
    // Draw the logo
    //  0: Just "ALICE" (for final data), to be added only if ALICE does not appear otherwise (e.g. in the legend)
    //  >0: ALICE Preliminary
    //DrawLogo(1, 0.58, 0.81);
    
    //end of layout configuration
    
    //alice fit
    alice *alic = alice::instance();
    
    TGraph *gShade = alic->fit();
    gShade->Draw("fsame");
    
    //theoretical models
    //STARLIGHT
    TF1* fStarlight = new TF1("fStarlight","((1.0-( (4.035*4.035)/(x*x) ))**2.0)*4.1*exp(0.65*log(x))",4.1,2000);
    fStarlight->SetLineStyle(3);
    fStarlight->SetLineWidth(2);
    
    //JMRT
    TGraphErrors* gJMRT_LO  = read_jmrt("JMRT-LO-JPSI.txt");
    TGraphErrors* gJMRT_NLO = read_jmrt("JMRT-NLO-JPSI.txt");
    SetGraphStyle(gJMRT_LO, 1, kDot, 1, kRed, 7, 2);
    SetGraphStyle(gJMRT_NLO, 1, kDot, 1, kBlack, 5, 2);
    //SetGraphStyle(gJMRT_NLO, 1, kDot, 1, kOrange, 5, 2);
    
    //model by Martin Hentschinski
    // http://inspirehep.net/record/1476661?ln=en
    TGraph *gMH_BFKL = read_mh_bfkl("martin_H_BFKL_refined.tab");
    
    // CCT
    TGraphErrors* gCCT = read_cct("Fig3_exc.txt");
    SetGraphStyle(gCCT, kBlue, kOpenCircle, 0, kBlue, 1, 2);
    gCCT->Draw("lsame");
    
    // JIMWLK
    TGraphErrors* gJIMWLK = read_jimwlk();
    SetGraphStyle(gJIMWLK, kBlue, kOpenCircle, 0, kBlack, 1, 2);
    //gJIMWLK->Draw("lsame");
    
    // Draw all models
    if (drawModels) {
        fStarlight->Draw("same");
        gJMRT_LO->Draw("cxsame");
        //gJMRT_NLO->Draw("cxsame");
        gMH_BFKL->Draw("fsame");
    }
    gJMRT_NLO->Draw("cxsame");
    
    //model by N. Armesto and A. Rezaeian
    // http://inspirehep.net/record/1282026?ln=en
    Bool_t pCGC = kTRUE;
    TGraph **gCGC = NULL, *gIP = NULL, *gbCGC = NULL;
    if( pCGC ) {
        //CGC
        gCGC = read_amir_cgc("data-jpsi-band-w_Amir.txt");
        if (drawModels) {
            gCGC[0]->Draw("fsame");
            gCGC[1]->Draw("lsame");
            gCGC[2]->Draw("lsame");
        }
    } else {
        //IP-Sat and bCGC
        TGraph **gIPb = read_amir_ipsat_bcgc("data-jpsi-w_Amir.txt");
        gIP = gIPb[0];
        gbCGC = gIPb[1];
        SetGraphStyle(gIP, 1, kDot, 1, kGreen+2, 9, 2);
        SetGraphStyle(gbCGC, 1, kDot, 1, kOrange-2, 10, 2);
        if (drawModels) {
            gIP->Draw("lsame");
            gbCGC->Draw("lsame");
        }
    }
    
    //experiments
    //H1
    TGraphErrors* gH1 = read_h1("d13-058.table_w_allH1_v1.txt");
    SetGraphStyle(gH1, kBlue, kFullCircle, 1, kBlue, 1, 1);
    TGraphErrors* gH1_2005 = read_h1_2005();
    SetGraphStyle(gH1_2005, kBlue, kFullCircle, 1, kBlue, 1, 1);
    if (drawPreviousMeas) {
        //gH1->Draw("pzsame");
        gH1_2005->Draw("pzsame");
    }
    
    //ZEUS
    TGraphErrors* gZeus_ee = read_zeus_ee();
    TGraphErrors* gZeus_mm = read_zeus_mm();
    SetGraphStyle(gZeus_ee, kBlue, kFullTriangleDown, 1, kBlue, 1, 1);
    SetGraphStyle(gZeus_mm, kBlue, kFullTriangleDown, 1, kBlue, 1, 1);
    if (drawPreviousMeas) {
        gZeus_ee->Draw("pzsame");
        gZeus_mm->Draw("pzsame");
    }
    
    //LHCb
    // @ 7TeV
    TGraphErrors* gLHCbPos = read_lhcb("lhcb_data2.txt",20,0);
    TGraphErrors* gLHCbNeg = read_lhcb("lhcb_data2.txt",20,1);
    //SetGraphStyle(gLHCbPos     ,kBlack,kOpenCircle      ,1.4,kBlack , 1,1);
    SetGraphStyle(gLHCbPos, kBlack, kOpenSquare, 1, kBlack, 1, 1);
    SetGraphStyle(gLHCbNeg, kBlack, kOpenSquare, 1, kBlack, 1, 1);
    // @ 13 TeV
    TGraphErrors* gLHCbPosNew = read_lhcb_pos_new();
    TGraphErrors* gLHCbNegNew = read_lhcb_neg_new();
    SetGraphStyle(gLHCbPosNew, kBlack, kOpenTriangleDown, 1, kBlack, 1, 1);
    SetGraphStyle(gLHCbNegNew, kBlack, kOpenTriangleDown, 1, kBlack, 1, 1);
    if (drawPreviousMeas) {
        gLHCbPos->Draw("pzsame");
        gLHCbNeg->Draw("pzsame");
        gLHCbPosNew->Draw("pzsame");
        gLHCbNegNew->Draw("pzsame");
    }
    
    // fixed target experiments
    TGraphAsymmErrors* gE687 = read_E687();
    TGraphErrors* gE401 = read_E401();
    TGraphErrors* gE516 = read_E516();
    SetGraphStyle(gE687, kBlue, kOpenDiamond, 1, kBlue, 1, 1);
    SetGraphStyle(gE516, kBlue, kOpenDiamond, 1, kBlue, 1, 1);
    SetGraphStyle(gE401, kBlue, kOpenDiamond, 1, kBlue, 1, 1);
    if (drawPreviousMeas) {
        gE687->Draw("pzsame");
        gE516->Draw("pzsame");
        gE401->Draw("pzsame");
    }
    
    //alice data
    TGraphErrors *aliceFwd      = alic->getForward();
    SetGraphStyle(aliceFwd, kOrange+1, kFullCircle, 1., kOrange+1, 1, 2);
    TGraphErrors *aliceNew      = alic->getNewData();
    SetGraphStyle(aliceNew, kOrange+1, kFullCircle, 1., kOrange+1, 1, 2);
    TGraphErrors *alice8T       = alic->getForwardEightTeV();
    SetGraphStyle(alice8T, kRed, kFullCircle, 1.5, kRed, 1, 2);
    if (drawPreviousMeas) {
        aliceFwd->Draw("psame");
        aliceNew->Draw("psame");
    }
    alice8T->Draw("psame");
    
    //legend for experimental data
    Double_t xl = 0.09, dxl = 0.4;
    Double_t yl = 0.51, dyl = 0.35;
    TLegend* leg1 = new TLegend(xl, yl, xl+dxl, yl+dyl);
    leg1->SetFillStyle(0);
    leg1->SetBorderSize(0);
    leg1->SetTextSize(siz-0.005);
    leg1->AddEntry(alice8T,"ALICE (#sqrt{#it{s}_{NN}} = 8.16 TeV)","p");
    if (drawPreviousMeas) {
        //leg1->AddEntry(aliceNew,"ALICE","p");
        //leg1->AddEntry(aliceFwd,"ALICE (PRL113 (2014) 232504)","p");
        leg1->AddEntry(aliceFwd,"ALICE (#sqrt{#it{s}_{NN}} = 5.02 TeV)","p");
        //leg1->AddEntry(gShade,"Power-law fit to ALICE data","lf");
        leg1->AddEntry(gE687,"Fixed target exp. (E401, E516, E687)","p");
        leg1->AddEntry(gH1,"H1","p");
        leg1->AddEntry(gZeus_ee,"ZEUS","p");
        /*
         leg1->AddEntry(gLHCbPos,"LHCb pp (W+ solutions)","p");
         leg1->AddEntry(gLHCbNeg,"LHCb pp (W- solutions)","p");
         */
        leg1->AddEntry(gLHCbNeg,"LHCb pp (#sqrt{#it{s}} = 7 TeV)","p");
        leg1->AddEntry(gLHCbNegNew,"LHCb pp (#sqrt{#it{s}} = 13 TeV)","p");
    }
    leg1->Draw("same");
    
    //legend with theoretical models
    if(pCGC) {
        xl = 0.61, dxl = 0.46;
        yl = 0.16, dyl = 0.30;
    } else {
        xl = 0.63, dxl = 0.46;
        yl = 0.16, dyl = 0.33;
    }
    if (!drawModels) {
        xl = 0.44, dyl = 0.21;
    }
    TLegend* leg2 = new TLegend(xl, yl, xl+dxl, yl+dyl);
    leg2->SetFillStyle(0);
    leg2->SetBorderSize(0);
    leg2->SetTextSize(siz-0.005);
    if (drawModels) {
        leg2->AddEntry(gJMRT_LO,"JMRT LO","l");
        //leg2->AddEntry(gJMRT_NLO,"JMRT NLO","l");
        leg2->AddEntry(fStarlight,"STARLIGHT param.","l");
        leg2->AddEntry(gMH_BFKL, "NLO BFKL", "f");
        if(pCGC) {
            leg2->AddEntry(gCGC[0], "CGC (IP-Sat, b-CGC)", "f");
        } else {
            leg2->AddEntry(gIP, "IP-Sat", "l");
            leg2->AddEntry(gbCGC, "b-CGC", "l");
        }
    }
    leg2->AddEntry(gJMRT_NLO,"JMRT NLO","l");
    leg2->AddEntry(gCCT,"CCT","l");
    //leg2->AddEntry(gJIMWLK, "JIMWLK", "l");
    leg2->AddEntry(gShade,"Power-law fit to ALICE data","lf");
    leg2->Draw("same");
    
    //put vertical scale atop the models
    gPad->RedrawAxis("Y");
    
    if (drawPreviousMeas) c1->SaveAs("fig-exc-final.eps");
    else c1->SaveAs("fig-exc-new-only-final.eps");
    
    return 0;
    
}//main
