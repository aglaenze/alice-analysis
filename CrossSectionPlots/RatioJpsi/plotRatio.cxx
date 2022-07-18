
/*

  plot of sigma(gamma-p) for semiforward and mid-rapidity paper, Alice data,
  other experiments and models

  compile:

g++ plotRatio.cxx `root-config --libs --cflags` -lEG -lMinuit -o plotRatio

  run:

./plotRatio

*/

// c++ headers
#include <iostream>
#include <fstream>
#include <string>
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
int main(void) {

  //configure layout of the plot
  gStyle->SetPadTickY(1);
  gStyle->SetTickLength(0.02,"Y");

  TCanvas* c1 = new TCanvas("c1","c1",800,700);

  Double_t pdiv = 0.4;
  TPad *pad1 = new TPad("p1", "p1", 0., pdiv, 1., 1.); // upper
  TPad *pad2 = new TPad("p2", "p2", 0., 0., 1., pdiv); // lower
  pad1->Draw();
  pad2->Draw();

  pad1->cd();

  //gPad = pad2;

  //gPad->Draw();

  TH1F* frame1 = gPad->DrawFrame(gxmin,gymin,gxmax,gymax);
  frame1->GetXaxis()->SetMoreLogLabels();
  gPad->SetLeftMargin(0.09);
  gPad->SetRightMargin(0.01);
  gPad->SetTopMargin(0.11);//0.025
  gPad->SetBottomMargin(0.);
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
  bxtit->DrawLatex(gxmax, gymax+400, "Bjorken-#it{x}");
  //end of layout configuration

  //alice fit
  alice *alic = alice::instance();

  TGraph *gShade = alic->fit();
  gShade->Draw("fsame");

  //theoretical models
  //STARLIGHT
  string starlight_param = "((1.0-( (4.035*4.035)/(x*x) ))**2.0)*4.1*exp(0.65*log(x))";
  TF1* fStarlight = new TF1("fStarlight",starlight_param.c_str(),4.1,2000);
  fStarlight->SetLineStyle(3);
  fStarlight->SetLineWidth(2);
  fStarlight->Draw("same");

  //JMRT
  TGraphErrors* gJMRT_LO  = read_jmrt("JMRT-LO-JPSI.txt");
  TGraphErrors* gJMRT_NLO = read_jmrt("JMRT-NLO-JPSI.txt");
  SetGraphStyle(gJMRT_LO  ,     1,kDot             ,  1,kRed   , 7,2);
  SetGraphStyle(gJMRT_NLO ,     1,kDot             ,  1,kBlue  , 5,2);
  //gJMRT_LO->Draw("cxsame");
  gJMRT_NLO->Draw("cxsame");

  //CGTT
  TGraphErrors *gCGTT = read_cgtt("guillermo_CGTT.txt");
  SetGraphStyle(gCGTT  ,     1,kDot             ,  1,kRed   , 7,2);
  gCGTT->Draw("cxsame");

  //model by Martin Hentschinski
  // http://inspirehep.net/record/1476661?ln=en
  TGraph *gMH_BFKL = read_mh_bfkl("martin_H_BFKL_refined.tab");
  gMH_BFKL->Draw("fsame");

  //model by N. Armesto and A. Rezaeian
  // http://inspirehep.net/record/1282026?ln=en
  Bool_t pCGC = kTRUE;
  TGraph **gCGC = NULL, *gIP = NULL, *gbCGC = NULL;
  if( pCGC ) {
    //CGC
    gCGC = read_amir_cgc("data-jpsi-band-w_Amir.txt");
    gCGC[0]->Draw("fsame");
    gCGC[1]->Draw("lsame");
    gCGC[2]->Draw("lsame");
  } else {
    //IP-Sat and bCGC
    TGraph **gIPb = read_amir_ipsat_bcgc("data-jpsi-w_Amir.txt");
    gIP = gIPb[0];
    gbCGC = gIPb[1];
    SetGraphStyle(gIP ,     1,kDot             ,  1,kGreen+2, 9,2);
    SetGraphStyle(gbCGC ,   1,kDot             ,  1,kOrange-2,10,2);
    gIP->Draw("lsame");
    gbCGC->Draw("lsame");
  }

  //experiments
  //H1
  TGraphErrors* gH1 = read_h1("d13-058.table_w_allH1_v1.txt");
  SetGraphStyle(gH1       ,kBlue,kOpenCircle      ,1.4,kBlue , 1,1);
  gH1->Draw("pzsame");

  //ZEUS
  TGraphErrors* gZeus_ee = read_zeus_ee();
  TGraphErrors* gZeus_mm = read_zeus_mm();
  SetGraphStyle(gZeus_ee  ,kBlue,kOpenTriangleDown      ,1.4,kBlue , 1,1);
  SetGraphStyle(gZeus_mm  ,kBlue,kOpenTriangleDown      ,1.4,kBlue , 1,1);
  gZeus_ee->Draw("pzsame");
  gZeus_mm->Draw("pzsame");

  //LHCb
  TGraphErrors* gLHCbPos = read_lhcb("lhcb_data2.txt",20,0);
  TGraphErrors* gLHCbNeg = read_lhcb("lhcb_data2.txt",20,1);
  SetGraphStyle(gLHCbPos     ,kBlack,kOpenCircle      ,1.4,kBlack , 1,1);
  SetGraphStyle(gLHCbNeg     ,kBlack,kOpenSquare      ,1.4,kBlack , 1,1);
  gLHCbPos->Draw("pzsame");
  gLHCbNeg->Draw("pzsame");

  //alice data
  TGraphErrors *aliceFwd = alic->getForward();
  aliceFwd->Draw("psame");
  TGraphErrors *aliceNew = alic->getNewData();
  alic->printData();
  aliceNew->Draw("psame");

  //legend for experimental data
  Double_t xl = 0.09, dxl = 0.4;
  Double_t yl = 0.51, dyl = 0.35;
  TLegend* leg1 = new TLegend(xl, yl, xl+dxl, yl+dyl);
  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);
  leg1->SetTextSize(siz-0.005);
  leg1->AddEntry(aliceNew,"ALICE","p");
  leg1->AddEntry(aliceFwd,"ALICE (PRL113 (2014) 232504)","p");
  leg1->AddEntry(gShade,"Power-law fit to ALICE data","lf");
  leg1->AddEntry(gH1,"H1","p");
  leg1->AddEntry(gZeus_ee,"ZEUS","p");
  leg1->AddEntry(gLHCbPos,"LHCb pp (W+ solutions)","p");
  leg1->AddEntry(gLHCbNeg,"LHCb pp (W- solutions)","p");
  leg1->Draw("same");

  //legend with theoretical models
  if(pCGC) {
    xl = 0.65, dxl = 0.46;
    yl = 0.1, dyl = 0.30;
  } else {
    xl = 0.65, dxl = 0.46;
    yl = 0.1, dyl = 0.33;
  }
  TLegend* leg2 = new TLegend(xl, yl, xl+dxl, yl+dyl);
  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);
  leg2->SetTextSize(siz-0.005);
  //leg2->AddEntry(gJMRT_LO,"JMRT LO","l");
  leg2->AddEntry(gCGTT,"CCT","l");
  leg2->AddEntry(gJMRT_NLO,"JMRT NLO","l");
  leg2->AddEntry(fStarlight,"STARLIGHT param.","l");
  leg2->AddEntry(gMH_BFKL, "NLO BFKL", "f");
  if(pCGC) {
    leg2->AddEntry(gCGC[0], "CGC (IP-Sat, b-CGC)", "f");
  } else {
    leg2->AddEntry(gIP, "IP-Sat", "l");
    leg2->AddEntry(gbCGC, "b-CGC", "l");
  }
  leg2->Draw("same");

  //put vertical scale atop the models
  gPad->RedrawAxis("Y");

  pad2->cd();

  TH1F* frame2 = gPad->DrawFrame(gxmin,0.6,gxmax,1.4-1.e-5); // 0.1, 2
  frame2->Draw();
  frame2->GetXaxis()->SetMoreLogLabels();
  gPad->SetLeftMargin(0.09);
  gPad->SetRightMargin(0.01);
  gPad->SetTopMargin(0.);//0.025
  gPad->SetBottomMargin(0.17);
  gPad->SetLogx();
  //gPad->SetLogy();

  frame2->SetTitle(";W_{#gammap} (GeV);Models / fit to data");
  siz = siz*(1.+(1.-pdiv));
  frame2->SetTitleSize(siz);       frame2->SetLabelSize(siz);
  frame2->SetTitleSize(siz, "Y");  frame2->SetLabelSize(siz, "Y");
  frame2->GetXaxis()->SetTitleOffset(1.1);
  frame2->GetYaxis()->SetTitleOffset(0.64);

  //ratios of data (alice fit) to models

  //line at one
  TLine *unity = new TLine(gxmin, 1., gxmax, 1.);
  unity->SetLineColor(kBlack);//same color as alice fit band
  unity->SetLineWidth(1);

  //alice fit
  TGraph *gShade_ratio = getRatio( gShade );
  gShade_ratio->Draw("fsame");
  unity->Draw();

  //Starlight ratio to data
  Double_t nfit, dfit;
  alic->getPowerLaw(nfit, dfit);
  string starlight_ratio = starlight_param + Form("/(%.3f*pow(x/90.,%.3f))", nfit, dfit);
  TF1* fStarlightRatio = new TF1("fStarlight",starlight_ratio.c_str(),4.1,2000);
  fStarlightRatio->SetLineStyle(3);
  fStarlightRatio->SetLineWidth(2);
  fStarlightRatio->Draw("same");

  //ratio JMRT LO and NLO
  TGraph *gJMRT_LO_ratio = getRatio( gJMRT_LO );
  TGraph *gJMRT_NLO_ratio = getRatio( gJMRT_NLO );
  //gJMRT_LO_ratio->Draw("cxsame");
  gJMRT_NLO_ratio->Draw("cxsame");

  //ratio CGTT
  TGraph *gCGTT_ratio = getRatio( gCGTT );
  gCGTT_ratio->Draw("cxsame");

  //ratio for model by Martin Hentschinski
  TGraph *gMH_BFKL_ratio = getRatio( gMH_BFKL );
  gMH_BFKL_ratio->Draw("fsame");

  //ratio for model by N. Armesto and A. Rezaeian
  TGraph *gCGC_ratio[3] = {NULL}, *gIP_ratio = NULL, *gbCGC_ratio = NULL;
  if( pCGC ) {
    //CGC
    for(Int_t i=0; i<3; i++) gCGC_ratio[i] = getRatio( gCGC[i] );
    gCGC_ratio[0]->Draw("fsame");
    gCGC_ratio[1]->Draw("lsame");
    gCGC_ratio[2]->Draw("lsame");
  } else {
    //IP-Sat and bCGC
    gIP_ratio = getRatio( gIP );
    gbCGC_ratio = getRatio( gbCGC );
    gIP_ratio->Draw("lsame");
    gbCGC_ratio->Draw("lsame");
  }

  //cout << "######### " << frame2->GetYaxis()->FindBin(1.4) << endl;;
  //frame2->GetYaxis()->SetBinLabel(2, "ddd");

  gPad->RedrawAxis("Y"); // axis atop models

  c1->SaveAs("fig3_ratio.pdf");

  return 0;

}//main

























