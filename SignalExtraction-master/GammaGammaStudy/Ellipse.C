#include <iostream>
#include <fstream>
#include <dirent.h>
#include <string>
#include <cmath>

#include "TSystem.h"
#include "TROOT.h"
#include "TMath.h"
#include "TString.h"
#include "TDatime.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TStyle.h"
#include "TColor.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TText.h"
#include "TLatex.h"
#include "TLegend.h"

using namespace std;

double Square(double x) {
    return x*x;
}

int Ellipse() {
    
    const int n = 100;       // number of points
    double xList[2*n], yList[2*n];
    
    // mean values
    const double mean1 = 0;
    const double mean2 = 0;
    // root mean square
    const double sigma1 = 20;
    const double sigma2 = 8;
    
    const double corr = 0.5;
    
    double size = 2;    // size of the TH2
    TH2* h = new TH2F("h", "h", size*n, mean1-n*size/2, mean1+n*size/2, size*n, mean2-n*size/2, mean2+n*size/2);
    TH2* h2 = new TH2F("h2", "h2", size*n, mean1-n*size/2, mean1+n*size/2, size*n, mean2-n*size/2, mean2+n*size/2);
    
    TH1* hx1 = new TH1F("hx1", "hx1", size*n, mean1-n*size/2, mean1+n*size/2);
    TH1* hx2 = new TH1F("hx2", "hx2", size*n, mean1-n*size/2, mean1+n*size/2);
    
    TH1* hy1 = new TH1F("hy1", "hy1", size*n, mean2-n*size/2, mean2+n*size/2);
    TH1* hy2 = new TH1F("hy2", "hy2", size*n, mean2-n*size/2, mean2+n*size/2);
    
    // normalisation coeff
    double tot1 = 0;
    double tot2 = 0;
    for (int i = int(mean1-n*size/2); i < int(mean1+n*size/2); i++) {
        for (int j = int(mean2-n*size/2); j < int(mean2+n*size/2); j++) {
            tot1 += TMath::Exp(-Square(mean1-i)/(2*Square(sigma1)))*TMath::Exp(-Square(mean2+(corr*i*sigma2/sigma1)-j)/(2*Square(sigma2*(1-abs(corr)))));
            tot2 += TMath::Exp(-Square(mean1+ (corr*j*sigma1/sigma2)-i)/(2*Square(sigma1*(1-abs(corr)))))*TMath::Exp(-Square(mean2-j)/(2*Square(sigma2)));
        }
    }
    
    std::cout << "tot1 = " << tot1 << std::endl;
    std::cout << "tot2 = " << tot2 << std::endl;
    
    // define points of ellipse
    for (int i = int(mean1-n*size/2); i < int(mean1+n*size/2); i++) {
        for (int j = int(mean2-n*size/2); j < int(mean2+n*size/2); j++) {
            double val1 = TMath::Exp(-Square(mean1-i)/(2*Square(sigma1)));
            double val2 = TMath::Exp(-Square(mean2+(corr*i*sigma2/sigma1)-j)/(2*Square(sigma2*(1-abs(corr)))));
            //double value = 100* TMath::Exp(-Square(mean1-i)/(2*Square(sigma1)))*TMath::Exp(-Square(mean2+(corr*i*sigma2/sigma1)-j)/(2*Square(sigma2*(1-abs(corr)))));
            double value = 100* val1*val2;
            if (value > 1)
            h->Fill(i, j, value);
            hx1->Fill(i, val1);
            hy1->Fill(j, val2);
            double val3 = TMath::Exp(-Square(mean1+ (corr*j*sigma1/sigma2)-i)/(2*Square(sigma1*(1-abs(corr)))));
            double val4 = TMath::Exp(-Square(mean2-j)/(2*Square(sigma2)));
            double value2 = 100* val3*val4;
            if (value2 > 1)
            h2->Fill(i, j, value2);
            hx2->Fill(i, val3);
            hy2->Fill(j, val4);
        }
    }
    
    // Draw the TGraph
    TCanvas* cv = new TCanvas("cv", "cv", 600, 300);
    cv->Divide(4,2);
    cv->cd(1);
    hy1->Draw("hist");
    cv->cd(2);
    h->Draw("colz");
    cv->cd(3);
    h2->Draw("colz");
    cv->cd(4);
    hy2->Draw("hist");
    cv->cd(6);
    hx1->Draw("hist");
    cv->cd(7);
    hx2->Draw("hist");
    //cv->SaveAs("ellipse.pdf");
    
    // real weights
    TH2* h3 = new TH2F("h3", "h3", size*n, mean1-n*size/2, mean1+n*size/2, size*n, mean2-n*size/2, mean2+n*size/2);
    
    for (int i = int(mean1-n*size/2); i < int(mean1+n*size/2); i++) {
        for (int j = int(mean2-n*size/2); j < int(mean2+n*size/2); j++) {
            double term1 = Square(mean1-i)/Square(sigma1);
            double term2 = Square(mean2-j)/Square(sigma2);
            double term3 = (mean1-i)*(mean2-j)/(sigma1*sigma2);
            double value = TMath::Exp(-1./2*(term1 + term2 - 2*corr *term3)/(1-corr*corr));
            //double value2 = TMath::Exp(-1./2*(term1 + term2));
            double norm = 2*TMath::Pi()*sigma1*sigma2*TMath::Sqrt(1-corr*corr);
            if (value/norm>0.000001) h3->Fill(i, j, value/norm);
        }
    }
    
    TCanvas* cv2 = new TCanvas("cv2", "cv2", 300, 300);
    h3->Draw("colz");
    cv2->SaveAs("ellipse2.pdf");
    
    return 0;
}
