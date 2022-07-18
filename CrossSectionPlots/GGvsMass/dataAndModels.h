
//global constants with range of the plot
const Double_t gxmin = 0.8, gymin = 0.03, gxmax = 2.7, gymax = 20;

//_____________________________________________________________________________
void SetGraphStyle(TGraph* g, Color_t mcolor, Style_t mstyle, Size_t msize, Color_t lcolor, Style_t lstyle, Width_t lwidth){
  g->SetMarkerColor(mcolor);
  g->SetMarkerSize(msize);
  g->SetMarkerStyle(mstyle);
  g->SetLineColor(lcolor);
  g->SetLineStyle(lstyle);
  g->SetLineWidth(lwidth);
}//SetGraphStyle

void ifcnChi2ModelA(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);//declaration of interface to chi2 function

//_____________________________________________________________________________
class alice {
  Int_t opt;//option for the fit
  static const Int_t ALICE_2013_n = 9;      // fw     fw   fw    sfw     sfw    mid      mid   sfw      fw
  const Double_t ALICE_2013_W[ALICE_2013_n] = {24.1, 30.9, 39.6, 50.4,  73.1,  129.9,  193.3,  391.2,  706.};
  const Double_t ALICE_2013_sig[ALICE_2013_n] = {27.6, 34.7, 38.5, 41.9,  61.9,  101.4,  126.2,  194.5,  284.};
  const Double_t ALICE_2013_sta[ALICE_2013_n] = { 3.6,  3.1,  5.6,  8.6,  10.8,    8.1,   10.0,   27.4,   36};
  // model A, considers MUON and trigger efficiencyes uncorrelated
  static const Int_t nUncA = 10; //  sources of uncorrelated error
  const Double_t ALICE_2013_unc_A[nUncA][ALICE_2013_n] = // relative uncertainty!
  {
    {0, 0, 0, 0.057, 0.012, 0.009, 0.008, 0.033, 0}, // TPC track selection
    {0, 0, 0, 0, 0, 0.013, 0.006, 0, 0}, // PID efficiency
    {0.04, 0.04, 0.04, 0.02, 0.02, 0, 0, 0.03, 0.06}, // muon track eff
    {0.01, 0.01, 0.01, 0.005, 0.005, 0, 0, 0.005, 0.01}, // muon matching eff
    {0.028, 0.028, 0.028, 0.02, 0.02, 0, 0 , 0.016, 0.032}, // muon trig eff
    {0, 0, 0, 0, 0, 0.101, 0.101, 0, 0}, // central trigger eff
    {0, 0, 0, 0, 0, 0, 0, 0.034, 0.035}, // V0C efficiency
    {0.066, 0.037, 0.08, 0.036, 0.022, 0.021, 0.019, 0.03, 0.06}, // signal extraction, last number should be -0,+0.06
    {0.031, 0.02, 0.02, 0.013, 0.01, 0.01, 0.01, 0.014, 0.031}, // feed-down
    //lumi at mid-rapidity by quadratic average:
    //sqrt((3.3%*2.1)**2 + (3%*4.8)**2) / (2.1 + 4.8) = 2.3%
    {0.033, 0.033, 0.033, 0.033, 0.033, 0.023, 0.023, 0.03, 0.03} // lumi uncorrelated
  };
  Double_t ALICE_2013_unc[ALICE_2013_n]; // to store total uncorrelated uncertainty
  static const Int_t nCorA = 3; // sources of correlated error
  const Double_t ALICE_2013_corr_A[nCorA][ALICE_2013_n] = // relative uncertainty!
  {
    {0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016}, // lumi
    {0.006, 0.006, 0.006, 0.006, 0.006, 0.004, 0.004, 0.006, 0.006}, // br
    {0.02, 0.02, 0.02, 0.027, 0.035, 0.021, 0.021, 0.012, 0.02}//, // v0a veto
    //{0, 0, 0, 0, 0, 0.09, 0.09, 0, 0} // TOF trigger, Eur. Phys. J. C (2013) 73:2617, now in central trigger eff
  };
  Double_t ALICE_2013_tot[ALICE_2013_n]; // to store total error
  Double_t ALICE_2013_statsys[ALICE_2013_n]; // to store stat + total sys error for the plot
  Double_t Nval, Dval; // power-law parameters for further use
  TGraph *fitGraph; // to store fit result
  void TotalUncA();
  void TotalErrorA();
  void StatSysError();
  TGraph* GetShade(Double_t N, Double_t d, Double_t cov00, Double_t cov11, Double_t cov01);
  static alice alic;
  alice(Int_t oo=0);
  alice(const alice&);
public:
  static alice* instance() {return &alic;}
  void fcnChi2ModelA(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
  TGraph *fit() const {return fitGraph;}
  Double_t getPowerLaw(Double_t wval) const {return Nval*pow(wval/90.,Dval);}
  void getPowerLaw(Double_t &nn, Double_t &dd) const {nn = Nval; dd = Dval;}
  TGraphErrors *getForward();
  TGraphErrors *getForwardEightTeV1();
    TGraphErrors *getForwardEightTeV2();
    TGraphErrors *getForwardEightTeV3();
    TGraphErrors *getForwardEightTeVSL1();
    TGraphErrors *getForwardEightTeVSL2();
    TGraphErrors *getForwardEightTeVSL3();
  void printData();
};//alice
alice alice::alic;

alice::alice(Int_t oo) : opt(oo), Nval(-1), Dval(-1), fitGraph(0) {

  for(Int_t i=0; i<ALICE_2013_n; i++) {ALICE_2013_unc[i]=0.; ALICE_2013_tot[i]=0.; ALICE_2013_statsys[i]=0.;}

  if (opt == 0) {
    TotalUncA();
    TotalErrorA();
    StatSysError();
  } else {
    cout << " Option not recognized. Bye!" << endl;
    exit(1);
  }

}//alice

void alice::TotalUncA()
{
// compute total uncorrelated uncertainty
  for(Int_t j=0;j<ALICE_2013_n;j++) {
    Double_t unc = 0;
    for(Int_t i=0;i<nUncA;i++) unc+=TMath::Power(ALICE_2013_unc_A[i][j],2);
    ALICE_2013_unc[j]=TMath::Sqrt(unc);
    //  cout << " Total uncorrelated uncertainty for measurement " << j << " is " << ALICE_2013_unc[j] << endl;
  }
}
void alice::TotalErrorA()
{
  // to print also sys err on differential cross section
  Double_t dSigmaDy[] = {6.9, 8.7, 10.6, 10.0, 7.1};
  // compute total  uncertainty
  for(Int_t j=0;j<ALICE_2013_n;j++) {
    Double_t tot = (ALICE_2013_unc[j]*ALICE_2013_unc[j]);
    for(Int_t i=0;i<nCorA;i++) tot+=TMath::Power(ALICE_2013_corr_A[i][j],2);
    ALICE_2013_tot[j]=TMath::Sqrt(tot)*ALICE_2013_sig[j];
    cout << " Total uncertainty for measurement " << j << " is ";
    cout << fixed << setprecision(1) << ALICE_2013_tot[j] << "  (" << 100*TMath::Sqrt(tot) << "%)";
    if(j > 2 && j < 8) {
      cout << ",   Delta dSigma/dy = " << TMath::Sqrt(tot)*dSigmaDy[j-3];
    }
    cout << endl;
  }
}
void alice::StatSysError()
{
// quadratic sum of statistical uncertainty and total systematic uncertainty, to be used for the plot
  for(Int_t i=0; i<ALICE_2013_n; i++) {
    Double_t err = 0.;
    err = TMath::Power(ALICE_2013_sta[i], 2) + TMath::Power(ALICE_2013_tot[i], 2);
    ALICE_2013_statsys[i] = TMath::Sqrt(err);
  }
}
TGraph* alice::GetShade(Double_t N, Double_t d, Double_t cov00, Double_t cov11, Double_t cov01)
{
// Define graph with error band from results of the fit
  Float_t wmin=gxmin;
  Float_t wmax=gxmax;
  Int_t nw = 200;
  Color_t color = 29;
  Double_t dw = (wmax-wmin)/nw;
  TGraph *gshade = new TGraph(2*nw+2);
  for (Int_t i=0;i<=nw;i++){
    Double_t w = wmin+i*dw;
    Double_t y = N*pow(w/90,d);
    Double_t dydN = y/N;
    Double_t dydd = y*log(w/90);
    Double_t dy = sqrt(dydN*dydN*cov00 + dydd*dydd*cov11 + 2*dydN*dydd*cov01);
    gshade->SetPoint(i       ,w,y+dy);
    gshade->SetPoint(2*nw-i+1,w,y-dy);
  }
  gshade->SetFillColor(color);
  gshade->SetLineColor(color);
  return gshade;
}
inline void alice::fcnChi2ModelA(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
// chi2 function to fit data according to prescription in
// Eur.Phys.J. C63 (2009) 625-678, section 9, eq (31)

  // ch2 = sum_i [m_i-mu_i-Sij]^2/D + Sbj
  // Sij = sum_j g_ij*m_i*b_j
  // Sbj = sum_j (b_j)^2
  // D = d_i_stat^2*mu_i*(m_i-Sij)+(d_i,unc*m_i)^2
  // mu_i is the measurement
  // g_ij relative normalization uncertainty at point i from source j
  // d_i = relative uncertainty (either stat or uncorr)

  //Double_t bj[nCorA] = {par[2],par[3],par[4],par[5]};
  Double_t bj[nCorA];
  for(Int_t i=0; i<nCorA; i++) {
    bj[i] = par[i+2];
  }

  Double_t Sbj = 0;
  for(Int_t j=0;j<nCorA;j++) Sbj += (bj[j]*bj[j]);
  Double_t chi2 = 0;
  for (Int_t i=0; i<ALICE_2013_n;i++) {
    Double_t mu_i = ALICE_2013_sig[i];
    Double_t m_i = par[0]*TMath::Power(ALICE_2013_W[i]/90,par[1]);
    Double_t Sij = 0;
    for(Int_t j=0;j<nCorA;j++) Sij += m_i*bj[j]*ALICE_2013_corr_A[j][i];
    Double_t d_stat_i =  ALICE_2013_sta[i]/ALICE_2013_sig[i];
    Double_t d_unc_i = ALICE_2013_unc[i];
    Double_t D = (d_stat_i*d_stat_i*mu_i*(m_i-Sij))+((d_unc_i*m_i)*(d_unc_i*m_i));
    chi2 += ((m_i-mu_i-Sij)*(m_i-mu_i-Sij)/D);
  }
  chi2 += Sbj;
  f = chi2;
}

TGraphErrors *alice::getForward() {
  //forward data
  const Int_t np = 4;
  const Double_t w[np] = {ALICE_2013_W[0], ALICE_2013_W[1], ALICE_2013_W[2], ALICE_2013_W[8]};
  const Double_t sig[np] = {ALICE_2013_sig[0], ALICE_2013_sig[1], ALICE_2013_sig[2], ALICE_2013_sig[8]};
  const Double_t err[np] = {ALICE_2013_statsys[0], ALICE_2013_statsys[1], ALICE_2013_statsys[2], ALICE_2013_statsys[8]};

  TGraphErrors *gr = new TGraphErrors(np, w, sig, NULL, err);
  SetGraphStyle(gr ,kBlack  ,kFullDiamond     ,2.,kBlack   , 1,2);

  return gr;

}//getForward


// Analysis with Run 2 settings
// 1.0-1.5
TGraphErrors *alice::getForwardEightTeV1() {
    //forward data
    const Int_t np = 3;
    
    const Double_t massMean[np]   = {1.25, 1.75, 2.25};
    const Double_t massInterval[np]   = {0.25, 0.25, 0.25};
    
    const Double_t rapInterval[np]   = {0.75, 0.75, 0.75};
    const Double_t rapMean[np]   = {3.625, 3.625, 3.625};
    const Double_t sig[np] = { 2.200, 0.820, 0.270};
    const Double_t err[np] = {0.200, 0.070, 0.040};
    Double_t dsig[np];
    for (int i = 0; i < np; i++) {dsig[i] = sig[i]/(2*massInterval[i]);}
    Double_t derr[np];
    for (int i = 0; i < np; i++) {derr[i] = err[i]/(2*massInterval[i]);}
    TGraphErrors *gr = new TGraphErrors(np, massMean, dsig, massInterval, derr);
    //TGraphErrors *gr = new TGraphErrors(np, massMean, sig, massInterval, err);
    SetGraphStyle(gr, kRed, kFullDiamond, 2., kRed, 1, 2);
    
    return gr;
    
}//getForward

// Analysis with Run 2 settings
// 1.5-2.0
TGraphErrors *alice::getForwardEightTeV2() {
    //forward data
    const Int_t np = 3;
    
    const Double_t massMean[np]   = {1.25, 1.75, 2.25};
    const Double_t massInterval[np]   = {0.25, 0.25, 0.25};
    
    const Double_t rapInterval[np]   = {0.75, 0.75, 0.75};
    const Double_t rapMean[np]   = {2.875, 2.875, 2.875};
    const Double_t sig[np] = {3.000, 1.100, 0.360};
    //const Double_t err[np] = {  3.300,                     2.500,                      4.700                      };
    const Double_t err[np] = {0.400, 0.100, 0.060};
    Double_t dsig[np];
    for (int i = 0; i < np; i++) {dsig[i] = sig[i]/(2*massInterval[i]);}
    Double_t derr[np];
    for (int i = 0; i < np; i++) {derr[i] = err[i]/(2*massInterval[i]);}
    TGraphErrors *gr = new TGraphErrors(np, massMean, dsig, massInterval, derr);
    //TGraphErrors *gr = new TGraphErrors(np, massMean, sig, massInterval, err);
    SetGraphStyle(gr, kBlue, kFullDiamond, 2., kBlue, 1, 2);
    
    return gr;
    
}//getForward




TGraphErrors *alice::getForwardEightTeVSL1() {
    //forward data
    const Int_t np = 3;
    
    const Double_t massMean[np]   = {1.25, 1.75, 2.25};
    const Double_t massInterval[np]   = {0.25, 0.25, 0.25};
    
    const Double_t rapInterval[np]   = {0.75, 0.75, 0.75};
    const Double_t rapMean[np]   = {3.625, 3.625, 3.625};
    const Double_t sig[np] = {2.000, 0.720, 0.340};
    //const Double_t err[np] = {  3.300,                     2.500,                      4.700                      };
    const Double_t err[np] = {0.010, 0.010, 0.010};
    Double_t dsig[np];
    for (int i = 0; i < np; i++) {dsig[i] = sig[i]/(2*massInterval[i]);}
    Double_t derr[np];
    for (int i = 0; i < np; i++) {derr[i] = err[i]/(2*massInterval[i]);}
    TGraphErrors *gr = new TGraphErrors(np, massMean, dsig, massInterval, derr);
    //TGraphErrors *gr = new TGraphErrors(np, massMean, sig, massInterval, err);
    SetGraphStyle(gr, kRed, kOpenDiamond, 2., kRed, 4, 2);
    
    return gr;
    
}//getForward

// Analysis with Run 2 settings
// StarLIGHT 1.5-2.0
TGraphErrors *alice::getForwardEightTeVSL2() {
    //forward data
    const Int_t np = 3;
    
    const Double_t massMean[np]   = {1.25, 1.75, 2.25};
    const Double_t massInterval[np]   = {0.25, 0.25, 0.25};
    
    const Double_t rapInterval[np]   = {0.75, 0.75, 0.75};
    const Double_t rapMean[np]   = {2.875, 2.875, 2.875};
    const Double_t sig[np] = {2.200, 0.780, 0.370};
    const Double_t err[np] = {0.010, 0.010, 0.010};
    
    Double_t dsig[np];
    for (int i = 0; i < np; i++) {dsig[i] = sig[i]/(2*massInterval[i]);}
    Double_t derr[np];
    for (int i = 0; i < np; i++) {derr[i] = err[i]/(2*massInterval[i]);}
    
    TGraphErrors *gr = new TGraphErrors(np, massMean, dsig, massInterval, derr);
    //TGraphErrors *gr = new TGraphErrors(np, massMean, sig, massInterval, err);
    SetGraphStyle(gr, kBlue, kOpenDiamond, 2., kBlue, 4, 2);
    
    return gr;
    
}//getForward

void DrawLogo (Int_t logo, Double_t xmin, Double_t ymin) {
    
    // Logo is not needed anymore, now we only write alice preliminary
    // Logo:
    // 0: Justr writes "ALICE" (for final data)
    // Anything eles: writes "ALICE Preliminary"
    
    TLatex *   tex = new TLatex(xmin,ymin, logo ? "ALICE Preliminary" : "ALICE");
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->Draw();
}

