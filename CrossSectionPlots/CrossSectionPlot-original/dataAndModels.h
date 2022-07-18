
//global constants with range of the plot
const Double_t gxmin = 17., gymin = 15., gxmax = 1990., gymax = 1000.;

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
  void makeFit();
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
  TGraphErrors *getForwardEightTeV();
  TGraphErrors *getNewData();
  void printData();
};//alice
alice alice::alic;

alice::alice(Int_t oo) : opt(oo), Nval(-1), Dval(-1), fitGraph(0) {

  for(Int_t i=0; i<ALICE_2013_n; i++) {ALICE_2013_unc[i]=0.; ALICE_2013_tot[i]=0.; ALICE_2013_statsys[i]=0.;}

  if (opt == 0) {
    TotalUncA();
    TotalErrorA();
    StatSysError();
    makeFit();
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
void alice::makeFit()
{

  cout << endl<< endl<< endl;
  if (opt == 0)   cout << " - o - o - o - o - o - o - o - o -   Model A   - o - o - o - o - o - o - o - o -" << endl;
  cout << endl<< endl<< endl;

  // initialize minuit with a maximum of parameters
  Int_t nPar = 0;
  if (opt == 0) nPar = 2+nCorA;
  TMinuit *myMinuit = new TMinuit(nPar);

  // set the function with the minimization process
  if (opt == 0) myMinuit->SetFCN(ifcnChi2ModelA);

  // define parameters
  myMinuit->DefineParameter(0,"N",70,1,0.0,1000.);
  myMinuit->DefineParameter(1,"#delta",0.7,0.1,0.1,10.);
  Char_t ParName[120];
  Int_t nCor = 0;
  if (opt == 0) nCor = nCorA;
  for(Int_t j=0;j<nCor;j++) {
    sprintf(ParName,"b_{%d}",j);
    myMinuit->DefineParameter(2+j,ParName,0.001,0.0001,-1,1.);
  }
  // migrad
  myMinuit->SetMaxIterations(500);
  myMinuit->Migrad();

  ///////////////////////////////////////////////////////////////////////////////
  // get results
  ///////////////////////////////////////////////////////////////////////////////
  Double_t Cov[nPar*nPar];
  myMinuit->mnemat(Cov,nPar);
  Double_t N = 0;  // normalization
  Double_t Nerr = 0;  // normalization error
  myMinuit->GetParameter(0,N,Nerr);
  Double_t D = 0;  // exponent
  Double_t Derr = 0;  // exponent error
  myMinuit->GetParameter(1,D,Derr);
  Double_t covNN = Cov[0];
  Double_t covDD = Cov[nPar+1];
  Double_t covND = Cov[1];
  Double_t rho = covND/TMath::Sqrt(covDD*covNN);
  //status of the fit
  Double_t fmin, fedm, errdef;
  Int_t npari, nparx, istat;
  myMinuit->mnstat(fmin, fedm, errdef, npari, nparx, istat);

  cout << "#####  Fit results  #####" << endl;
  cout << " N     = " << fixed << setprecision(2) << N << " +/- " << Nerr << endl;
  cout << " delta = " << D << " +/- " << Derr << endl;
  cout << " rho   = " << fixed << setprecision(2) << rho << endl;
  cout << " chi2  = " << fixed << setprecision(2) << fmin << endl;
  cout << " stat: " << istat << endl;

  delete myMinuit;

  Nval = N;
  Dval = D;

  fitGraph = GetShade(N,D,covNN,covDD,covND);
}

void ifcnChi2ModelA(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  //interface to chi2 function

  alice::instance()->fcnChi2ModelA(npar, gin, f, par, iflag);

}//interface to chi2 function

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

// Analysis with Run 1 settings
TGraphErrors *alice::getForwardEightTeV() {
  //forward data
  const Int_t np = 3;
  const Double_t w[np]   = { 42.000 /* full rapidity*/, 33.000 /* 3.25 < y < 4. */, 48.000 /* 2.5 < y < 3.25 */ };
  const Double_t sig[np] = { 40.400,                    33.700,                     48.500                      };
  const Double_t err[np] = {  3.300,                     2.500,                      4.700                      };

  TGraphErrors *gr = new TGraphErrors(np, w, sig, NULL, err);
  SetGraphStyle(gr ,kMagenta  ,kFullDiamond     ,2.,kMagenta   , 1,2);

  return gr;

}//getForward






TGraphErrors *alice::getNewData() {
  //semiforward and mid-rapidity data
  const Int_t np = 5;
  const Double_t *w = &ALICE_2013_W[3];
  const Double_t *sig = &ALICE_2013_sig[3];
  const Double_t *err = &ALICE_2013_statsys[3];
/*
  cout << "#####  SFW + MID data for plot ######" << endl;
  cout.width(6);
  cout << "Wgp";
  cout.width(8);
  cout << "sigma" << endl;
  cout << fixed << setprecision(1);
  for(Int_t i=0; i<np; i++) {
    cout.width(6);
    cout << w[i];
    cout.width(8);
    cout << sig[i];
    cout << " +/- ";
    cout.width(4);
    cout << err[i] << endl;
  }
  cout << "#####################################" << endl;
*/
  TGraphErrors *gr = new TGraphErrors(np, w, sig, NULL, err);
  SetGraphStyle(gr ,kRed  ,kFullCircle     ,1.3,kRed   , 1,2);

  return gr;

}//getForward

void alice::printData() {
  //print all alice data data to stdout
  const Int_t np = 9;
  const Double_t *w = ALICE_2013_W;
  const Double_t *sig = ALICE_2013_sig;
  const Double_t *err = ALICE_2013_statsys;

  cout << "#####  Alice data for plot ######" << endl;
  cout.width(6);
  cout << "Wgp";
  cout.width(8);
  cout << "sigma" << endl;
  cout << fixed << setprecision(1);
  for(Int_t i=0; i<np; i++) {
    cout.width(6);
    cout << w[i];
    cout.width(8);
    cout << sig[i];
    cout << " +/- ";
    cout.width(4);
    cout << err[i];
    if( i == 0 ) cout << "   Forward p-Pb";
    else if( i == 3 ) cout << "   Semi-forward p-Pb";
    else if( i == 5 ) cout << "   Mid-rapidity";
    else if( i == 7 ) cout << "   Semi-forward Pb-p";
    else if( i == 8 ) cout << "   Forward Pb-p";
    cout << endl;
  }
  cout << "#####################################" << endl;

}//printData

//_____________________________________________________________________________
TGraphErrors* read_h1(const char* fileName){
  // read H1 data from http://www-h1.desy.de/psfiles/figures/d13-058.table_w_allH1_v1.txt
  ifstream f;
  f.open(fileName);
  string b;
  Double_t sigm[29];
  Double_t dsig[29];
  Double_t wavg[29];
  Double_t n,pd,wmin,wmax,PhiT;
  // header
  for (Int_t i=0;i<65;i++) getline(f,b);
  Int_t nExclusive=0;
  for (Int_t i=0;i<30;i++){
    f >> n >> pd;
    if (pd==1) { getline(f,b); continue; }
    f >> wmin >> wmax >> wavg[nExclusive] >> PhiT >> sigm[nExclusive] >> dsig[nExclusive];
    getline(f,b);
    nExclusive++;
  }
  return new TGraphErrors(nExclusive,wavg,sigm,NULL,dsig);
}//read_h1

//_____________________________________________________________________________
TGraphErrors* read_zeus_ee(){
  const int nZEUSee=14;
  Double_t wavgZEUSee[]    = { 27.5, 42.5, 55.0, 65.0, 75.0, 85.0, 100.0, 117.5, 132.5, 155.0, 185.0, 215.0, 245.0, 275.0 };
  Double_t sigmZEUSee[]    = { 33.6, 43.8, 57.2, 62.5, 68.9, 72.1, 81.9, 95.7, 103.9, 115.0, 129.1, 141.7, 140.3, 189.0 };
  Double_t statZEUSee[]    = { 1.6, 2.0, 1.8, 2.3, 2.6, 2.9, 2.3, 3.2, 3.6, 3.3, 4.7, 6.1, 7.4, 13.0 };
  Double_t systZEUSee[]    = { 2.4, 2.8, 3.5, 3.9, 4.5, 4.5, 4.8, 5.4, 5.8, 6.7, 7.7, 8.7, 9.9, 26.0 };
  Double_t dsigZEUSee[nZEUSee];
  for (Int_t i=0;i<nZEUSee;i++) dsigZEUSee[i]=sqrt(pow(statZEUSee[i],2.)+pow(systZEUSee[i],2.));
  return new TGraphErrors(nZEUSee,wavgZEUSee,sigmZEUSee,NULL,dsigZEUSee);
}//read_zeus_ee

//_____________________________________________________________________________
TGraphErrors* read_zeus_mm(){
  // ZEUS points
  const int nZEUSmm=8;
  Double_t wavgZEUSmm[]    = { 25.0, 40.0, 60.0, 80.0, 100.0, 120.0, 140.0, 160.0 };
  Double_t sigmZEUSmm[]    = { 32.6, 41.5, 55.8, 66.6,  73.4,  86.7, 104.0, 110.0 };
  Double_t statZEUSmm[]    = { 5.4, 1.1, 1.5, 2.0, 2.3, 3.2, 5.0, 11.0 };
  Double_t systZEUSmm[]    = { 5.2, 3.3, 4.6, 7.0, 6.0, 6.5,11.0, 12.0 };
  Double_t dsigZEUSmm[nZEUSmm];
  for (Int_t i=0;i<nZEUSmm;i++) dsigZEUSmm[i]=sqrt(pow(statZEUSmm[i],2.)+pow(systZEUSmm[i],2.));
  return new TGraphErrors(nZEUSmm,wavgZEUSmm,sigmZEUSmm,NULL,dsigZEUSmm);
}//read_zeus_mm

//_____________________________________________________________________________
TGraphErrors* read_lhcb(const char* fileName, const int nPoints, Int_t opt){
  ifstream f;
  f.open(fileName);
  Double_t w[nPoints],cs[nPoints],dcs[nPoints];
  for (Int_t i=0;i<nPoints;i++)  f >> w[i] >> cs[i] >> dcs[i];
  f.close();
  if (opt == 1)  return new TGraphErrors(10,w,cs,NULL,dcs);
  return new TGraphErrors(10,&(w[10]),&(cs[10]),NULL,&(dcs[10]));
}//read_lhcb

//_____________________________________________________________________________
TGraphErrors* read_jmrt(const char* fileName, const Int_t nPoints = 147){
  ifstream f;
  f.open(fileName);
  string b;
  getline(f,b);
  getline(f,b);
  Double_t q,w[nPoints],cs[nPoints],dcs[nPoints];
  for (Int_t i=0;i<nPoints;i++)  f >> q >> w[i] >> cs[i] >> dcs[i];
  f.close();
  return new TGraphErrors(nPoints,w,cs,NULL,dcs);
}//read_jmrt

//_____________________________________________________________________________
TGraphErrors* read_cgtt(const char *filenam, const Int_t nPoints = 41) {

  ifstream in;
  in.open(filenam);

  Double_t w[nPoints], cs[nPoints];
  //file loop
  for(Int_t i=0; i<nPoints; i++) in >> w[i] >> cs[i];
  in.close();

  return new TGraphErrors(nPoints, w, cs);

}//read_cgtt

//_____________________________________________________________________________
TGraph* read_mh_bfkl(const string filenam, const Int_t nPoints = 199) {

  Double_t wval, smin, smax;
  string com;

  ifstream f(filenam.c_str());
  for(Int_t i=0; i<4; i++) getline(f, com); // skip first 4 lines

  TGraph *g = new TGraph(2*nPoints);

  //file loop
  for(Int_t i=0; i<nPoints; i++) {
    f >> wval >> smin >> smax;
    g->SetPoint(i,   wval, smin);
    g->SetPoint(2*nPoints-i-1, wval, smax);
  }

  g->SetFillColor(kRed-9);
  g->SetLineColor(kRed-9);
  g->SetLineWidth(0);

  //fill style is in:
  //  https://root.cern.ch/doc/master/classTAttFill.html
  g->SetFillStyle(3207);

  f.close();

  return g;

}//read_mh_bfkl

//_____________________________________________________________________________
TGraph **read_amir_cgc(const string filenam, const Int_t nPoints = 196) {

  Double_t wval, sigma;
  string com;

  ifstream f(filenam.c_str());
  for(Int_t i=0; i<7; i++) getline(f, com); // skip first 7 lines

  TGraph *g = new TGraph(nPoints);
  TGraph *gup = new TGraph(nPoints/2);
  TGraph *glo = new TGraph(nPoints/2);

  //file loop
  for(Int_t i=0; i<nPoints; i++) {
    f >> wval >> sigma;
    g->SetPoint(i, wval, sigma);
    if(i < nPoints/2) gup->SetPoint(i, wval, sigma);
    if(i >= nPoints/2) glo->SetPoint(i-nPoints/2, wval, sigma);
  }//file loop

  f.close();

  g->SetFillColorAlpha(kOrange-2, 0.6);
  //g->SetFillColor(kOrange-2);
  g->SetLineColor(kOrange-2);
  g->SetLineWidth(1);

  g->SetFillStyle(3008);

  gup->SetLineColor(kOrange-2);
  gup->SetLineWidth(1);
  glo->SetLineColor(kOrange-2);
  glo->SetLineWidth(1);

  static TGraph *gcoll[3] = {g, glo, gup};

  return gcoll;

}//read_amir_cgc


//_____________________________________________________________________________
TGraph **read_amir_ipsat_bcgc(const string filenam, const Int_t nip = 98, const Int_t nb = 100) {

  Double_t wval, sigma;
  string com;

  ifstream f(filenam.c_str());
  for(Int_t i=0; i<8; i++) getline(f, com); // skip first 8 lines

  TGraph *gip = new TGraph(nip);
  TGraph *gb = new TGraph(nb);

  //IP-Sat loop
  for(Int_t i=0; i<nip; i++) {
    f >> wval >> sigma;
    gip->SetPoint(i, wval, sigma);
  }//IP-Sat loop

  for(Int_t i=0; i<7; i++) getline(f, com); // skip following 7 lines

  //b-CGC loop
  for(Int_t i=0; i<nb; i++) {
    f >> wval >> sigma;
    gb->SetPoint(i, wval, sigma);
  }//b-CGC loop

  f.close();

  static TGraph *gcoll[] = {gip, gb};

  return gcoll;

}//read_amir_ipsat_bcgc

//_____________________________________________________________________________
TGraph *getRatio(TGraph *g) {

  //ratio of model from TGraph argument over the data

  Double_t wval, sigma;
  Int_t i = 0;

  //output graph
  TGraph *gout = (TGraph*) g->Clone();

  //alice data
  alice *alic = alice::instance();

  //input graph loop
  while( g->GetPoint(i, wval, sigma) >=0 ) {

    gout->SetPoint(i, wval, sigma/alic->getPowerLaw(wval) );

    i++;
  }//input graph loop

  return gout;

}//getRatio
