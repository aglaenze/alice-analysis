
//global constants with range of the plot
const Double_t gxmin = 17., gymin = 0.3, gxmax = 1990., gymax = 1.7;

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
    static const Int_t nCorA = 3; // sources of correlated error
    Double_t Nval, Dval; // power-law parameters for further use
    TGraph *fitGraph; // to store fit result
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
    TGraphErrors *getForwardEightTeVDiss(bool stat, bool syst);
    void printData();
};//alice
alice alice::alic;

alice::alice(Int_t oo) : opt(oo), Nval(-1), Dval(-1), fitGraph(0) {
    
    
}//alice


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


// Analysis with Run 2 settings dissociative
TGraphErrors *alice::getForwardEightTeVDiss(bool stat, bool syst) {
    //forward data
    const Int_t np = 2;
    const Double_t w[np]   = { 32.8 /* 3.25 < y < 4. */, 47.7 /* 2.5 < y < 3.25 */ };
    const Double_t sig[np] = { 1.29, 1.21};
    const Double_t errStat[np] = {0.23, 0.18};
    const Double_t errSystV0[np] = {0.055, 0.052};
    const Double_t errSystSig[np] = {0.10, 0.09};
    
    Double_t errSyst[np] = {0, 0};
    Double_t err[np] = {0, 0};
    
    for (int i = 0; i < np; ++i) {
        errSyst[i] = sig[i] * TMath::Sqrt((errSystV0[i]/sig[i])*(errSystV0[i]/sig[i])
                                          + (errSystSig[i]/sig[i])*(errSystSig[i]/sig[i]));
        if (stat && syst) {
            err[i] = sig[i] * TMath::Sqrt( (errStat[i]/sig[i])*(errStat[i]/sig[i])
                                          + (errSyst[i]/sig[i])*(errSyst[i]/sig[i]));
        }
        else if (stat) {err[i] = errStat[i];}
        else if (syst) {err[i] = errSyst[i];}
    }
    
    TGraphErrors *gr = new TGraphErrors(np, w, sig, NULL, err);
    
    return gr;
    
}//getForward

//_____________________________________________________________________________
TGraphErrors* read_cct(){
    // read H1 data from https://hep.fjfi.cvut.cz/projects/hotspots/paper-1/Fig3_dis.txt
    const char* fileNameExc = "Fig3_exc.txt";
    const char* fileNameDis = "Fig3_dis.txt";
    // first get diss x section
    ifstream fExc;
    fExc.open(fileNameExc);
    string bExc;
    Double_t wavExc, sigExc;
    std::vector<Double_t> sigExcVec = {}, wExcVec = {};
    while (getline(fExc,bExc)) {
        //std::cout << b << endl;
        if (bExc.at(0) == '#') {
            //std::cout << b << endl;
            continue;
        }
        std::stringstream stream(bExc);
        stream >> wavExc >> sigExc;
        sigExcVec.push_back(sigExc);
        wExcVec.push_back(wavExc);
    }
    fExc.close();
    // second get exc x section
    ifstream fDis;
    fDis.open(fileNameDis);
    string bDis;
    Double_t wavDis, sigDis;
    std::vector<Double_t> sigDisVec = {}, wDisVec = {};
    while (getline(fDis,bDis)) {
        //std::cout << b << endl;
        if (bDis.at(0) == '#') {
            //std::cout << b << endl;
            continue;
        }
        std::stringstream stream(bDis);
        stream >> wavDis >> sigDis;
        sigDisVec.push_back(sigDis);
        wDisVec.push_back(wavDis);
    }
    fDis.close();
    
    const int size = sigExcVec.size();
    Double_t sigm[size];
    Double_t wavg[size];
    for (int i = 0; i < size; i ++) {
        sigm[i] = sigDisVec[i]/sigExcVec[i];
        wavg[i] = wExcVec[i];
    }
    return new TGraphErrors(size,wavg,sigm, NULL, NULL);
}//read_h1_mm


TGraphErrors* read_h1(){
    // H1
    const Int_t nH1DisW = 8;
    const Double_t WH1DisW[nH1DisW] = {43.5, 50, 57.3, 65.3, 73.9, 83.2, 93.3, 104.2};
    const Double_t sigH1DisW[nH1DisW] = {46.0, 52.1, 58.7, 58.7, 61.5, 67.7, 69.8, 68.9};
    const Double_t errH1DisW[nH1DisW] = {6.0, 6.5, 7.4, 7.5, 8.0, 8.7, 9.0, 9.0};
    const Double_t sigH1ExcW[nH1DisW] = {50.7, 59.5, 61.8, 67.6, 72.4, 79.9, 84.4, 86.7};
    const Double_t errH1ExcW[nH1DisW] = {4.9, 5.8, 6.2, 6.2, 6.4, 7.0, 7.0, 7.3};
    Double_t sigH1Ratio[nH1DisW] = {};
    Double_t errH1Ratio[nH1DisW] = {};
    for (int i = 0; i < nH1DisW; ++i) {
        sigH1Ratio[i] = sigH1DisW[i]/sigH1ExcW[i];
        errH1Ratio[i] = sigH1Ratio[i] * TMath::Sqrt((errH1DisW[i]/sigH1DisW[i])*(errH1DisW[i]/sigH1DisW[i]) + (errH1ExcW[i]/sigH1ExcW[i])*(errH1ExcW[i]/sigH1ExcW[i]) );
    }
    TGraphErrors *grH1DisW = new TGraphErrors(nH1DisW,WH1DisW,sigH1Ratio,NULL,errH1Ratio);
    grH1DisW->SetMarkerStyle(24);
    grH1DisW->SetMarkerColor(2);
    grH1DisW->SetLineColor(2);
    return grH1DisW;
}//read_h1


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

//_____________________________________________________________________________
TGraphErrors* read_jimwlk(const char* fileName){
    ifstream f;
    Int_t nPoints = 0;
    string b;
    f.open(fileName);
    getline(f,b);
    getline(f,b);
    getline(f,b);
    getline(f,b);
    while (getline(f,b)) {
        nPoints++;
    }
    f.close();
    f.open(fileName);
    getline(f,b);
    getline(f,b);
    getline(f,b);
    getline(f,b);
    Double_t w[nPoints],cs[nPoints],dcs[nPoints];
    for (Int_t i=0;i<nPoints;i++)  {
        f >> w[i] >> cs[i] >> dcs[i];
        //cout << "w  = " << w[i] << " cs = " << cs[i] << " dcs = " << dcs[i] << endl;
    }
    f.close();
    TGraphErrors* g = new TGraphErrors(nPoints,w,cs,NULL,dcs);
    //g->SetFillColor(kBlack);
    g->SetFillColorAlpha(kGray, 0.1);
    g->SetLineColor(kBlack);
    g->SetLineWidth(2);
    //fill style is in:
    //  https://root.cern.ch/doc/master/classTAttFill.html
    //g->SetFillStyle(3207);
    
    return g;
}//read_jimwlk

//_____________________________________________________________________________
TGraphErrors* read_jimwlk2(){
    const char* fileName = "cgc_coh_incoh_xs_wdep.txt";
    ifstream f;
    f.open(fileName);
    string b;
    Double_t wav, sigCoh, sigIncoh;
    std::vector<Double_t> sigVec = {}, wVec = {};
    while (getline(f,b)) {
        //std::cout << b << endl;
        if (b.at(0) == '#') {
            //std::cout << b << endl;
            continue;
        }
        std::stringstream stream(b);
        stream >> wav >> sigCoh >> sigIncoh;
        sigVec.push_back(sigIncoh/sigCoh);
        wVec.push_back(wav);
    }
    const int size = sigVec.size();
    Double_t sigm[size];
    Double_t wavg[size];
    for (int i = 0; i < size; i ++) {
        sigm[i] = sigVec[i];
        wavg[i] = wVec[i];
    }
    return new TGraphErrors(size,wavg,sigm, NULL, NULL);
}//read_jimwlk2


TGraph* GetShade(Double_t N, Double_t d, Double_t cov00, Double_t cov11, Double_t cov01)
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

