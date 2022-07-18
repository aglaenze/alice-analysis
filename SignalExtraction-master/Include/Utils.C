#include <iostream>
#include <fstream>
#include <dirent.h>
#include <string>
#include <cmath>


using namespace std;

bool Initiate(string period, Double_t& mMin, Double_t& mMax, Double_t& ptMin, Double_t& ptMax, Double_t& rapMin, Double_t& rapMax, bool& useCuts, bool& logScale, bool& drawPulls, bool& exp, bool& exclusiveOnly) {
	ifstream file("input-" + period + ".txt", ios::in);
	string a, b;
	bool mMinFound = false, mMaxFound = false, ptMinFound = false, ptMaxFound = false, rapMinFound = false, rapMaxFound = false, useCutsFound = false, logScaleFound = false, drawPullsFound = false, expFound = false, exclusiveOnlyFound = false;
	if (file) {
		string line;
		getline(file, line);	//first line does not contains info
		cout << "#################################" << endl;
		while(getline(file, line)) {
			if (line.find("mMin") != string::npos) {
				stringstream stream(line);
				stream >> a >> b >> mMin;
				mMinFound = true;
				cout << "#\tmMin = " << mMin << "\t\t#" << endl;
			}
			else if (line.find("mMax") != string::npos) {
				stringstream stream(line);
				stream >> a >> b >> mMax;
				mMaxFound = true;
				cout << "#\tmMax = " << mMax << "\t\t#" << endl;
			}
			else if (line.find("ptMin") != string::npos) {
				stringstream stream(line);
				stream >> a >> b >> ptMin;
				ptMinFound = true;
				cout << "#\tptMin = " << ptMin << "\t\t#" << endl;
			}
			else if (line.find("ptMax") != string::npos) {
				stringstream stream(line);
				stream >> a >> b >> ptMax;
				ptMaxFound = true;
				cout << "#\tptMax = " << ptMax << "\t\t#" << endl;
			}
			else if (line.find("rapMin") != string::npos) {
				stringstream stream(line);
				stream >> a >> b >> rapMin;
				rapMinFound = true;
				cout << "#\trapMin = " << rapMin << "\t\t#" << endl;
			}
			else if (line.find("rapMax") != string::npos) {
				stringstream stream(line);
				stream >> a >> b >> rapMax;
				rapMaxFound = true;
				cout << "#\trapMax = " << rapMax << "\t\t#" << endl;
			}
			else if (line.find("useCuts") != string::npos) {
				stringstream stream(line);
				stream >> a >> b >> useCuts;
				useCutsFound = true;
				cout << "#\tuseCuts = " << useCuts << "\t\t#" << endl;
			}
			else if (line.find("logScale") != string::npos) {
				stringstream stream(line);
				stream >> a >> b >> logScale;
				logScaleFound = true;
				cout << "#\tlogScale = " << logScale << "\t\t#" << endl;
			}
			else if (line.find("drawPulls") != string::npos) {
				stringstream stream(line);
				stream >> a >> b >> drawPulls;
				drawPullsFound = true;
				cout << "#\tdrawPulls = " << drawPulls << "\t\t#" << endl;
			}
			else if (line.find("exp") != string::npos) {
				stringstream stream(line);
				stream >> a >> b >> exp;
				expFound = true;
				cout << "#\texp = " << exp << "\t\t\t#" << endl;
			}
			else if (line.find("exclusiveOnly") != string::npos) {
				stringstream stream(line);
				stream >> a >> b >> exclusiveOnly;
				exclusiveOnlyFound = true;
				cout << "#\texclusiveOnly = " << exclusiveOnly << "\t#" << endl;
			}
		}
		cout << "#################################" << endl;
		cout << endl;
		return (mMinFound && mMaxFound && ptMinFound && ptMaxFound && rapMinFound && rapMaxFound && useCutsFound && logScaleFound && drawPullsFound && expFound && exclusiveOnlyFound);
	}
	else cout << "Error: not possible to open input.txt file in reading mode" << endl;
	return false;
}


bool GetParameters(string period, Double_t& bExc, Double_t& gammaPbYield, Double_t& muLandau, Double_t& sigmaLandau, Double_t& pt0, Double_t& nInc) {
	ifstream file("input-" + period + ".txt", ios::in);
	string a, b;
	bool bExcFound = false, gammaPbYieldFound = false, muFound = false, sigmaFound = false, pt0Found = false, nIncFound = false;
	if (file) {
		string line;
		getline(file, line);	//first line does not contains info
		while(getline(file, line)) {
			if (line.find("bExc") != string::npos) {
				stringstream stream(line);
				stream >> a >> b >> bExc;
				//cout << "\n\n\n\n\n\nb = " << bExc << endl << endl << endl << endl << endl;
				bExcFound = true;
			}
			else if (line.find("gammaPbYield") != string::npos) {
				stringstream stream(line);
				stream >> a >> b >> gammaPbYield;
				gammaPbYieldFound = true;
			}
            else if (line.find("muLandau") != string::npos) {
                stringstream stream(line);
                stream >> a >> b >> muLandau;
                muFound = true;
                cout << "#\tmu(Landau) = " << muLandau << "\t\t#" << endl;
            }
            else if (line.find("sigmaLandau") != string::npos) {
                stringstream stream(line);
                stream >> a >> b >> sigmaLandau;
                sigmaFound = true;
                cout << "#\tsigma(Landau) = " << sigmaLandau << "\t\t#" << endl;
            }
            else if (line.find("pt0") != string::npos) {
                stringstream stream(line);
                stream >> a >> b >> pt0;
                pt0Found = true;
                cout << "#\tpt0 = " << pt0 << "\t\t#" << endl;
            }
            else if (line.find("nInc") != string::npos) {
                stringstream stream(line);
                stream >> a >> b >> nInc;
                nIncFound = true;
                cout << "#\tnInc = " << nInc << "\t\t#" << endl;
            }
		}
		return (bExcFound && gammaPbYieldFound && muFound && sigmaFound && pt0Found && nIncFound);
	}
	else cout << "Error: not possible to open input.txt file in reading mode" << endl;
	return false;
}

bool GetParameters(string period, Double_t& bExc, Double_t& gammaPbYield) {
    Double_t muLandau = 0, sigmaLandau = 0;
    Double_t pt0 = 0, nInc = 0;
    return GetParameters(period, bExc, gammaPbYield, muLandau, sigmaLandau, pt0, nInc);
}

bool GetParametersLandau(string period, Double_t& muLandau, Double_t& sigmaLandau, Double_t& pt0, Double_t& nInc) {
    Double_t bExc = 0, gammaPbYield = 0;
    return GetParameters(period, bExc, gammaPbYield, muLandau, sigmaLandau, pt0, nInc);
}
