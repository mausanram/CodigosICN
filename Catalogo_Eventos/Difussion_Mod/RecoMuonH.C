#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include "TMinuit.h"

using namespace std;

//#include <TROOT.h>
//#include "Minuit2/Minuit2Minimizer.h"
//#include "Math/Functor.h"
// #include "TMinuit.h"

const float kD          = 725;
const int kMaxNPixels   = 50000;

//Parameters
#define dsigmad		0.01	//Resolution of the seep over sigmad.
#define SDTOP		0.7
#define SDBOTTOM	0.1
#define dmu		0.02	//Resolution of the sweep over mux and muy.
#define MUTop		0.5
#define MUBottom	-0.5

//Constants
#define RAIZ2PI 2.506628274631000
#define LOG_RAIZ2PI 0.918938533204673
#define RAIZ2   1.414213562373095
#define MAXHITPIXELS 1000
#define PI 3.14159265359

typedef struct {
    double mux;
    double charge;
    double sigmad;
    double maxEventPDF;
} event1D_param; //Event parameters from ML algorithm.

//Global Variables
double ampnoise;		//Noise of the amplifier in ADU.
double ampnoiseVar;	//Variance of the amplifier noise.
vector<double> x_pix;
vector<double> q_pix;
uint npix_event=0;
long double *cumlog; //cumlog vector have the result of ln(n!)
event1D_param parameters;
double QTii;
double gain = 0;
int flag=0;

//******************
//Return the ln(n!)
//******************
int loadcumlog(uint qtmax)
{
    uint i;
    cumlog = (long double *) malloc(qtmax*sizeof(long double));
    cumlog[0]=0;
    for(i=1 ; i<qtmax ; i++) {
        cumlog[i]= cumlog[i-1] + log(i);
    }
    return 0;
}

// ***********
// Apply de ML
// ***********
//double applyml1D(const double *xx)
void fcn(int& npar, double* deriv, double& f, double par[], int flag)
{

    const double mux = par[0];//valor de la media en columns a evaluar
    const double sd = par[1]; //valor de la difusion a evaluar
    // const int QT = (int)floor(par[2]);
    const double QT = par[2];

    //Variables
		double xmax[kMaxNPixels];
    double xmin[kMaxNPixels];
    double lambda[kMaxNPixels];
    double log_lambda[kMaxNPixels];
    double log_oneMinuslambda[kMaxNPixels];
    double pdf_only_noise[kMaxNPixels];
    double q_pix_round[kMaxNPixels];
    long iimin[kMaxNPixels];
    
		double eventPDFConstant = npix_event*log(1/(ampnoise*RAIZ2PI));
    long iimax;
    long dq = (long) ceil(3*ampnoise);
    uint ii;	// Used during loops.
    long double eventPDF;
    long double exponente;
    long double pixelPDF;

    //initialize q_pix_round, iimin, pdf_only_noise
    uint pix;
    for (pix=0; pix<npix_event; pix++) {
        q_pix_round[pix] = round(q_pix[pix]);
        iimin[pix] = ((q_pix_round[pix]-dq)>0)? (q_pix_round[pix]-dq) : 0;
        pdf_only_noise[pix] = exp(-(q_pix[pix]*q_pix[pix])/(2*ampnoiseVar));
    }

    //calculate lambda for all pixels
    for (pix=0; pix<npix_event; pix++) {
        xmin[pix]=x_pix[pix]-0.5-mux;
        xmax[pix]=x_pix[pix]+0.5-mux;
        lambda[pix]=fabs(0.5*(erf(xmax[pix]/(RAIZ2*sd))-erf(xmin[pix]/(RAIZ2*sd))));
        log_lambda[pix] = log(lambda[pix]);
        log_oneMinuslambda[pix] = log(1.0-lambda[pix]);
    }

    eventPDF=0;
    for(pix=0; pix<npix_event; pix++) { //Sweep over event pixels
        // PIXEL PDF
        pixelPDF = 0;

        double safe_lambda = lambda[pix];
        if (safe_lambda < 1e-9) safe_lambda = 1e-9;
        if (safe_lambda > 1.0 - 1e-9) safe_lambda = 1.0 - 1e-9;

        double logL = log(safe_lambda);
        double logOneMinusL = log(1.0 - safe_lambda);

        // iimax = ((q_pix_round[pix]+dq) < QT) ? (q_pix_round[pix]+dq) : QT;
        // int iimax = ((q_pix_round[pix]+dq) < QT) ? (q_pix_round[pix]+dq) : (int)floor(QT);

        iimax = (int)floor(QT); // Use floor for the loop limit
        if ((q_pix_round[pix] + dq) < iimax) iimax = (int)floor(q_pix_round[pix] + dq);


        if (lambda[pix]<0.0001) {  //si lambda=0, entonces no hay contribucion de carga del evento y solo habra ruido gaussiano
            pixelPDF = pdf_only_noise[pix];
        } else {
            for( ii=iimin[pix] ; ii<=iimax ; ii++ ) {
                // exponente=cumlog[QT]+(ii*log_lambda[pix])+(QT-ii)*log_oneMinuslambda[pix]-pow(q_pix[pix]-ii,2)/(2*ampnoiseVar)-cumlog[ii]-cumlog[QT-ii];
                // pixelPDF=pixelPDF + exp(exponente);

                double arg3 = QT - (double)ii + 1.0;
                if (arg3 < 1e-9) arg3 = 1e-9;

                exponente=TMath::LnGamma(QT + 1.0)+ (ii * logL)+(QT-ii)*logOneMinusL
                         -pow(q_pix[pix]-ii,2)/(2*ampnoiseVar)-TMath::LnGamma((double)ii + 1.0)-TMath::LnGamma(arg3);

                pixelPDF+= exp(exponente);
            }
        }
        //EVENT PDF
        if (pixelPDF < 1e-30) pixelPDF = 1e-30;
        eventPDF += log(pixelPDF);
    }
    eventPDF += eventPDFConstant;

    f=-1.0*eventPDF;
		return;
}

//*************************************************
//Funcion para ajustar el sigma del slice del muon
//*************************************************
double fitSpread1D(double qInit, double muxInit, double sdInit)
{
    //Initialize minuit
    //http://www.desy.de/~rosem/flc_statistics/data/likelihood_example.cpp
    double arglist[10];
    int ierflg;
    double fmin,fedm;
    double errdef;
    int nparv,nparx;
    int fstat;
    TString pname;
    double pvalue, perror;
    double plbound, pubound;
    int pvari;
    const int npar=3;			//number of FCN parameters (sigma,mux,muy,Q)
    TMinuit minuit(npar);
    minuit.SetFCN(fcn); //Function to be minimized
    // minuit.SetPrintLevel(2); // Enable to see the fit status
    minuit.SetPrintLevel(-1);

    double step[4] = {0.01,0.01,1.0};
    double muXLowLim = muxInit - 0.5;
    double muXUpLim = muxInit + 0.5;
    double sdLowLim = 0.0;
    double sdUpLim = 1.7;
    double dQT = 3*sqrt(npix_event*ampnoiseVar);
    double QTmin = (floor(qInit-dQT)>=1)? floor(qInit-dQT) : 1;
    double QTmax = ceil(qInit+dQT);

    minuit.mnparm(0,"y",muxInit,step[0],muXLowLim,muXUpLim,ierflg);
    minuit.mnparm(1,"sd",sdInit,step[1],sdLowLim,sdUpLim,ierflg);
    minuit.mnparm(2,"q",qInit,step[2],QTmin,QTmax,ierflg);

    //******Monte Carlo
    arglist[0]=5000;
    // arglist[1] = 3;
    minuit.mnexcm("SEEK",arglist,1,ierflg); //execute minimization.

    if(ierflg==0) {
        minuit.mnstat(fmin,fedm,errdef,nparv,nparx,fstat); //read current status of the minimization.
        parameters.maxEventPDF = fmin/npix_event;
        minuit.mnpout(0,pname,pvalue,perror,plbound,pubound,pvari);
        parameters.mux = pvalue;
        minuit.mnpout(1,pname,pvalue,perror,plbound,pubound,pvari);
        parameters.sigmad = pvalue;
        minuit.mnpout(2,pname,pvalue,perror,plbound,pubound,pvari);
        parameters.charge = pvalue;
    } else {
        double errorCode = 4000000;
        parameters.mux = errorCode;
        parameters.sigmad = errorCode;
        parameters.maxEventPDF = errorCode;
    }

    return parameters.sigmad;
}

//Mass:
double bsq()
{
    double bsTotq=0;
    for(uint i=0; i<q_pix.size(); i++) {
        bsTotq=bsTotq+q_pix[i];
    }
    return bsTotq;
}

//Center of mass:
double bsmeanx(double totq)
{
    double meany=0;
    for(uint i=0; i<x_pix.size(); i++) {
        meany=meany+(x_pix[i]*q_pix[i]);
    }
    return meany/totq;
}

//Moment of inertia:
double bssigma(double totq, double meanx)
{
    double sigma=0;
    for(uint i=0; i<x_pix.size(); i++) {
        sigma=sigma+((x_pix[i]-meanx)*(x_pix[i]-meanx)*q_pix[i]);
    }
    return sqrt(sigma/totq);
}


void RecoMuonH(const char *muonfile, const char *ofilename)
{
    // std::cout<<"I'm in recoDiff"<<endl;
    //Create file with the muon event:
    std::ifstream mfile(muonfile);

    if (!mfile.is_open()) {
        std::cerr << "Error opening file!" << std::endl;
        // std::exit();
    }

    string line;
    vector<double> x, y, e;
    double x_val, y_val, e_val;
    double ee, eemax=0;

    while (getline(mfile, line)) {
        // std::cout<< line <<endl;
        if (line.empty() || line[0] == '#') continue;

        stringstream ss(line);
        if (ss >> x_val >> y_val >> e_val) {
            x.push_back(x_val);
            y.push_back(y_val);
            e.push_back(e_val);
            if(ee>eemax){eemax=ee;}
        }
    }
    // std::cout<<"xye file read succesfuly!"<<endl;
    mfile.close();

    // ampnoise = 0.00109; // EXT1 keV
    // ampnoise = 0.000907; // EXT2 keV

    ampnoise = 0.298; // EXT1 e-
    // ampnoise = 0.247; // EXT2 e-

    ampnoiseVar = pow(ampnoise,2);

    // std::cout<< "I'm ampnoisevar: "<< ampnoiseVar<< endl;
    
    loadcumlog((uint)round(eemax+1000));    //Load cumlog vector. cumlog vector have the result of ln(n!).
    
    auto max_it_x = std::max_element(std::cbegin(x), std::cend(x));
    double max_value_x = *max_it_x;

    auto min_it_x = std::min_element(std::cbegin(x), std::cend(x));
    double min_value_x = *min_it_x;

    auto max_it_y = std::max_element(std::cbegin(y), std::cend(y));
    double max_value_y = *max_it_y;

    auto min_it_y = std::min_element(std::cbegin(y), std::cend(y));
    double min_value_y = *min_it_y;

    double y_size = max_value_y - min_value_y;

    // std::cout << "The X maximum value is: " << max_value_x << std::endl;
    // std::cout << "The X minimum value is: " << min_value_x << std::endl;
    // std::cout << "The Y maximum value is: " << max_value_y << std::endl;
    // std::cout << "The Y minimum value is: " << min_value_y << std::endl;
    // std::cout<< "Y size: " << y_size << std::endl;

    double dl_depth = kD/y_size;

    //Lateral spread of each muon slice
    double fitSigma;
    vector<double> recSigma;
    vector<double> SigmaBS; // Vector with the sigma estimate from basic statistics.
    vector<double> sliceQ; // This array will contain the charge of each slice.
    vector<double> sliceY;
    for(Int_t i=min_value_x; i<max_value_x+1; i++) { //loop over the slices
        //Get One Slice
        x_pix.clear(); //This is a global variable used during fitSpread1D
        q_pix.clear();
        for(uint j=0; j<x.size(); ++j) {
            if(x[j]==i) {
                x_pix.push_back(y[j]);
                q_pix.push_back(e[j]);
	        }
        }
	    npix_event=x_pix.size();

        //Basic Statistics of the slice
        double bsTotq=bsq();
        double bsMeanx=bsmeanx(bsTotq);
        double bsSigma=bssigma(bsTotq,bsMeanx);
        SigmaBS.push_back(bsSigma);

        // Fit spread
        fitSigma=fitSpread1D(bsTotq,bsMeanx,bsSigma); 
        // std::cout<< "I'm ampnoisevar: "<< ampnoiseVar<< endl;
        // std::cout<<bsSigma << " " <<fitSigma<<endl;
        recSigma.push_back(fitSigma); // We use this value for the depth-spread plot
        sliceQ.push_back(bsTotq); // Total charge of the line
        sliceY.push_back(i); // N-Line, NOT depth
        // sliceY.push_back((i - min_value_y)*dl_depth); // Depth DON'T WORK WELL!!!
    }

    // Muon data for fit
    ofstream ofile(ofilename);
    for(uint i=0; i<sliceQ.size(); i++) {
        // ofile << sliceQ[i] << " " << recSigma[i] << " " << SigmaBS[i]<< " " << sliceY[i] << endl;
        ofile << recSigma[i] << " " << sliceY[i] << endl;
    }
    ofile.close();
}

int RecoMuons(int argc, char argv[])
{
    // The first argument to the program is the xye file
    // The second argument is the muonID in the catalog.
    // This program process only one muon at a time.
    // argv[0] program name
    // argv[1] output xye file
    // argv[2] output txt file

    std::cout<<argv<< endl;

    // flag=0;
    // if(argc==7) {
    //   gain= argv[5];
    //   ampnoise= argv[6];
    //   flag=1;
    // } 

    std::cout<<argv[1]<< endl;
    // RecoMuon2(argv[1]);

    free(cumlog);
    return 0;
}
