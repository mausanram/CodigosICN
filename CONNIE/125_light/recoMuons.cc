#include <sstream>
#include <string>
#include <fstream>
#include <vector>
#include <cmath>
#include <iostream>
using namespace std;
#include "TFile.h"
#include "TNtuple.h"
//#include <TROOT.h>
//#include "Minuit2/Minuit2Minimizer.h"
//#include "Math/Functor.h"
#include "TMinuit.h"

const float kD          = 675;
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

    const Double_t mux = par[0];//valor de la media en columns a evaluar
    const Double_t sd = par[1]; //valor de la difusion a evaluar
    const int QT = (int)floor(par[2]);

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
        iimax = ((q_pix_round[pix]+dq) < QT)? (q_pix_round[pix]+dq) : QT;
        if (lambda[pix]<0.0001) {  //si lambda=0, entonces no hay contribucion de carga del evento y solo habra ruido gaussiano
            pixelPDF = pdf_only_noise[pix];
        } else {
            for( ii=iimin[pix] ; ii<=iimax ; ii++ ) {
                exponente=cumlog[QT]+(ii*log_lambda[pix])+(QT-ii)*log_oneMinuslambda[pix]-pow(q_pix[pix]-ii,2)/(2*ampnoiseVar)-cumlog[ii]-cumlog[QT-ii];
                pixelPDF=pixelPDF + exp(exponente);
            }
        }
        //EVENT PDF
        eventPDF = eventPDF + log(pixelPDF);
    }
    eventPDF = eventPDF + eventPDFConstant;

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
    minuit.SetPrintLevel(-1);

    double step[4] = {0.01,0.01,1.0};
    double muXLowLim = muxInit - 0.5;
    double muXUpLim = muxInit + 0.5;
    double sdLowLim = 0.0;
    double sdUpLim = 1.5;
    double dQT = 3*sqrt(npix_event*ampnoiseVar);
    double QTmin = (floor(qInit-dQT)>=1)? floor(qInit-dQT) : 1;
    double QTmax = ceil(qInit+dQT);

    minuit.mnparm(0,"x",muxInit,step[0],muXLowLim,muXUpLim,ierflg);
    minuit.mnparm(1,"sd",sdInit,step[1],sdLowLim,sdUpLim,ierflg);
    minuit.mnparm(2,"q",qInit,step[2],QTmin,QTmax,ierflg);

    //******Monte Carlo
    arglist[0]=5000;
    minuit.mnexcm("SEE",arglist,1,ierflg); //execute minimization.

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
    double meanx=0;
    for(uint i=0; i<x_pix.size(); i++) {
        meanx=meanx+(x_pix[i]*q_pix[i]);
    }
    return meanx/totq;
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

void recoDiff(TTree *data, const int eventID, const char muonfile[],const char ofilename[])
{

    data->SetBranchStatus("*",0);
    
    Int_t n;
    data->SetBranchStatus("nSavedPix", 1);
    data->SetBranchAddress("nSavedPix", &n);

    Int_t xPix[kMaxNPixels];
    data->SetBranchStatus("xPix", 1);
    data->SetBranchAddress("xPix", &xPix);

    Int_t yPix[kMaxNPixels];
    data->SetBranchStatus("yPix", 1);
    data->SetBranchAddress("yPix", &yPix);

    Float_t ePix[kMaxNPixels];
    data->SetBranchStatus("ePix", 1);
    data->SetBranchAddress("ePix", &ePix);

    Float_t level[kMaxNPixels];
    data->SetBranchStatus("level", 1);
    data->SetBranchAddress("level", &level);

    Int_t runID;
    data->SetBranchStatus("runID", 1);
    data->SetBranchAddress("runID", &runID);

    Int_t ohdu;
    data->SetBranchStatus("ohdu", 1);
    data->SetBranchAddress("ohdu", &ohdu);

    Int_t yMin;
    data->SetBranchStatus("yMin", 1);
    data->SetBranchAddress("yMin", &yMin);

    Int_t yMax;
    data->SetBranchStatus("yMax", 1);
    data->SetBranchAddress("yMax", &yMax);

    Int_t error;

    //Float_t gainCu;
    //data->SetBranchStatus("gainCu", 1);
    //error = data->SetBranchAddress("gainCu", &gainCu);
    //
    //Float_t osMadSCN;
    //data->SetBranchStatus("osMadSCN", 1);
    //error = data->SetBranchAddress("osMadSCN", &osMadSCN);

    Double_t gainCu;
    data->SetBranchStatus("GainI", 1);
    error = data->SetBranchAddress("GainI", &gainCu);
    
    //Double_t osMadSCN;
    Float_t osMadSCN;
    data->SetBranchStatus("Noise", 1);
    error = data->SetBranchAddress("Noise", &osMadSCN);

    data->GetEntry(eventID); //**** Here the muon is read *****

    if(error==0) {
      //gain= (double) gainCu/(1000/3.745); //Gain in ADU/e.
      //ampnoise=((double) osMadSCN*1.48)/gain;     // Value of amp noise
      gain= 1.0;                   //AA: skp images are calibrated [pixval in in e-]
      ampnoise= (double) osMadSCN; //AA: skp images noise is in header and given in e-
    } else {
      cout << endl << "No gain and noise in the input ROOT file." << endl;
      if(flag==0){
        cout<<endl<<"Gain and noise were not provided as input arguments. Bye."<<endl;
        exit(1);
      }
    }

    //cout<<endl<<"Gain= "<<gain<<"adu/e-"<<endl;
    //cout<<endl<<"Noise= "<<ampnoise<<"e-"<<endl;
    // AA: skp images are calibrated (pixval in e-) and Noise is in header, given in e-
    cout<<endl<<"Gain= "<<gain<<" (images are calibrated)"<<endl;
    cout<<endl<<"Noise= "<<ampnoise<<"e-"<<endl;
		
    ampnoiseVar=ampnoise*ampnoise;	//Variance of the noise.

    //Create file with the muon event:
    ofstream mfile;
    mfile.open(muonfile);
    mfile<<"# runID: "<<runID<<endl;
    mfile<<"# ohdu: "<<ohdu<<endl;
    mfile<<"# Gain: "<<gain<<endl;
    mfile<<"# Noise: "<<ampnoise<<endl;
    vector<double> x;
    vector<double> y;
    vector<double> e;
    double ee, eemax=0;
    for(Int_t i=0; i<n; ++i) {
        if(level[i]==0) {
	    ee=(double) ePix[i]/gain;
            x.push_back((double) xPix[i]);
            y.push_back((double) yPix[i]);
            e.push_back(ee);
            mfile << xPix[i] << " " << yPix[i] <<" "<< ee <<endl;
	    if(ee>eemax){eemax=ee;}
        }
    }
    mfile.close();
    
    loadcumlog((uint)round(eemax+1000));    //Load cumlog vector. cumlog vector have de result of ln(n!).

    //Lateral spread of each muon slice
    double fitSigma;
    vector<double> recSigma;
    vector<double> SigmaBS; // Vector with the sigma estimate from basic statistics.
    vector<double> sliceQ; // This array will contain the charge of each slice.
    vector<double> sliceY;
    for(Int_t i=yMin; i<yMax+1; i++) { //loop over the slices
        //Get One Slice
        x_pix.clear(); //This is a global variable used during fitSpread1D
        q_pix.clear();
        for(uint j=0; j<y.size(); ++j) {
            if(y[j]==i) {
                x_pix.push_back(x[j]);
                q_pix.push_back(e[j]);
	    }
        }
	npix_event=x_pix.size();
        //Basic Statistics of the slice
        double bsTotq=bsq();
        double bsMeanx=bsmeanx(bsTotq);
        double bsSigma=bssigma(bsTotq,bsMeanx);
        SigmaBS.push_back(bsSigma);
        //Fit spread
	fitSigma=fitSpread1D(bsTotq,bsMeanx,bsSigma); 
	recSigma.push_back(fitSigma);
        sliceQ.push_back(bsTotq);
        sliceY.push_back(i);
    }
    //Muon data for fit
    ofstream ofile;
    ofile.open(ofilename);
    for(uint i=0; i<sliceQ.size(); i++) {
        ofile << sliceQ[i] << " " << recSigma[i] << " " << SigmaBS[i]<< " " << sliceY[i] << endl;
    }
    ofile.close();
}

int main(int argc, char const *argv[])
{
    // The first argument to the program is the catalog.root with the muons.
    // The second argument is the muonID in the catalog.
    // This program process only one muon at a time.
    // argv[0] program name
    // argv[1] ROOT catalog
    // argv[2] muon ID
    // argv[3] output xye file
    // argv[4] output txt file
    // argv[5] gain
    // argv[6] noise

    TFile *file0 = new TFile(argv[1],"READ"); //Read ROOT catalog.
    TNtuple *data = (TNtuple*)file0->Get("hitSumm");

    flag=0;
    if(argc==7) {
      gain= atof(argv[5]);
      ampnoise= atof(argv[6]);
      flag=1;
    } 

    recoDiff(data, atoi(argv[2]),argv[3],argv[4]);

    free(cumlog);
    return 0;
}
