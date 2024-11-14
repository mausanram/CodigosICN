Double_t comp_poisson(Int_t k, Double_t lambda, Double_t pgeom) {

    //Double_t k = int(x[0]);
    //Double_t lambda = par[0];
    //Double_t pgeom  = par[1];

    const int KMAX = 1000;
    Double_t P[KMAX];
    Double_t q = 1.-pgeom;

    // Double_t q = 1.;

    Double_t z=(lambda*q/(1.0-q));
    //printf("Valor de z=%f\n",z);

    P[0]=TMath::Exp(-lambda);
    P[1]=P[0]*(1.0-q)*z;

    if (k>1) {
       P[k]=P[k-1] *( (1.-q) * ( (2.*k) - 2. + z )/k ) +  P[k-2] * ((2.-k)/k) * (1.-q) * (1.-q) ;
       //for (int i=2;i<k;i++){
       //  double xx = (i);
       //  P[i]=P[i-1] *( (1.-q) * ( (2.*xx) - 2. + z )/xx ) +  P[i-2] * ((2.-xx)/xx) * (1.-q) * (1.-q) ;
       //} //for i
    } // if k

   return P[k];

}//comp_poisson

Double_t gauss_comppoisson_fit(Double_t *x, Double_t *par) {

    Int_t k =1000;
    //Int_t m = 4;
    //Double_t ydata = 0;
    Double_t xval = x[0];
    Double_t a     = par[0];
    Double_t mu    = par[1];
    Double_t sigma = par[2];    
    Double_t lambda_poisson = par[3];
    Double_t pgeom = par[4];
    Double_t gain  = par[5];
    // Double_t b = par[4];
    // Double_t c = par[5];
    // Double_t d = par[6];
    // Double_t e = par[7];
    // Double_t f = par[8];
    // Double_t g = par[9];
    // Double_t h = par[10];
    // Double_t i = par[11];
    // Double_t j = par[12];
    // Double_t k2 = par[13];
    
    //Double_t comppoisson(Int_t, Double_t, Double_t);
    

    Double_t fitval = 0.0;
    for (Int_t p = 0; p <= k; p++){
        // fitval += a * TMath::Gaus(xval*gain,p-mu,sigma,1) * TMath::PoissonI(p,lambda_poisson);
          fitval += a * TMath::Gaus(xval*gain,p-mu,sigma,1) * comp_poisson(1.*p,lambda_poisson, pgeom);
               //+ b * TMath::Gaus(xval,mu+1,sigma,1) * TMath::PoissonI(p,1)
               //+ c * TMath::Gaus(xval,mu+2,sigma,1) * TMath::PoissonI(p,2)
               //+ d * TMath::Gaus(x,mu+3,sigma,1) * TMath::PoissonI(p,3)
               //+ e * TMath::Gaus(x,mu+4,sigma,1) * TMath::PoissonI(p,4)
               //+ f * TMath::Gaus(x,mu+5,sigma,1) * TMath::PoissonI(p,5)
               //+ g * TMath::Gaus(x,mu+6,sigma,1) * TMath::PoissonI(p,6)
               //+ h * TMath::Gaus(x,mu+7,sigma,1) * TMath::PoissonI(p,7)
               //+ i * TMath::Gaus(x,mu+8,sigma,1) * TMath::PoissonI(p,8)
               //+ j * TMath::Gaus(x,mu+9,sigma,1) * TMath::PoissonI(p,9)
               //+ k2 * TMath::Gaus(x,mu+10,sigma,1) * TMath::PoissonI(p,10);

    }//for Int_t p

    return fitval;

}//gauss_comppoisson_fit


//----------------------------------------------------
// ***************************************************


void fit_gausCONVcomppois()
{//begin

  //------------- Style --------------
  //gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  //------------- Style --------------

  //----------  Text Style  ---------
  //ft = 10 * fontID + precision
  Int_t ft = 10 * 4 + 2;
  Double_t sz = 0.04;
  //---------------------------------


//----------------------------------------------------
//- csv file description
//----------------------------------------------------

//- File has 1011 rows with two doubles separated by a comma ','
//- Each row is intepreted as the lower bin edge and bin content.
//- Last row is interpreted as the upper edge of the last bin and 
//- the overflow bin content.

//----------------------------------------------------
//- Read csv file into Histogram
//----------------------------------------------------
//const int NBinsPlus1 = 1011;
//int       NBins      = 1010;
const int NBinsPlus1 = 990;
int       NBins      = 989;
double    XBins[NBinsPlus1]; 
double    BinCont[NBinsPlus1]; 

ifstream inFile;
//inFile.open("histograma_ext1.csv");
inFile.open("histograma_sin_SRE.csv");
double val1, val2;
char del;
int ibin = 0;
while ( (inFile >> val1 >> del >> val2) && (del==',') )
{
  //printf("%4d: %2.16f ,%2.19f\n",ibin,val1,val2);
  XBins[ibin]   = val1;
  BinCont[ibin] = val2;
  ibin++;
} //while inFile

TH1F  *histo = new TH1F("histo","",NBins, XBins);
histo->GetXaxis()->SetTitle("charge (e^{-})");
for (int ibin = 0; ibin<NBins; ibin++){
   histo->SetBinContent(ibin,BinCont[ibin]);
}
//double sf = 1/0.0003232822194619644;
//histo->Scale(sf);
double nevents = histo->Integral();
cout <<"events: "<< nevents  << endl;

//----------------------------------------------------
//- Fit function definition
//----------------------------------------------------

int npar = 6;
double norm = 3091.7;
double offs = 0.020;
double sigm = 0.211;
double lamb = 0.163;
double pgeo = 0.083;
double gain = 1.014;

double lofit = -0.5;
double hifit =  4.5;

cout <<"HERE 1!"<<endl;

TF1 *fitf = new TF1("fitf",gauss_comppoisson_fit,lofit,hifit,npar);
fitf->SetParameter(0,norm);
fitf->SetParameter(1,offs);
fitf->SetParameter(2,sigm);
fitf->SetParameter(3,lamb);
fitf->SetParameter(4,pgeo);
fitf->SetParameter(5,gain);
fitf->SetNpx(400);
fitf->SetMinimum(1e-3);
fitf->SetLineWidth(1);

cout <<"HERE 2!"<<endl;

//----------------------------------------------------
//- Fit to histogram
//----------------------------------------------------

histo->Fit("fitf","R");
norm = fitf->GetParameter(0);
offs = fitf->GetParameter(1);
sigm = fitf->GetParameter(2);
lamb = fitf->GetParameter(3);
pgeo = fitf->GetParameter(4);
gain = fitf->GetParameter(5);

double norme = fitf->GetParError(0);
double offse = fitf->GetParError(1);
double sigme = fitf->GetParError(2);
double lambe = fitf->GetParError(3);
double pgeoe = fitf->GetParError(4);
double gaine = fitf->GetParError(5);

double chisq = fitf->GetChisquare();
int    ndegf = fitf->GetNDF();
double proba = fitf->GetProb();

TF1 *fitfcp = new TF1("fitfcp",gauss_comppoisson_fit,-1,10,npar);
fitfcp->SetNpx(400);
fitfcp->SetLineWidth(1);
fitfcp->SetLineStyle(2);
fitfcp->SetLineColor(11);
fitfcp->SetParameter(0,norm);
fitfcp->SetParameter(1,offs);
fitfcp->SetParameter(2,sigm);
fitfcp->SetParameter(3,lamb);
fitfcp->SetParameter(4,pgeo);
fitfcp->SetParameter(5,gain);

cout <<"HERE 3!"<<endl;

TF1 *fitfcl = (TF1*) fitf->Clone();
fitfcl->SetRange(-1,10);
fitfcl->SetLineStyle(2);

cout <<"HERE 4!"<<endl;

//----------------------------------------------------
//-Drawing section
//----------------------------------------------------

double loframe = -1.0;
double hiframe =  10.0;
TH1F *frame = new TH1F("frame","",100, loframe,hiframe);
frame->GetXaxis()->SetTitle("charge (e^{-})");
frame->GetYaxis()->SetTitle("events");
frame->SetMinimum(1e-5);
frame->SetMaximum(histo->GetMaximum()*5);

TLatex *lat = new TLatex();
lat->SetNDC();
lat->SetTextSize(0.033);

TCanvas *canv = new TCanvas("canv","",600,400);
gPad->SetLogy();
frame->Draw();
histo->Draw("h same");
fitf->Draw("same");
fitfcl->Draw("same");
fitfcp->Draw("same");
lat->DrawLatex(0.7,0.84,Form("N = %4.1f #pm %3.1f",norm, norme));
lat->DrawLatex(0.7,0.81,Form("#delta = (%2.4f #pm %2.4f) e^{-}",offs, offse));
lat->DrawLatex(0.7,0.78,Form("#sigma = (%2.4f #pm %2.4f) e^{-}",sigm, sigme));
lat->DrawLatex(0.7,0.75,Form("#lambda = %2.3f #pm %2.3f",lamb, lambe));
lat->DrawLatex(0.7,0.72,Form("p = %2.3f #pm %2.3f",pgeo, pgeoe));
lat->DrawLatex(0.7,0.69,Form("g = %2.3f #pm %2.3f ",gain, gaine));
lat->DrawLatex(0.7,0.65,Form("#chi^{2}/NDF = %4.2f / %d ",chisq, ndegf));
lat->DrawLatex(0.7,0.62,Form("Prob = %2.4f ",proba));

lat->SetTextSize(0.05);
lat->DrawLatex(0.2, 0.91, "Convolution: Gauss #otimes (Poisson #otimes Geometric)");

canv->Print("fit_gausCONVcomppois.pdf");

}// end
