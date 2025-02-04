
double Smith_Dull(double *lx , double *lpar){

    double E_mu = lx[0];
    double theta = TMath::Pi() * lpar[0]/180;

    //## ------------- Constantes físicas -------------- ##
    double k = 8.0 / 3.0;
    double b = 0.771;
    double lambda_pi = 120;    // g/cm^2
    double y_0 = 1000;      // g/cm^2
    double r = 0.76;
    double b_mu = 0.8;   

    // # if units == 0:
    double a = 2.5;      //## MeV cm^2/g
    double m_mu = 105.7;    //## MeV/c^2 
    double m_pi = 139.6;    //## MeV/c^2

    double tau_mu_0 = 2.2 * pow(10,-6);   //## s
    double tau_0 = 2.6 * pow(10, -8);     //## s
    double rho_0 = 0.00129; //## g/cm^3
    double c = TMath::C()* 100; //## cm/s

    //### ---------------------- Parámetros ---------------------- ###
    double E_pi = (1 / r) * (E_mu + a * y_0 * ((1/TMath::Cos(theta)) - 0.1));
    double B_mu = (b_mu * m_mu * y_0)/(tau_mu_0 * rho_0 * c);
    // double B_mu = (b_mu * m_mu * y_0 * c)/(tau_mu_0 * rho_0);
    double P_mu = pow((0.1 * TMath::Cos(theta)) * (1 - (a * (y_0 *(1/TMath::Cos(theta)) - 100))/( r * E_pi)), ((B_mu)/((r * E_pi + 100 * a) * TMath::Cos(theta))));
    double j_pi = (m_pi * y_0)/(c * tau_0 * rho_0);
    // double j_pi = (m_pi * y_0 * c)/(tau_0 * rho_0);

    //## Intensidad diferencial
    //# C_1 = E_pi ** (-k) * P_mu * lambda_pi * b * j_pi
    double C_1 = pow(E_pi, -k)* P_mu * lambda_pi * b * j_pi;
    double C_2 = E_pi * TMath::Cos(theta);
    double C_3 = b * j_pi;

    //# return (C_1 * np.sin(theta)) / (C_2 + C_3)
    double En_k = (C_1 ) / (C_2 + C_3);
    return En_k;
    }
/// ============================================================= ///



void Test_energy(){
//TFile *file = new TFile("Sim_ab_initio_NMUONS_300000.root");
// TFile *file = new TFile("Sim_ab_initio_NMUONS_400000.root");
// TFile *file = new TFile("Sim_ab_initio_NMUONS_500000_PLANES_3.0x3.0_RADIO_12_0.root");

// TFile *file = new TFile("Sim_ab_initio_NMUONS_300000_PLANES_1.5x1.5_RADIO_8_CCDSIZE_400x600_SIGMA_LV_1.0_.root");
TFile *file = new TFile("Sim_ab_initio_Barra_NMUONS_1000000_PLANES_150x150_RADIO_450.root");
TTree *tree = (TTree*) file->Get("tree");


int NB = 100;
double tlow = 1;
double thi = pow(10,5);
double xbins[NB+1];
double low_bin;

double Delta_log = (log10(thi) - log10(tlow))/(NB -1);
for (int i = 0; i<NB;i++){
    low_bin = pow(10, log10(tlow) + i * Delta_log);
    // cout<<low_bin <<endl;
    xbins[i] = low_bin;
};

xbins[NB] = thi;


TH1F *energy_all = new TH1F("energy_all", "", NB, xbins);
TH1F *energy_in = new TH1F("energy_in", "", NB, xbins);
// TH1F *energy_in = new TH1F("energy_in", "", NB, tlow, thi);

TH1F *thet_0_6 = new TH1F("thet_0_6", "", NB, xbins);
TH1F *thet_43_47= new TH1F("thet_43_47", "", NB, xbins);
TH1F *thet_73_77 = new TH1F("thet_73_77", "", NB, xbins);

// Fill histograms //
// tree->Draw("epri>>energy_all", "thet>3*TMath::Pi()/180 && thet<4*TMath::Pi()/180");
// tree->Draw("epri>>energy_in", "l>0");
tree->Draw("epri>>thet_0_6", "thet>0 && thet<3*TMath::Pi()/180");
tree->Draw("epri>>thet_43_47", "thet>44*TMath::Pi()/180 && thet<46*TMath::Pi()/180");
tree->Draw("epri>>thet_73_77", "thet>74*TMath::Pi()/180 && thet<76*TMath::Pi()/180");


// Define fuctions //
// double thet = 0;
// double thet = 45;
double thet = 75;

TF1 *Smith = new TF1("Smith", Smith_Dull, 1, pow(10,5), 1);
Smith->SetParameter(0, thet);
// TF1 *func1 = new TF1("func1", "[0]*sin(x)*(cos(x))^2", 0, 90);
// func1->SetParameter(0, 5000);

// TF1 *func2 = new TF1("func2", "[0]*sin(x)*(cos(x))^3+[1]*(sin(x))^2*(cos(x))^2", 0, 90);
// func2->SetParameter(0, 2000);
// func2->SetParameter(1, 1000);

// Fit functions //
// theta_all->Fit("func1", "R");
// theta_in->Fit("func2", "R");

double t0 = 0; // Degrees
double t45 = 45;
double t75 = 75;

TF1 *Smith0 = new TF1("Smith0", Smith_Dull, 1, pow(10,5), 1);
Smith0->SetParameter(0, t0);
double max0= Smith0->GetMaximum();
Smith0->SetLineColor(1);

TF1 *Smith45 = new TF1("Smith45", Smith_Dull, 1, pow(10,5), 1);
Smith45->SetParameter(0, t45);
double max45= Smith45->GetMaximum();
Smith45->SetLineColor(2);

TF1 *Smith75 = new TF1("Smith75", Smith_Dull, 1, pow(10,5), 1);
Smith75->SetParameter(0, t75);
double max75= Smith75->GetMaximum();
Smith75->SetLineColor(8);

double maxh0 = thet_0_6->GetMaximum();
thet_0_6->Scale((max0/maxh0));
double maxh45 = thet_43_47->GetMaximum();
thet_43_47->Scale(max45/maxh45);
double maxh75 = thet_73_77->GetMaximum();
thet_73_77->Scale(max75/maxh75);




// Create Canvas //
TCanvas *canv = new TCanvas("canv","", 2*700, 600);
// canv->Divide(2,1);

canv->cd(1);
thet_0_6->Draw("h0");
Smith0->Draw("L same");

// thet_43_47->Draw("h0");
// Smith45->Draw("L same");

// thet_73_77->Draw("h0");
// Smith75->Draw("L same");

// Smith0->Draw();
// Smith45->Draw("same");
// Smith75->Draw("same");

// thet_43_47->Draw("same");
// thet_43_47->Draw("same");

gPad->SetLogy(1);
gPad->SetLogx(1);
// func1->Draw("same");

// canv->cd(2);
// energy_in->Draw();
// gPad->SetLogy(1);
// gPad->SetLogx(1);
// func2->Draw("same");
canv->Print("Dis_enpri.pdf");

}
