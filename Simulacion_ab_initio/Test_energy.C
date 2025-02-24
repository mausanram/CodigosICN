
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
    return En_k * 1000000000000;
    }
/// ============================================================= ///



void Test_energy(){

    TChain *tree = new TChain("tree");
tree->Add("Sim_ab_initio_NMUONS_1000000_PLANES_1.5_RADIO_8_CCDSIZE_400X600_C.root");
tree->Add("Sim_ab_initio_NMUONS_1000000_PLANES_1.5_RADIO_8_CCDSIZE_400X600_C_0.root");
tree->Add("Sim_ab_initio_NMUONS_1000000_PLANES_1.5_RADIO_8_CCDSIZE_400X600_C_1.root");
tree->Add("Sim_ab_initio_NMUONS_1000000_PLANES_1.5_RADIO_8_CCDSIZE_400X600_C_2.root");
tree->Add("Sim_ab_initio_NMUONS_1000000_PLANES_1.5_RADIO_8_CCDSIZE_400X600_C_3.root");
tree->Add("Sim_ab_initio_NMUONS_1000000_PLANES_1.5_RADIO_8_CCDSIZE_400X600_C_4.root");
tree->Add("Sim_ab_initio_NMUONS_1000000_PLANES_1.5_RADIO_8_CCDSIZE_400X600_C_5.root");
    // //TFile *file = new TFile("Sim_ab_initio_NMUONS_300000.root");
// // TFile *file = new TFile("Sim_ab_initio_NMUONS_400000.root");
// // TFile *file = new TFile("Sim_ab_initio_NMUONS_500000_PLANES_3.0x3.0_RADIO_12_0.root");

// // TFile *file = new TFile("Sim_ab_initio_NMUONS_300000_PLANES_1.5x1.5_RADIO_8_CCDSIZE_400x600_SIGMA_LV_1.0_.root");
// // TFile *file = new TFile("Sim_ab_initio_Barra_NMUONS_1000000_PLANES_150x150_RADIO_450.root");
// // TFile *file = new TFile("MuonGen_NMUONS_10000_pyth.root");
// TFile *file = new TFile("Sim_ab_initio_NMUONS_100000_PLANES_1.5_RADIO_8_CCDSIZE_400X600_C_0.root");
// TTree *tree = (TTree*) file->Get("tree");


int NB = 100;
double tlow = 0.1;
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


TH1F *thet_0_6 = new TH1F("thet_0_6", "", NB, xbins);
TH1F *thet_43_47= new TH1F("thet_43_47", "", NB, xbins);
TH1F *thet_73_77 = new TH1F("thet_73_77", "", NB, xbins);


// Fill histograms //
tree->Draw("epri>>thet_0_6", "thet>0 && thet<(6*TMath::Pi()/180)");
tree->Draw("epri>>thet_43_47", "thet>(43*TMath::Pi()/180) && thet<(47*TMath::Pi()/180)");
tree->Draw("epri>>thet_73_77", "thet>(73*TMath::Pi()/180) && thet<(77*TMath::Pi()/180)");

double cont = 0;
double w = 0;
double val = 0;

// for (int i = 1; i < NB; i++){
//     cont = thet_0_6->GetBinContent(i);
//     w = thet_0_6->GetBinWidth(i);
//     val = cont / w;
//     thet_0_6->SetBinContent(i, val);
// }

// for (int i = 1; i < NB; i++){
//     cont = thet_43_47->GetBinContent(i);
//     w = thet_43_47->GetBinWidth(i);
//     val  = cont / w;
//     thet_43_47->SetBinContent(i, val);
// }

// for (int i = 1; i < NB; i++){
//     cont = thet_73_77->GetBinContent(i);
//     w = thet_73_77->GetBinWidth(i);
//     val  = cont / w;
//     thet_73_77->SetBinContent(i, val);
// }

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
// thet_0_6->Scale((max0/maxh0));
thet_0_6->Scale(2.0);

double maxh45 = thet_43_47->GetMaximum();
// thet_43_47->Scale(max45/maxh45);
thet_43_47->Scale(10.0);


double maxh75 = thet_73_77->GetMaximum();
// thet_73_77->Scale(max75/maxh75);
thet_73_77->Scale(2.0);




// Create Canvas //
TCanvas *canv = new TCanvas("canv","", 2*700, 600);
// canv->Divide(2,1);

canv->cd(1);
Smith0->Draw("L");
thet_0_6->Draw("hist same");

// Smith45->Draw("L same");
// thet_43_47->Draw("hist same");

// Smith75->Draw("L same");
// thet_73_77->Draw("hist same");

// Smith0->Draw();
// Smith45->Draw("same");
// Smith75->Draw("same");

// gPad->SetLogy(1);
// gPad->SetLogx(1);
// func1->Draw("same");

// canv->cd(2);
// energy_in->Draw();
// gPad->SetLogy(1);
// gPad->SetLogx(1);
// func2->Draw("same");
canv->Print("Dis_enpri.pdf");

}
