void Test_L(){
//TFile *file = new TFile("Sim_ab_initio_NMUONS_300000.root");
// TFile *file = new TFile("Sim_ab_initio_NMUONS_400000.root");
TFile *file = new TFile("Sim_ab_initio_NMUONS_100000_PLANES_3.0x3.0_RADIO_12_.root");

// TFile *file = new TFile("Sim_ab_initio_Barra_NMUONS_300000_PLANES_150x150_RADIO_450_1.root");
// TFile *file = new TFile("Sim_ab_initio_Barra_NMUONS_300000_PLANES_150x150_RADIO_450_0.root");
TTree *tree = (TTree*) file->Get("tree");


int NB = 160;
// int NB = 120;
double tlow = 0;

double thi = 0.4;
// double thi = 30;
TH1F *L = new TH1F("L", "", NB, tlow, thi);
TH1F *Lcut = new TH1F("Lcut", "", NB, tlow, thi);

// Fill histograms //
tree->Draw("l>>L", "l>0");
tree->Draw("l>>Lcut", "l>0 & thet>22*TMath::Pi()/180");
// tree->Draw("l>>Lcut", "l>0");

// Define fuctions //
// TF1 *func1 = new TF1("func1", "[0]*sin(x)*(cos(x))^2", 0, 90);
// func1->SetParameter(0, 5000);

// TF1 *func2 = new TF1("func2", "[0]*sin(x)*(cos(x))^3+[1]*(sin(x))^2*(cos(x))^2", 0, 90);
// func2->SetParameter(0, 2000);
// func2->SetParameter(1, 1000);

// Fit functions //
// theta_all->Fit("func1", "R");
// theta_in->Fit("func2", "R");

// Create Canvas //
TCanvas *canv = new TCanvas("canv","", 2*700, 600);
canv->Divide(2,1);
canv->cd(1);
L->Draw();
// Lcut->Draw("same");
// func1->Draw("same");

canv->cd(2);
Lcut->Draw();
// func2->Draw("same");
canv->Print("Dis_L.pdf");

}
