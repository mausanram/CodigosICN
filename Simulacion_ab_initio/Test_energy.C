void Test_energy(){
//TFile *file = new TFile("Sim_ab_initio_NMUONS_300000.root");
// TFile *file = new TFile("Sim_ab_initio_NMUONS_400000.root");
// TFile *file = new TFile("Sim_ab_initio_NMUONS_500000_PLANES_3.0x3.0_RADIO_12_0.root");

TFile *file = new TFile("Sim_ab_initio_Barra_NMUONS_300000_PLANES_150x150_RADIO_450_0.root");
TTree *tree = (TTree*) file->Get("tree");


int NB = 70;
double tlow = 0;
double thi = pow(10,5);
TH1F *energy_all = new TH1F("energy_all", "", NB, tlow, thi);
TH1F *energy_in = new TH1F("energy_in", "", NB, tlow, thi);

// Fill histograms //
tree->Draw("epri>>energy_all");
tree->Draw("epri>>energy_in", "l>0");

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
energy_all->Draw();
canv->SetLogy(1);
// func1->Draw("same");

canv->cd(2);
energy_in->Draw();
// func2->Draw("same");
canv->Print("Dis_enpri.pdf");

}
