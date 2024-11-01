void Test_L(){
//TFile *file = new TFile("Sim_ab_initio_NMUONS_300000.root");
TFile *file = new TFile("Sim_ab_initio_NMUONS_400000.root");
TTree *tree = (TTree*) file->Get("tree");


int NB = 100;
double tlow = 0;
double thi = 0.4;
TH1F *L = new TH1F("L", "", NB, tlow, thi);
TH1F *Lcut = new TH1F("Lcut", "", NB, tlow, thi);

// Fill histograms //
tree->Draw("l>>L", "edep>0");
tree->Draw("l>>Lcut", "edep>0 & thet>22*TMath::Pi()/180");

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
// func1->Draw("same");

canv->cd(2);
Lcut->Draw();
// func2->Draw("same");
canv->Print("Dis_L.pdf");

}
