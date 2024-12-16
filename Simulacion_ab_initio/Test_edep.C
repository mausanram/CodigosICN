void Test_edep(){
//TFile *file = new TFile("Sim_ab_initio_NMUONS_300000.root");
// TFile *file = new TFile("Sim_ab_initio_NMUONS_400000.root");
// TFile *file = new TFile("Sim_ab_initio_NMUONS_50000_PLANES_1x1_RADIO_12_CCDSIZE_400x600_.root");
TFile *file = new TFile("Sim_ab_initio_NMUONS_50000_PLANES_1x1_RADIO_5_CCDSIZE_400x600_.root");
TTree *tree = (TTree*) file->Get("tree");


int NB = 90;
double tlow = 0;
// double thi = 620;
double thi = 1000;
TH1F *edep = new TH1F("edep", "", NB, tlow, thi);
TH1F *edep_cut = new TH1F("edep_cut", "", NB, tlow, thi);

// Fill histograms //
tree->Draw("edep>>edep", "l>0");
tree->Draw("edep>>edep_cut", "thet>22*TMath::Pi()/180 & edep>0");

// // Define fuctions //
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
edep->Draw();
// func1->Draw("same");

canv->cd(2);
edep_cut->Draw();
// func2->Draw("same");
canv->Print("Dis_edep.pdf");

}
