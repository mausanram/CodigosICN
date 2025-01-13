void Test_edep(){
//TFile *file = new TFile("Sim_ab_initio_NMUONS_300000.root");
// TFile *file = new TFile("Sim_ab_initio_NMUONS_400000.root");
// TFile *file = new TFile("Sim_ab_initio_NMUONS_50000_PLANES_1x1_RADIO_12_CCDSIZE_400x600_.root");
// TFile *file = new TFile("Sim_ab_initio_NMUONS_50000_PLANES_1x1_RADIO_5_CCDSIZE_400x600_.root");

//TFile *file = new TFile("../build/muons_100K_vacuum_file.root");
TFile *file = new TFile("../build/B02ntuples.root");

TTree *tree = (TTree*) file->Get("B02Evts");
//TTree *tree = (TTree*) file->Get("B02Hits");

TFile *file0 = new TFile("Edep_allclusters_NSAMP324_MeV.root");
TTree *tree0 = (TTree*) file0->Get("tree");

TFile *file1 = new TFile("Sim_ab_initio_NMUONS_100000_PLANES_1x1_RADIO_5_CCDSIZE_400x600_SIGMA_LV_1_.root");
TTree *tree1 = (TTree*) file1->Get("tree");


int NB = 90;
double tlow = 0;
// double thi = 620;
double thi = 1;
TH1F *edep = new TH1F("edep", "", NB, tlow, thi);
edep->GetXaxis()->SetTitle("Energy (MeV)");


TH1F *edep_cut = new TH1F("edep_cut", "", NB, tlow, thi);

TH1F *edep0 = new TH1F("edep0", "", NB, tlow, thi);
TH1F *edep1 = new TH1F("edep1", "", NB, tlow, thi);

// Fill histograms //
//tree->Draw("WevtBar>>edep"); // GEANT4 INFO
tree->Draw(" EevtBar>>edep"); // GEANT4 INFO
//tree->Draw("Ehitbar>>edep_cut", "thet>22*TMath::Pi()/180 & edep>0");

tree0->Draw("edep>>edep0"); // EXPERIMENTAL INFO
tree1->Draw("edep/1000>>edep1"); // SIM_AB_INITIO INFO

// ========== Scale histograms =========== //
edep->Scale(0.75, "");
edep->SetLineColor(2);

edep0->Scale(1);
edep0->SetLineColor(4);

edep1->Scale(0.33);
edep1->SetLineColor(1);


// ======================================= //


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
edep->Draw("H");
edep0->Draw("H same");
edep1->Draw("H same");
// func1->Draw("same");

canv->cd(2);
//edep_cut->Draw();
// func2->Draw("same");
canv->Print("Dis_edep.pdf");

}
