void Test_edep(){
//TFile *file = new TFile("Sim_ab_initio_NMUONS_300000.root");
// TFile *file = new TFile("Sim_ab_initio_NMUONS_400000.root");
// TFile *file = new TFile("Sim_ab_initio_NMUONS_50000_PLANES_1x1_RADIO_12_CCDSIZE_400x600_.root");
// TFile *file = new TFile("Sim_ab_initio_NMUONS_50000_PLANES_1x1_RADIO_5_CCDSIZE_400x600_.root");

// TFile *file = new TFile("Sim_ab_initio_NMUONS_100000_PLANES_1x1_RADIO_5_CCDSIZE_400x600_SIGMA_LV_0.1_.root");


// =============== Archivos mas nuevos y que si funciona el ajuste ======================== //
// TFile *file = new TFile("Sim_ab_initio_NMUONS_20000_PLANES_2.4x2.4_RADIO_12_CCDSIZE_1058x1278_SIGMA_LV_1_.root");
// TFile *file = new TFile("Sim_ab_initio_NMUONS_30000_PLANES_1.5x1.5_RADIO_8_CCDSIZE_400x600_SIGMA_LV_1_.root");
// TFile *file = new TFile("Sim_ab_initio_NMUONS_100000_PLANES_1.5x1.5_RADIO_8_CCDSIZE_400x600_SIGMA_LV_1_.root");
// TFile *file = new TFile("Sim_ab_initio_NMUONS_100000_PLANES_1.5x1.5_RADIO_8_CCDSIZE_400x600_SIGMA_LV_0.3_.root");
TFile *file = new TFile("Sim_ab_initio_NMUONS_100000_PLANES_1.5x1.5_RADIO_8_CCDSIZE_400x600_SIGMA_LV_0.6_.root");


TTree *tree = (TTree*) file->Get("tree");

// TFile *file0 = new TFile("Edep_allclusters_NSAMP324_MeV.root");
// TTree *tree0 = (TTree*) file->Get("tree");

int NB = 90;
double tlow = 0;
// double thi = 620;
double thi = 1000;
TH1F *edep = new TH1F("edep", "", NB, tlow, thi);
edep->GetXaxis()->SetTitle("Energy (KeV)");


TH1F *edep_cut = new TH1F("edep_cut", "", NB, tlow, thi);

// TH1F *edep0 = new TH1F("edep0", "", NB, tlow, thi);

// Fill histograms //
tree->Draw("edep>>edep", "l>0");
tree->Draw("edep>>edep_cut", "thet>22*TMath::Pi()/180 & edep>0");

// tree0->Draw("edep>>edep0");

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
// edep_cut->Draw("same");
// func1->Draw("same");

canv->cd(2);
edep_cut->Draw();
// edep0->Draw("same");
// func2->Draw("same");
canv->Print("Dis_edep.pdf");

}
