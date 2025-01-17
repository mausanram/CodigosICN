void Test_edep(){
//TFile *file = new TFile("./root_files/muons_100K_vacuum_file.root");
TFile *file = new TFile("./root_files/muons_100K_vacuum_file.root");
TTree *tree = (TTree*) file->Get("B02Evts");
//TTree *tree = (TTree*) file->Get("B02Hits");

TFile *file0 = new TFile("Edep_allclusters_NSAMP324_MeV.root");
TTree *tree0 = (TTree*) file0->Get("tree");

TFile *file1 = new TFile("Sim_ab_initio_NMUONS_100000_PLANES_1x1_RADIO_5_CCDSIZE_400x600_SIGMA_LV_1_.root");
TTree *tree1 = (TTree*) file1->Get("tree");

TFile *file2 = new TFile("Edep_CONNIE_NSAMP400_MeV.root");
TTree *tree2 = (TTree*) file2->Get("tree");


int NB = 90;
double tlow = 0;
// double thi = 620;
double thi = 1;
TH1F *edep = new TH1F("edep", "", NB, tlow, thi);
edep->GetXaxis()->SetTitle("Energy (MeV)");
edep->SetStats(0);


TH1F *edep_cut = new TH1F("edep_cut", "", NB, tlow, thi);
edep_cut->SetStats(0);

TH1F *edep0 = new TH1F("edep0", "", NB, tlow, thi);
edep0->SetStats(0);

TH1F *edep1 = new TH1F("edep1", "", NB, tlow, thi);
edep1->SetStats(0);

TH1F *edep2 = new TH1F("edep2", "", NB, tlow, thi);
edep2->SetStats(0);

// ============= Fill histograms =========== //
tree->Draw("WevtBar>>edep"); // GEANT4 INFO
tree->Draw(" EevtBar>>edep1"); // GEANT4 INFO
//tree->Draw("Ehitbar>>edep_cut", "thet>22*TMath::Pi()/180 & edep>0");
// ========================================= //

tree0->Draw("edep>>edep0"); // EXPERIMENTAL INFO
//tree1->Draw("edep/1000>>edep1"); // SIM_AB_INITIO INFO
tree2->Draw("edep>>edep2"); // EXPERIMENTAL INFO CONNIE

// ========== Scale histograms =========== //
//edep->Scale(0.65, "");
edep->Scale(4.5, "");
edep->SetLineColor(2);

edep0->Scale(1);
edep0->SetLineColor(4);

edep1->Scale(8.);
edep1->SetLineColor(1);

edep2->Scale(0.86);
edep2->SetLineColor(7);
// ======================================= //


// ======== Fits ========== //

// ======================== //

// Create Canvas //
TCanvas *canv = new TCanvas("canv","", 2*700, 600);
canv->Divide(2,1);
canv->cd(1);
edep->Draw("H");
edep0->Draw("H same");
edep1->Draw("H same");
edep2->Draw("H same");

TLegend *leg = new TLegend(0.5, 0.7, 0.9, 0.9);
leg->AddEntry(edep, "Sim-Wedep", "lep");
leg->AddEntry(edep1, "Sim-Edep", "LEP");
leg->AddEntry(edep0, "Datos ICN-NSAMP324", "LEP");
leg->AddEntry(edep2, "Datos CONNIE-NSAMP400", "LEP");
leg->Draw();


//canv->cd(2);
//edep_cut->Draw();
// func2->Draw("same");
canv->Print("Dis_edep_GEANT-PP-EXP.pdf");

}
