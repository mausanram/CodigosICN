void Test_edep_new(){

// TFile *f_geant = new TFile("./root_files/muons_2M_vacuum_file.root");
// TFile *f_geant = new TFile("./root_files/muons_1M_vacuum_250x529_file_m.root");
TFile *f_geant = new TFile("./root_files/muons_1M_vacuum_250x529_file_m_old_SDLog_nHG_10.root");
// TFile *file = new TFile("./root_files/muons_1M_vacuum_file.root");
TTree *tree_geant = (TTree*) f_geant->Get("B02Evts");
//TTree *tree = (TTree*) file->Get("B02Hits");

TLatex lat;

int NB = 150;
double tlow = -100;
double thi = 1000;

//TH1F *edep0 = new TH1F("edep0", "Lab. Det. ", NB, tlow, thi);
TH1F *edep_icn = new TH1F("edep_icn", "Espectro de Energ#acute{i}as Depositada (ICN)", NB, tlow, thi);
// TH1F *edep_icn = new TH1F("edep_icn", "Espectro de Energ#acute{i}as Depositada (#theta > 20)", NB, tlow, thi);
edep_icn->GetXaxis()->SetTitle("Energ#acute{i}a (KeV)");
edep_icn->SetLineStyle(1);
edep_icn->SetLineColor(1);
edep_icn->SetStats(0);

TH1F *edep_g4 = new TH1F("edep_g4", "Distribuci#acute{o}n de Energ#acute{i}a Depositada", NB, tlow, thi);
edep_g4->SetStats(0);
edep_g4->SetLineStyle(1);
edep_g4->SetLineColor(2);
edep_g4->GetXaxis()->SetTitle("Energ#acute{i}a (KeV)");

TH1F *edep_geant_scale = new TH1F("edep_geant_scale", "Distribuci#acute{o}n de Energ#acute{i}as Depositadas (#theta > 20^{o})", NB, tlow, thi);
edep_geant_scale->SetStats(0);
edep_geant_scale->SetLineStyle(1);
edep_geant_scale->SetLineColor(2);
edep_geant_scale->GetXaxis()->SetTitle("Energ#acute{i}a (KeV)");

// ============= GEANT4 histograms =========== //
tree_geant->Draw("EevtBar>>edep_g4"); // GEANT4 INFO (NO BIRKS)

// ============ Create Canvas ============== //
TCanvas *canv = new TCanvas("canv","Edep", 2*800, 600);
canv->Divide(2,1);
canv->cd(1);
edep_g4->Draw("hist");

TLegend *leg = new TLegend(0.5, 0.7, 0.9, 0.9);
leg->AddEntry(edep_icn, "Datos de muones ICN (NSAMP324)", "LP");
leg->AddEntry(edep_geant_scale, "Simulaci#acute{o}n de Geant4 (escalada: 0.904)", "LP");
leg->Draw();


//canv->cd(2);
//edep_cut->Draw();
// func2->Draw("same");
canv->Print("Dis_edep_GEANT-PP.pdf");

}
