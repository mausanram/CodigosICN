void Test_dedl(){

TFile *f_geant = new TFile("./root_files/muons_1M_vacuum_250x529_file_m_old_SDLog_nHG_1.root");
TTree *tree_geant = (TTree*) f_geant->Get("B02Evts");


TFile *f_pp = new TFile("../../../../Simulacion_ab_initio/Sim_ab_initio_NMUONS_1000000_PLANES_1.5_RADIO_8_CCDSIZE_250X529_SIGMA_1.0_C_0.root");
TTree *tree_pp = (TTree*) f_pp->Get("tree");

TLatex lat;

int NB = 150;
double tlow = 0;
double thi = 10000;

TH1F *dedl_g4 = new TH1F("dedl_g4", "Distribuci#acute{o}n de dE/dL", NB, tlow, thi);
dedl_g4->SetStats(0);
dedl_g4->SetLineStyle(1);
dedl_g4->SetLineColor(2);
dedl_g4->GetXaxis()->SetTitle("dE/dL (KeV/cm)");

TH1F *dedl_pp = new TH1F("dedl_pp", "Simulation PP", NB, tlow, thi);
dedl_pp->SetStats(0);
dedl_pp->SetLineStyle(1);
dedl_pp->SetLineColor(4);

TH1F *dedl_geant_scale = new TH1F("dedl_geant_scale", "Distribuci#acute{o}n de Energ#acute{i}as Depositadas (#theta > 20^{o})", NB, tlow, thi);
dedl_geant_scale->SetStats(0);
dedl_geant_scale->SetLineStyle(1);
dedl_geant_scale->SetLineColor(1);
dedl_geant_scale->GetXaxis()->SetTitle("Energ#acute{i}a (KeV)");

// ============= GEANT4 histograms =========== //
tree_geant->Draw("(EevtBar*1000)/LengthMuLAr>>dedl_g4", "LengthMuLAr>0"); 
tree_geant->Draw("(EevtBar*0.904*1000)/LengthMuLAr>>dedl_geant_scale", "LengthMuLAr>0"); 

// ============= SIM AB INITIO histograms =========== //
tree_pp->Draw("(edep*1000)/l>>dedl_pp", "l>0"); // SIM_AB_INITIO INFO


// double val_max_g4 = dedl_g4->FindFirstBinAbove(4000.);
double binmax_g4 = dedl_g4->GetMaximumBin();
double binmax_pp = dedl_pp->GetMaximumBin();
double binmax_g4scale = dedl_geant_scale->GetMaximumBin();

double max_value_g4 = dedl_g4->GetXaxis()->GetBinCenter(binmax_g4);
double max_value_pp = dedl_pp->GetXaxis()->GetBinCenter(binmax_pp);
double max_value_g4scale = dedl_geant_scale->GetXaxis()->GetBinCenter(binmax_g4scale);

cout<< "Max value G4: " << max_value_g4 << endl;
cout<< "Max value PP: " << max_value_pp << endl;
cout<< "Max value scaleG4: " << max_value_g4scale << endl;

// ============ Create Canvas ============== //
TCanvas *canv = new TCanvas("canv","Edep", 2*800, 600);
// canv->Divide(2,1);
canv->cd(1);

dedl_g4->Draw("hist same");
dedl_pp->Draw("hist same");
dedl_geant_scale->Draw("hist same");

TLegend *leg = new TLegend(0.5, 0.7, 0.9, 0.9);
leg->AddEntry(dedl_g4, "Simulaci#acute{o}n Geant4", "LP");
leg->AddEntry(dedl_pp, "Simulaci#acute{o}n ab initio", "LP");
leg->AddEntry(dedl_geant_scale, "Simulaci#acute{o}n de Geant4 (escalada: 0.904)", "LP");
leg->Draw();


//canv->cd(2);
//edep_cut->Draw();
// func2->Draw("same");
canv->Print("Dis_edep_GEANT-PP.pdf");

}
