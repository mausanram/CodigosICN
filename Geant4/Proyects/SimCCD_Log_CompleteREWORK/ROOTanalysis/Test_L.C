void Test_L(){

TChain *tree = new TChain("B02Evts");
//tree->Add("./root_files/muons_100K_SimCOMPLETE_100R_5P_0.root");
tree->Add("../../../../../Data_Geant4/Version_R100_P5/muons_500K_SimCOMPLETE_*");
//tree->Add("../../../../../Data_Geant4/Version_R20_P2/muons_500K_SimCOMPLETE_2.root");
//tree->Add("../../../../../Data_Geant4/Version_R20_P2/muons_500K_SimCOMPLETE_1.root");
//tree->Add("../../../../../Data_Geant4/Version_R20_P2/muons_500K_SimCOMPLETE_3.root");
//tree->Add("../../../../../Data_Geant4/Version_R20_P2/muons_500K_SimCOMPLETE_4.root");
//tree->Add("../../../../../Data_Geant4/Version_R20_P2/muons_500K_SimCOMPLETE_5.root");
//tree->Add("../../../../../Data_Geant4/Version_R20_P2/muons_500K_SimCOMPLETE_6.root");

// TTree *tree = (TTree*) file->Get("B02Evts");

TChain *tree_g4CCD = new TChain("B02Evts");
tree_g4CCD->Add("../../SimCCD_SDLog/ROOTanalysis/root_files/muons_1M_vacuum_250x529_file_m_old_SDLog_nHG_10.root");


TFile *file0 = new TFile("../../../../Simulacion_ab_initio/Sim_ab_initio_NMUONS_1000000_PLANES_1.5_RADIO_8_CCDSIZE_250X529_SIGMA_1.0_C_0.root");
TTree *tree0 = (TTree*) file0->Get("tree");


int NB = 150;
double tlow = 0.;
double thi = 0.4;

// TH1F *L = new TH1F("L", "Distance Distribution (CCD size: 0.9x0.6x0.0725 cm)", NB, tlow, thi);
// L->GetXaxis()->SetTitle("Distance (cm)");
// TH1F *L = new TH1F("L", "Distribuci#acute{o}n de Longitudes con corte angular de 25^{o}", NB, tlow, thi);
TH1F *L = new TH1F("L", "Distribuci#acute{o}n de Longitudes", NB, tlow, thi);
L->GetXaxis()->SetTitle("Distancia (cm)");
// L->SetStats(0);
L->SetLineColor(2);

TH1F *L_ccd = new TH1F("L_ccd", "Distribuci#acute{o}n de Longitudes", NB, tlow, thi);
L_ccd->GetXaxis()->SetTitle("Distancia (cm)");
// L->SetStats(0);
L_ccd->SetLineColor(1);

// TH1F *LPP = new TH1F("LPP", "Distribuci#acute{o}n de Longitudes con corte angular de 25^{o}", NB, tlow, thi);
TH1F *LPP = new TH1F("LPP", "Distribuci#acute{o}n de Longitudes", NB, tlow, thi);
LPP->GetXaxis()->SetTitle("Distancia (cm)");
LPP->SetStats(0);
LPP->SetLineStyle(1);
LPP->SetLineColor(4);

TH1F *Lcut = new TH1F("Lcut", "", NB, tlow, thi);

// Fill histograms //
// tree->Draw("LengthMuLAr>>L", "LengthMuLAr>0");  
// tree->Draw("LengthMuLAr>>L", "nHitBar>0");
// tree->Draw("LengthMuLAr>>L", "nHitBar>0 && LengthMuLAr>0");

tree_g4CCD->Draw("LengthMuLAr>>L_ccd" , "LengthMuLAr>0");
tree0->Draw("l>>LPP", "l > 0");
tree->Draw("LengthMuLAr>>L" , "LengthMuLAr>0");
// tree->Draw("LengthMuLAr>>L", "LengthMuLAr>0 && thetaPri>20*TMath::Pi()/180");
// tree0->Draw("l>>LPP", "l > 0 && thet>25*TMath::Pi()/180");



//TLine *line = new TLine(0.0725,0,0.0725,30000);
TLine *line = new TLine(0.0725,0,0.0725,8000);
line->SetLineStyle(2);
line->SetLineWidth(2);

// LPP->Scale(0.94);
// LPP->Scale(0.1);
// L->Scale(10.);
// L_ICN->Scale(6);

// Create Canvas //
TCanvas *canv = new TCanvas("canv","", 2*700, 600);
canv->Divide(2,1);
canv->cd(1);
L->Draw("hist");
L_ccd->Draw("hist same");
// LPP->Draw("hist same");
// line->Draw("same");
// Lcut->Draw("same");

TLegend *leg = new TLegend(0.5, 0.7, 0.8, 0.8);
leg->AddEntry(L, "Simulaci#acute{o}n Geant4", "lp");
leg->AddEntry(LPP, "Simulaci#acute{o}n ab initio", "lp");
leg->AddEntry(L_ccd, "Simulaci#acute{o}n G4 (solo CCD)", "lp");
// leg->AddEntry(L_ICN, "Datos ICN", "l");
leg->AddEntry(line, "Grosor de CCD: 0.0725 cm", "l");
leg->Draw();

//canv->cd(2);
//Lcut->Draw();
// func2->Draw("same");
canv->Print("Dis_L.pdf");

}
