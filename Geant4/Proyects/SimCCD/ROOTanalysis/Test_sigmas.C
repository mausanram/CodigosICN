void Test_sigmas(){

TFile *file = new TFile("../../../../Simulacion_ab_initio/Sim_ab_initio_NMUONS_100000_PLANES_1.5x1.5_RADIO_8_CCDSIZE_400x600_SIGMA_LV_1.0_.root");
TTree *tree = (TTree*) file->Get("tree");

TFile *file0 = new TFile("../../../../Simulacion_ab_initio/Sim_ab_initio_NMUONS_100000_PLANES_1.5x1.5_RADIO_8_CCDSIZE_400x600_SIGMA_LV_0.9_.root"); 
TTree *tree0 = (TTree*) file0->Get("tree");

TFile *file1 = new TFile("../../../../Simulacion_ab_initio/Sim_ab_initio_NMUONS_100000_PLANES_1.5x1.5_RADIO_8_CCDSIZE_400x600_SIGMA_LV_0.8_.root");
TTree *tree1 = (TTree*) file1->Get("tree");

TFile *file2 = new TFile("../../../../Simulacion_ab_initio/Sim_ab_initio_NMUONS_100000_PLANES_1.5x1.5_RADIO_8_CCDSIZE_400x600_SIGMA_LV_0.7_.root");
TTree *tree2 = (TTree*) file2->Get("tree");

TFile *file3 = new TFile("../../../../Simulacion_ab_initio/Sim_ab_initio_NMUONS_100000_PLANES_1.5x1.5_RADIO_8_CCDSIZE_400x600_SIGMA_LV_0.6_.root");
TTree *tree3 = (TTree*) file3->Get("tree");


int NB = 90;
double tlow = 0.001;
double thi = 1000;

TH1F *edep = new TH1F("edep", "Energy Distribution (Sim. PP)", NB, tlow, thi);
edep->GetXaxis()->SetTitle("Energy (MeV)");
edep->SetLineStyle(1);
edep->SetLineColor(1);
edep->SetStats(0);

TH1F *edep0 = new TH1F("edep0", "", NB, tlow, thi);
edep0->SetLineStyle(1);
edep0->SetLineColor(2);
edep0->SetStats(0);

TH1F *edep1 = new TH1F("edep1", "", NB, tlow, thi);
edep1->SetStats(0);
edep1->SetLineStyle(1);
edep1->SetLineColor(3);

TH1F *edep2 = new TH1F("edep2", "", NB, tlow, thi);
edep2->SetStats(0);
edep2->SetLineStyle(1);
edep2->SetLineColor(4);

TH1F *edep3 = new TH1F("edep3", "Energy Distribution (Sim. PP)", NB, tlow, thi);
edep3->GetXaxis()->SetTitle("Energy (KeV)");
edep3->SetStats(0);
edep3->SetLineStyle(1);
edep3->SetLineColor(5);

// ============= Fill histograms =========== //
tree->Draw("edep>>edep"); 
tree0->Draw("edep>>edep0"); // 
tree1->Draw("edep>>edep1"); 
tree2->Draw("edep>>edep2"); // 
tree3->Draw("edep>>edep3"); 

// ========== Scale histograms =========== //
/*edep->Scale(4.5);
edep0->Scale(9.1);
edep1->Scale(6.95);
edep2->Scale(0.86);
edep3->Scale(6.);*/
// ======================================= //


// ============ Create Canvas ============== //
TCanvas *canv = new TCanvas("canv","Edep", 2*800, 600);
canv->Divide(1,1);
canv->cd(1);
edep3->Draw("h"); 		// edepG4 with birks
edep0->Draw("he0 same"); 	// ICN data 
edep1->Draw("he0 same"); 	// edepG4 no Birks
edep2->Draw("he0 same"); 	// CONNIE data
edep->Draw("he0 same"); 	// edepPP

TLegend *leg = new TLegend(0.5, 0.7, 0.9, 0.9);
leg->AddEntry(edep, "Sigma: 1.0", "lep");
leg->AddEntry(edep1, "Sigma: 0.9", "LEP");
leg->AddEntry(edep0, "Sigma: 0.8", "LEP");
leg->AddEntry(edep2, "Sigma: 0.7", "LEP");
leg->AddEntry(edep3, "Sigma: 0.6", "LEP");
leg->Draw();


//canv->cd(2);
//edep_cut->Draw();
// func2->Draw("same");
canv->Print("Dis_edep_PP_SIGMAS.pdf");

}
