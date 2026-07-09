void Test_edep(){

TFile *fileEdepCatalog = new TFile("../../../../../CONNIE/edepTreeCatalog.root"); // INFO ALL_CLUSTERS
TH1F *histEdepCat = (TH1F*)fileEdepCatalog->Get("edep");
histEdepCat->SetStats(0);

TFile *fileEdepE1Catalog = new TFile("../../../../../CONNIE/edepE1TreeCatalog.root"); // INFO ALL_CLUSTERS
TH1F *histEdepE1Cat = (TH1F*)fileEdepE1Catalog->Get("edepE1");
histEdepE1Cat->SetStats(0);
histEdepE1Cat->SetLineColor(4);

TFile *f_conn = new TFile("Edep_CONNIE_NSAMP400_MeV.root"); // INFO ALL_CLUSTERS
TTree *tree_conn = (TTree*) f_conn->Get("tree");

TLatex lat;

int NB = 150;
double tlow = 0;
double thi = 1000;


TH1F *edep_conn = new TH1F("edep_conn", "Energy Deposition Spectrum (CONNIE 2021-2022)", NB, tlow, thi);
edep_conn->GetXaxis()->SetTitle("Energ#acute{i}a (KeV)");
edep_conn->SetStats(0);
edep_conn->SetLineStyle(1);
edep_conn->SetLineColor(1);

tree_conn->Draw("edep*1000>>edep_conn", "edep>0"); 

edep_conn->Scale(1);
edep_conn->SetMaximum(2000);

double areaEdepConn = edep_conn->Integral(tlow, thi);
double areaEdepCat = histEdepCat->Integral(tlow, thi);
double areaEdepE1Cat = histEdepE1Cat->Integral(tlow, thi);
// edep_conn->Scale(10000);
// ============ Create Canvas ============== //
TCanvas *canv = new TCanvas("canv","Edep", 2*800, 600);
edep_conn->Draw("hist");
histEdepCat->Draw("hist same");
histEdepE1Cat->Draw("hist same");

TLegend *leg = new TLegend(0.5, 0.7, 0.9, 0.9);
leg->AddEntry(histEdepCat, "Catalog E0 Energy", "LEP");
leg->AddEntry((TObject*)0, Form("Events: %.0f \n", areaEdepCat), "");
leg->AddEntry(histEdepE1Cat, "Catalog E1 Energy", "LEP");
leg->AddEntry((TObject*)0, Form("Events: %.0f", areaEdepCat), "");
leg->AddEntry(edep_conn, "Clustering Python Method Energy", "LEP");
leg->AddEntry((TObject*)0, Form("Events: %.0f", areaEdepConn), "");
leg->Draw();


//canv->cd(2);
//edep_cut->Draw();
// func2->Draw("same");
canv->Print("Dis_edep_GEANT-PP.pdf");

}
