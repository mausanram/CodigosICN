void Test_dedl(){

TFile *f_icn = new TFile("../tree_muons_ICN_NSAMP324_250x539_SIGMAS_13_MeV.root"); // INFO ALL_CLUSTERS
TTree *tree_ICN = (TTree*) f_icn->Get("tree");


int NB = 200;
double tlow = 0;
double thi = 10000;

TH1F *L_ICN = new TH1F("L_ICN", "Distribucion de dE/dL", NB, tlow, thi);
L_ICN->GetXaxis()->SetTitle("dE/dL (KeV/cm)");
L_ICN->SetStats(0);
L_ICN->SetLineColor(1);

// Fill histograms //
tree_ICN->Draw("dedl>>L_ICN");
// tree->Draw("l>>Lcut", "l>0");

// Create Canvas //
TCanvas *canv = new TCanvas("canv","", 2*700, 600);
//canv->Divide(2,1);
canv->cd(1);
L_ICN->Draw("hist");

TLegend *leg = new TLegend(0.7, 0.8, 0.9, 0.9);
// leg->SetTextAlign(11);
leg->SetFillStyle(0);
leg->AddEntry(L_ICN, "Datos ICN", "l");
leg->Draw();

// canv->cd(2);
// Lcut->Draw();
// // func2->Draw("same");
canv->Print("Dis_L.pdf");

}
