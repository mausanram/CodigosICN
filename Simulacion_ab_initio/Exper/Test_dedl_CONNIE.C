void Test_dedl_CONNIE(){

TFile *f_connie = new TFile("../tree_muons_CONNIE_NSAMP400_1022x420_SIGMAS_5_MeV.root"); // INFO ALL_CLUSTERS
TTree *tree_CON = (TTree*) f_connie->Get("tree");


int NB = 200;
double tlow = 0;
double thi = 10000;

TH1F *L_CON = new TH1F("L_CON", "Distribucion de dE/dL", NB, tlow, thi);
L_CON->GetXaxis()->SetTitle("dE/dL (KeV/cm)");
L_CON->SetStats(0);
L_CON->SetLineColor(1);

// Fill histograms //
tree_CON->Draw("dedl>>L_CON");
// tree->Draw("l>>Lcut", "l>0");

// Create Canvas //
TCanvas *canv = new TCanvas("canv","", 2*700, 600);
// canv->Divide(2,1);
canv->cd(1);
L_CON->Draw("hist");

TLegend *leg = new TLegend(0.7, 0.8, 0.9, 0.9);
// leg->SetTextAlign(11);
leg->SetFillStyle(0);
leg->AddEntry(L_CON, "Datos CONNIE", "l");
leg->Draw();

// canv->cd(2);
// Lcut->Draw();
// // func2->Draw("same");
canv->Print("Dis_L.pdf");

}
