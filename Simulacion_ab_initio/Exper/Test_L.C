void Test_L(){

TFile *f_connie = new TFile("..//Edep_muons_CONNIE_NSAMP400_700x420_SIGMAS_5_MeV.root"); // INFO ALL_CLUSTERS
TTree *tree_connie = (TTree*) f_connie->Get("tree");

TFile *f_ICN = new TFile("../Edep_NSAMP324_MeV.root"); // INFO ALL_CLUSTERS
TTree *tree_ICN = (TTree*) f_ICN->Get("tree");


int NB = 165;
// int NB = 120;
double tlow = 0;

double thi = 0.4;
// double thi = 30;

TH1F *L_ICN = new TH1F("L_ICN", "Distribucion de Longitudes", NB, tlow, thi);
L_ICN->GetXaxis()->SetTitle("Distancia (cm)");
L_ICN->SetStats(0);
L_ICN->SetLineColor(1);

TH1F *L_CON = new TH1F("L_CON", "Distribucion de Longitudes", NB, tlow, thi);
L_CON->GetXaxis()->SetTitle("Distancia (cm)");
L_CON->SetStats(0);
L_CON->SetLineColor(2);


TH1F *Lcut = new TH1F("Lcut", "", NB, tlow, thi);

// Fill histograms //
// tree_connie->Draw("l+0.0045>>L_CON", "l>0");
tree_connie->Draw("l>>L_CON", "l>0");
tree_ICN->Draw("l>>L_ICN", "l>0");
// tree->Draw("l>>Lcut", "l>0");

// Define fuctions //
// TF1 *func1 = new TF1("func1", "[0]*sin(x)*(cos(x))^2", 0, 90);
// func1->SetParameter(0, 5000);

// TF1 *func2 = new TF1("func2", "[0]*sin(x)*(cos(x))^3+[1]*(sin(x))^2*(cos(x))^2", 0, 90);
// func2->SetParameter(0, 2000);
// func2->SetParameter(1, 1000);

// Fit functions //
// theta_all->Fit("func1", "R");
// theta_in->Fit("func2", "R");

TLine *line = new TLine(0.0725,0,0.0725,1000);
// TLine *line = new TLine(0.068,0,0.068,30000);
line->SetLineStyle(2);
line->SetLineWidth(2);


// Create Canvas //
TCanvas *canv = new TCanvas("canv","", 2*700, 600);
// canv->Divide(2,1);

L_CON->Scale(1.6);
canv->cd(1);
L_ICN->Draw("hist");
L_CON->Draw("hist same");
line->Draw("same");
// Lcut->Draw("same");
// func1->Draw("same");

TLegend *leg = new TLegend(0.7, 0.8, 0.9, 0.9);
// leg->SetTextAlign(11);
leg->SetFillStyle(0);
leg->AddEntry(L_ICN, "Datos ICN", "l");
leg->AddEntry(L_CON, "Datos CONNIE", "l");
leg->AddEntry(line, "Grosor de la CCD: 0.0725 cm", "l");
// leg->AddEntry(line, "Grosor de la CCD: 0.068 cm", "l");
leg->Draw();

// canv->cd(2);
// Lcut->Draw();
// // func2->Draw("same");
canv->Print("Dis_L.pdf");

}
