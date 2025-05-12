void Test_SolvsElip(){

TFile *f_connie = new TFile("../tree_muons_ICN_NSAMP324_250x539_SIGMAS_20_MeV.root"); // INFO ALL_CLUSTERS
TTree *tree_connie = (TTree*) f_connie->Get("tree");

TFile *f_connie_nonmuons = new TFile("../tree_nonmuons_ICN_NSAMP324_250x539_SIGMAS_20_MeV.root"); // INFO ALL_CLUSTERS
TTree *tree_connie_nonmuons = (TTree*) f_connie_nonmuons->Get("tree");

int NBE = 150;
double elo = 0.;
double ehi = 1.;
int NBS = 150;
double slo = 0.;
double shi = 1.;

TH2F *histLT = new TH2F("histLT", "Elip vs Sol", NBE, elo, ehi, NBS, slo, shi);
histLT->GetXaxis()->SetTitle("Elipticity");
histLT->GetYaxis()->SetTitle("Solidity");
histLT->SetStats(0);

TH2F *histLTnonmuons = new TH2F("histLTnonmuons", "Elip vs Sol", NBE, elo, ehi, NBS, slo, shi);
histLTnonmuons->GetXaxis()->SetTitle("Elipticity");
histLTnonmuons->GetYaxis()->SetTitle("Solidity");
histLTnonmuons->SetStats(0);

tree_connie->Draw("sol:elip>>histLT", "l > 0");
tree_connie_nonmuons->Draw("sol:elip>>histLTnonmuons");

TLine *line = new TLine(0.65,0.7,0.65,1);
// TLine *line = new TLine(0.068,0,0.068,30000);
line->SetLineStyle(2);
line->SetLineWidth(2);

TLine *line_0 = new TLine(0.65,0.7, 1,0.7);
// TLine *line = new TLine(0.068,0,0.068,30000);
line_0->SetLineStyle(2);
line_0->SetLineWidth(2);


// Create Canvas //
TCanvas *canv = new TCanvas("canv","", 2*700, 600);
canv->Divide(2,1);

canv->cd(1);
histLT->Draw("colZ");
// histLTnonmuons->Draw("same colZ");
line->Draw("same");
line_0->Draw("same");

canv->cd(2);
histLTnonmuons->Draw("same colZ");
line->Draw("same");
line_0->Draw("same");

TLegend *leg = new TLegend(0.7, 0.8, 0.9, 0.9);
// leg->SetTextAlign(11);
leg->SetFillStyle(0);
leg->AddEntry(histLT, "Datos CONNIE", "l");
leg->AddEntry(line, "Grosor de la CCD: 0.0725 cm", "l");
// leg->AddEntry(line, "Grosor de la CCD: 0.068 cm", "l");
// leg->Draw();

// canv->cd(2);
// Lcut->Draw();
// // func2->Draw("same");
canv->Print("Dis_L.pdf");

}
