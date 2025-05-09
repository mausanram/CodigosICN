void Test_thet(){
TFile *f_connie = new TFile("../Edep_muons_CONNIE_NSAMP400_700x420_SIGMAS_5_MeV.root"); // INFO ALL_CLUSTERS
TTree *tree_connie = (TTree*) f_connie->Get("tree");

TFile *f_ICN = new TFile("../Edep_NSAMP324_MeV.root"); // INFO ALL_CLUSTERS
TTree *tree_ICN = (TTree*) f_ICN->Get("tree");

TLatex lat;

// int NB = 90;
int NB = 60;
double tlow = 0;
double thi = TMath::Pi()/2.0;

TH1F *theta_icn = new TH1F("theta_icn", "Distribuci#acute{o}n angular #theta de los muones que impactaron la CCD", NB, tlow, thi);
theta_icn->GetXaxis()->SetTitle("#theta(rad)");
// theta_icn->SetStats(0);
theta_icn->SetLineColor(1);

TH1F *theta_con = new TH1F("theta_con", "Distribuci#acute{o}n angular #theta de los muones que impactaron la CCD", NB, tlow, thi);
theta_con->GetXaxis()->SetTitle("#theta(rad)");
// theta_con->SetStats(0);
theta_con->SetLineColor(2);


// TH1F *theta_incut = new TH1F("theta_incut", "", NB, tlow, thi);

// Fill histograms //
tree_ICN->Draw("thet>>theta_icn", "l>0");
tree_connie->Draw("thet*TMath::Pi()/180>>theta_con", "l>0");
// tree->Draw("thet>>theta_incut", "thet>22*TMath::Pi()/180");


// Create Canvas //
TCanvas *canv = new TCanvas("canv","", 2*900, 700);
// canv->Divide(2,1);
canv->cd(1);
canv->SetGrid(1);

theta_con->Scale(1.95);
theta_icn->Draw("hist same");
theta_con->Draw("hist same");


// ==== CONNIE ===   ///
TLegend *leg = new TLegend(0.7, 0.8, 0.9, 0.9);
leg->AddEntry(theta_icn, "Datos ICN", "f");
leg->AddEntry(theta_con, "Datos CONNIE", "f");
leg->Draw();

canv->Print("Dis_thet.pdf");

}
