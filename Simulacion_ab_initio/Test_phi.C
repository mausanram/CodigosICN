void Test_phi(){
TFile *file = new TFile("Sim_ab_initio_NMUONS_100000.root");
TTree *tin = (TTree*) file->Get("tree");
TTree *tall = (TTree*) file->Get("tree_1");


int NB = 70;
double tlow = 0;
double thi = 3.1415/2.0;
TH1F *theta_all = new TH1F("theta_all", "", NB, tlow, thi);
TH1F *theta_in = new TH1F("theta_in", "", NB, tlow, thi);
TH1F *theta_incut = new TH1F("theta_incut", "", NB, tlow, thi);

// Fill histograms //
tall->Draw("thetall*3.1415/180>>theta_all");
tin->Draw("thet*3.1415/180>>theta_in");
tin->Draw("thet*3.1415/180>>theta_incut", "thet>22*3.1415/180");

// Define fuctions //
TF1 *func1 = new TF1("func1", "[0]*sin(x)*(cos(x))^2", 0, 90);
func1->SetParameter(0, 5000);

TF1 *func2 = new TF1("func2", "[0]*sin(x)*(cos(x))^3+[1]*(sin(x))^2*(cos(x))^2", 0, 90);
func2->SetParameter(0, 2000);
func2->SetParameter(1, 1000);

// Fit functions //
theta_all->Fit("func1", "R");
theta_in->Fit("func2", "R");

// Create Canvas //
TCanvas *canv = new TCanvas("canv","", 2*700, 600);
canv->Divide(2,1);
canv->cd(1);
theta_all->Draw();
func1->Draw("same");

canv->cd(2);
theta_in->Draw();
func2->Draw("same");
canv->Print("plot.pdf");

}
