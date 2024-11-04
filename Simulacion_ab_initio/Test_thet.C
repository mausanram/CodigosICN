void Test_thet(){
//TFile *file = new TFile("Sim_ab_initio_NMUONS_300000.root");
TFile *file = new TFile("Sim_ab_initio_NMUONS_400000.root");
TTree *tree = (TTree*) file->Get("tree");


int NB = 70;
double tlow = 0;
double thi = TMath::Pi()/2.0;
TH1F *theta_all = new TH1F("theta_all", "", NB, tlow, thi);
TH1F *theta_in = new TH1F("theta_in", "", NB, tlow, thi);
// TH1F *theta_incut = new TH1F("theta_incut", "", NB, tlow, thi);

// Fill histograms //
tree->Draw("thet>>theta_all");
tree->Draw("thet>>theta_in", "edep>0");
// tree->Draw("thet>>theta_incut", "thet>22*TMath::Pi()/180");

// Define fuctions //
TF1 *func1 = new TF1("func1", "[0]*sin(x)*(cos(x))^2", 0, 90*TMath::Pi()/180);
func1->SetParameter(0, 5000);

// TF1 *func2 = new TF1("func2", "[0]*(1899639*sin(x)*(cos(x))^3+(161472/TMath::Pi())*(sin(x))^2*(cos(x))^2)", 0, 90*TMath::Pi()/180);

// TF1 *func2 = new TF1("func2", "([0]*sin(x)*(cos(x))^3+([1])*(sin(x))^2*(cos(x))^2)", 0, 90*TMath::Pi()/180);
// TF1 *func2 = new TF1("func2", "[0]*((12725.6)*sin(x)*(cos(x))^3+(373.769)*(sin(x))^2*(cos(x))^2)", 0, 90*TMath::Pi()/180);
// TF1 *func2 = new TF1("func2", "[0]*(7598556*sin(x)*(cos(x))^3+(645888/TMath::Pi())*(sin(x))^2*(cos(x))^2)", 0, 90*TMath::Pi()/180);

TF1 *func2 = new TF1("func2", "[0]*((3.29909/1)*sin(x)*(cos(x))^3+(2.8041616/(8*TMath::Pi()))*(sin(x))^2*(cos(x))^2)", 0, 90*TMath::Pi()/180);
// TF1 *func2 = new TF1("func2", "[0]*((7.598556/4)*sin(x)*(cos(x))^3+(0.80896/(4*TMath::Pi()))*(sin(x))^2*(cos(x))^2)", 0, 90*TMath::Pi()/180); // El mejor 
func2->SetParameter(0, 2000);

// Fit functions //
theta_all->Fit("func1", "R");
theta_in->Fit("func2", "R");

double Prob1 = func1->GetProb();
double chi1 = func1->GetChisquare();

double Prob2 = func2->GetProb();
double chi2 = func2->GetChisquare();

std::cout << "ChiSquare1 = "<< chi1 << std::endl;
std::cout << "Prob1 = "<< Prob1 << std::endl;

std::cout << "ChiSquare2 = "<< chi2 << std::endl;
std::cout << "Prob2 = "<< Prob2 << std::endl;


// Create Canvas //
TCanvas *canv = new TCanvas("canv","", 2*700, 600);
canv->Divide(2,1);
canv->cd(1);
theta_all->Draw();
func1->Draw("same");

canv->cd(2);
theta_in->Draw();
func2->Draw("same");
canv->Print("Dis_thet.pdf");

}
