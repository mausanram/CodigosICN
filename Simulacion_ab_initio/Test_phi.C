void Test_phi(){

TChain *tree = new TChain("tree");
// tree->Add("Sim_ab_initio_NMUONS_1000000_PLANES_1.5_RADIO_8_CCDSIZE_400X600_C.root");
tree->Add("Sim_ab_initio_NMUONS_2000000_PLANES_1.5_RADIO_8_CCDSIZE_400X600_C_0.root");
// tree->Add("Sim_ab_initio_NMUONS_1000000_PLANES_1.5_RADIO_8_CCDSIZE_400X600_C_0.root");
// tree->Add("Sim_ab_initio_NMUONS_1000000_PLANES_1.5_RADIO_8_CCDSIZE_400X600_C_1.root");
// tree->Add("Sim_ab_initio_NMUONS_1000000_PLANES_1.5_RADIO_8_CCDSIZE_400X600_C_2.root");
// tree->Add("Sim_ab_initio_NMUONS_1000000_PLANES_1.5_RADIO_8_CCDSIZE_400X600_C_3.root");
// tree->Add("Sim_ab_initio_NMUONS_1000000_PLANES_1.5_RADIO_8_CCDSIZE_400X600_C_4.root");
// tree->Add("Sim_ab_initio_NMUONS_1000000_PLANES_1.5_RADIO_8_CCDSIZE_400X600_C_5.root");


// TFile *file = new TFile("Sim_ab_initio_NMUONS_300000_PLANES_1.5x1.5_RADIO_8_CCDSIZE_400x600_SIGMA_LV_1.0_.root");
// TFile *file = new TFile("Sim_ab_initio_NMUONS_200000_PLANES_1.5_RADIO_8_CCDSIZE_400X600_C.root");
// TFile *file = new TFile("Sim_ab_initio_NMUONS_1000000_PLANES_1.5_RADIO_8_CCDSIZE_400X600_C_0.root");
// TTree *tree = (TTree*) file->Get("tree");


int NB = 70;
double tlow = 0;
double thi = 2*TMath::Pi() + 0.1;
TH1F *theta_all = new TH1F("phi_all", "", NB, tlow, thi);
TH1F *theta_in = new TH1F("phi_in", "", NB, tlow, thi);
// TH1F *theta_incut = new TH1F("theta_incut", "", NB, tlow, thi);

// Fill histograms //
tree->Draw("phi>>phi_all");
tree->Draw("phi>>phi_in", "edep > 0");
// tree->Draw("phi>>theta_incut", "thet>22");

// Define fuctions //
TF1 *func1 = new TF1("func1", "[0]", 0.1, 2*TMath::Pi() - 0.1);
func1->SetParameter(0, 100000);

double Ah = 0.54;
double Al = 0.06525;
double Ac = 0.0435;

// ======================================= Funcionde ajuste para la CCD  =============================================== ###
// TF1 *func2 = new TF1("func2", Form("[0]*(((%f)/(1 * TMath::Pi())) + ((%f)/2)*abs(cos(x)) + ((%f)/2)*abs(sin(x)))", Ah, Al, Ac), 0.05,2*TMath::Pi() - 0.1); 
TF1 *func2 = new TF1("func2", Form("[0]*(((%f)/(1 * TMath::Pi())) + ((%f)/4)*abs(cos(x)) + ((%f)/4)*abs(sin(x)))", Ah, Al, Ac), 0.05,2*TMath::Pi() - 0.1); 

// ======================================= Funcionde ajuste para la barra =============================================== ###
// TF1 *func2 = new TF1("func2", "[0]*((20/TMath::Pi()) + (5)*abs(cos(x)) + (1/2)*abs(sin(x)))",0.01,2*TMath::Pi() - 0.1);


// TF1 *func2 = new TF1("func2", "[0]*((8 * 0.8247739/(2 * TMath::Pi())) + (4 * 0.0499529/4)*(cos(x)) + (4 * 0.0376767715/4)*(sin(x)))", 0.01,2*TMath::Pi() - 0.1); 
// TF1 *func2 = new TF1("func2", "[0]*((8 * 0.8247739/(2 * TMath::Pi())) + (4 * 0.0499529/4)*abs(cos(x)))", 0.01,2*TMath::Pi() - 0.1); // EL mejor
// func2->SetParameter(0, 1600);
//func2->SetParameter(1, 1000);

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
theta_in->SetMinimum(0);
theta_in->Draw();
func2->Draw("same");
canv->Print("Dis_phi.pdf");

}
