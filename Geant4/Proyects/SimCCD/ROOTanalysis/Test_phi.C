void Test_phi(){
//TFile *file = new TFile("Sim_ab_initio_NMUONS_300000.root");
// TFile *file = new TFile("Sim_ab_initio_NMUONS_400000.root");
// TFile *file = new TFile("Sim_ab_initio_NMUONS_10000_PLANES_3.0x3.0_RADIO_12_.root");
TFile *file = new TFile("./root_files/muons_2M_vacuum_file.root");

// TFile *file = new TFile("Sim_ab_initio_Barra_NMUONS_200000_PLANES_150x150_RADIO_100.root");
// TFile *file = new TFile("Sim_ab_initio_Barra_NMUONS_500000_PLANES_150x150_RADIO_100(1).root");
// TFile *file = new TFile("Sim_ab_initio_Barra_NMUONS_100000_PLANES_150x150_RADIO_450_0.root");
// TFile *file = new TFile("Sim_ab_initio_Barra_NMUONS_300000_PLANES_150x150_RADIO_450_1.root");
// TFile *file = new TFile("Sim_ab_initio_Barra_NMUONS_300000_PLANES_150x150_RADIO_450_0.root");
TTree *tree = (TTree*) file->Get("B02Evts");


int NB = 40;
double tlow = - TMath::Pi() - 0.1;
double thi = TMath::Pi() + 0.1;
TH1F *theta_all = new TH1F("phi_all", "", NB, tlow, thi);
theta_all->GetXaxis()->SetTitle("Angles (Rad)");


TH1F *theta_in = new TH1F("phi_in", "", NB, tlow, thi);
// TH1F *theta_incut = new TH1F("theta_incut", "", NB, tlow, thi);

// Fill histograms //
tree->Draw("phiPri>>phi_all");
tree->Draw("phiPri>>phi_in", " LengthMuLAr >0");
// tree->Draw("phi>>theta_incut", "thet>22");

// Define fuctions //
TF1 *func1 = new TF1("func1", "[0]", 0.01, 2*TMath::Pi() - 0.1);
// func1->SetParameter(0, 7500);

// ======================================= Funcionde ajuste para la CCD  =============================================== ###
TF1 *func2 = new TF1("func2", "[0]*((0.8568933/(2 * TMath::Pi())) + (0.078292082/4)*abs(cos(x)) + (0.064814572/4)*abs(sin(x)))", 0.01,2*TMath::Pi() - 0.05); 

// ======================================= Funcionde ajuste para la barra =============================================== ###
// TF1 *func2 = new TF1("func2", "[0]*((20/TMath::Pi()) + (5)*abs(cos(x)) + (1/2)*abs(sin(x)))",0.01,2*TMath::Pi() - 0.1);


// TF1 *func2 = new TF1("func2", "[0]*((8 * 0.8247739/(2 * TMath::Pi())) + (4 * 0.0499529/4)*(cos(x)) + (4 * 0.0376767715/4)*(sin(x)))", 0.01,2*TMath::Pi() - 0.1); 
// TF1 *func2 = new TF1("func2", "[0]*((8 * 0.8247739/(2 * TMath::Pi())) + (4 * 0.0499529/4)*abs(cos(x)))", 0.01,2*TMath::Pi() - 0.1); // EL mejor
// func2->SetParameter(0, 1600);
//func2->SetParameter(1, 1000);

// Fit functions //
/*theta_all->Fit("func1", "R");
//theta_in->Fit("func2", "R");

double Prob1 = func1->GetProb();
double chi1 = func1->GetChisquare();

double Prob2 = func2->GetProb();
double chi2 = func2->GetChisquare();

std::cout << "ChiSquare1 = "<< chi1 << std::endl;
std::cout << "Prob1 = "<< Prob1 << std::endl;

std::cout << "ChiSquare2 = "<< chi2 << std::endl;
std::cout << "Prob2 = "<< Prob2 << std::endl;
*/
// Create Canvas //
TCanvas *canv = new TCanvas("canv","", 2*700, 600);
canv->Divide(2,1);
canv->cd(1);
theta_all->Draw();
//func1->Draw("same");

canv->cd(2);
theta_in->Draw();
func2->Draw("same"); 


canv->Print("Dis_phi.pdf");

}
