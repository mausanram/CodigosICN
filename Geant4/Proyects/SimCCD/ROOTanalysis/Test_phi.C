void Test_phi(){
//TFile *file = new TFile("Sim_ab_initio_NMUONS_300000.root");
// TFile *file = new TFile("Sim_ab_initio_NMUONS_400000.root");
// TFile *file = new TFile("Sim_ab_initio_NMUONS_10000_PLANES_3.0x3.0_RADIO_12_.root");
TFile *file = new TFile("./root_files/muons_1M_vacuum_250x529_file_m.root");

// TFile *file = new TFile("Sim_ab_initio_Barra_NMUONS_200000_PLANES_150x150_RADIO_100.root");
// TFile *file = new TFile("Sim_ab_initio_Barra_NMUONS_500000_PLANES_150x150_RADIO_100(1).root");
// TFile *file = new TFile("Sim_ab_initio_Barra_NMUONS_100000_PLANES_150x150_RADIO_450_0.root");
// TFile *file = new TFile("Sim_ab_initio_Barra_NMUONS_300000_PLANES_150x150_RADIO_450_1.root");
// TFile *file = new TFile("Sim_ab_initio_Barra_NMUONS_300000_PLANES_150x150_RADIO_450_0.root");
TTree *tree = (TTree*) file->Get("B02Evts");


int NB = 50;
// double tlow = - TMath::Pi() - 0.1;
double tlow = 0;
double thi = 2 * TMath::Pi() + 0.1;
TH1F *theta_all = new TH1F("phi_all", "Distribuci#acute{o}n angular #phi de todos los muones simulados", NB, tlow, thi);
theta_all->GetXaxis()->SetTitle("#acute{A}ngulos (Rad)");
theta_all->SetStats(0);
theta_all->SetLineColor(2);


TH1F *theta_in = new TH1F("phi_in", "Distribuci#acute{o}n angular #phi de los muones que impactaron la CCD", NB, tlow, thi);
theta_in->SetStats(0);
theta_in->GetXaxis()->SetTitle("#acute{A}ngulos (Rad)");
theta_in->SetLineColor(2);
// TH1F *theta_incut = new TH1F("theta_incut", "", NB, tlow, thi);

// Fill histograms //
tree->Draw("(phiPri + TMath::Pi())>>phi_all");
tree->Draw("(phiPri + TMath::Pi())>>phi_in", " EevtBar > 0 ");
// tree->Draw("phi>>theta_incut", "thet>22");

// Define fuctions //
TF1 *func1 = new TF1("func1", "[0]", 0.01, 2*TMath::Pi() - 0.1);
func1->SetLineColor(4);
// func1->SetParameter(0, 7500);

double Ah = 0.2975625;
double Al = 0.05752875;
double Ac = 0.0271875;

// ======================================= Funcionde ajuste para la CCD  =============================================== ###
// TF1 *func2 = new TF1("func2", Form("[0]*(((%f)/(1 * TMath::Pi())) + ((%f)/2)*abs(cos(x)) + ((%f)/2)*abs(sin(x)))", Ah, Al, Ac), 0.05,2*TMath::Pi() - 0.1); 
TF1 *func2 = new TF1("func2", Form("[0]*(((%f)/(1 * TMath::Pi())) + ((%f)/4)*abs(cos(x)) + ((%f)/4)*abs(sin(x)))", Ah, Al, Ac), 0.05,2*TMath::Pi() - 0.1); 
func2->SetLineColor(4);
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

TLegend *leg = new TLegend(0.5, 0.8, 0.9, 0.9);
// leg->SetTextAlign(11);
leg->AddEntry(func1, "A_{2}", "l");
leg->AddEntry(theta_all, "Datos Simulados", "f");
leg->Draw();


canv->cd(2);
theta_in->SetMinimum(0);
theta_in->Draw();
func2->Draw("same"); 

leg = new TLegend(0.7, 0.8, 0.9, 0.9);
// leg->SetTextAlign(11);
leg->SetHeader("", "C");
// leg->SetFillStyle(0);
leg->AddEntry(func1, "A_{3} [0.2975625 + #left(#frac{0.05752875}{4}#right)#left|cos#theta #right| ", "l");
leg->AddEntry((TObject*)0, "", "");
leg->AddEntry((TObject*)0, "+ #left(#frac{0.0271875}{4}#right) #left|sin#theta #right| ]", " ");
leg->AddEntry((TObject*)0, "", "");
leg->AddEntry(theta_all, "Datos Simulados", "f");
leg->Draw();

canv->Print("Dis_phi.pdf");

}
