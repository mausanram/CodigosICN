void Test_phi(){

TChain *tree = new TChain("tree");
// tree->Add("Sim_ab_initio_NMUONS_1000000_PLANES_1.5_RADIO_8_CCDSIZE_400X600_C.root");
// tree->Add("Sim_ab_initio_NMUONS_2000000_PLANES_1.5_RADIO_8_CCDSIZE_400X600_C_0.root");
// tree->Add("Sim_ab_initio_NMUONS_1000000_PLANES_1.5_RADIO_8_CCDSIZE_400X600_C_0.root");
// tree->Add("Sim_ab_initio_NMUONS_1000000_PLANES_1.5_RADIO_8_CCDSIZE_400X600_C_1.root");
// tree->Add("Sim_ab_initio_NMUONS_1000000_PLANES_1.5_RADIO_8_CCDSIZE_400X600_C_2.root");
// tree->Add("Sim_ab_initio_NMUONS_1000000_PLANES_1.5_RADIO_8_CCDSIZE_400X600_C_3.root");
// tree->Add("Sim_ab_initio_NMUONS_1000000_PLANES_1.5_RADIO_8_CCDSIZE_400X600_C_4.root");
// tree->Add("Sim_ab_initio_NMUONS_1000000_PLANES_1.5_RADIO_8_CCDSIZE_400X600_C_5.root");

// tree->Add("Sim_ab_initio_NMUONS_1000000_PLANES_1.5_RADIO_8_CCDSIZE_250X529_C_0.root");
// tree->Add("Sim_ab_initio_CONNIE_NMUONS_1000000_PLANES_1.7_RADIO_7_CCDSIZE_420X1022_C_0.root");

tree->Add("Sim_ab_initio_NMUONS_1000000_PLANES_1.5x1.5_RADIO_8_CCDSIZE_250x529_SIGMA_LV_1.0_.root");


int NB = 70;
double tlow = 0;
double thi = 2*TMath::Pi() + 0.1;
TH1F *theta_all = new TH1F("phi_all", "Distribuci#acute{o}n angular #phi de todos los muones simulados", NB, tlow, thi);
theta_all->GetXaxis()->SetTitle("#acute{A}ngulo (Rad)");
theta_all->SetStats(0);


TH1F *theta_in = new TH1F("phi_in", "Distribuci#acute{o}n angular #phi de los muones que impactaron la CCD", NB, tlow, thi);
theta_in->GetXaxis()->SetTitle("#acute{A}ngulo (Rad)");
theta_in->SetStats(0);

// TH1F *theta_incut = new TH1F("theta_incut", "", NB, tlow, thi);

// Fill histograms //
tree->Draw("phi>>phi_all");
tree->Draw("phi>>phi_in", "edep > 0");
// tree->Draw("phi>>theta_incut", "thet>22");

// Define fuctions //
TF1 *func1 = new TF1("func1", "[0]", 0.1, 2*TMath::Pi() - 0.1);
func1->SetParameter(0, 100000);

// double Ah = 0.54;
// double Al = 0.06525;
// double Ac = 0.0435;

double Ah = 0.2975625;
double Al = 0.05752875;
double Ac = 0.0271875;

// === CONNIE === //
// double Ah = 0.96579;
// double Al = 0.104244;
// double Ac = 0.04284;

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
leg->SetFillStyle(0);
leg->AddEntry(func1, "A_{3} [0.96579 + #left(#frac{ 0.104244}{4}#right)#left|cos#theta #right| ", "l");
leg->AddEntry((TObject*)0, "", "");
leg->AddEntry((TObject*)0, "+ #left(#frac{0.04284}{4}#right) #left|sin#theta #right| ]", " ");
leg->AddEntry((TObject*)0, "", "");
leg->AddEntry(theta_all, "Datos Simulados", "f");
leg->Draw();

canv->Print("Dis_phi.pdf");

}
