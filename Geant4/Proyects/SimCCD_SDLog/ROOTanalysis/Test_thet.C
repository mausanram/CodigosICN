void Test_thet(){
//TFile *file = new TFile("Sim_ab_initio_NMUONS_300000.root");
// TFile *file = new TFile("Sim_ab_initio_NMUONS_400000.root");
// TFile *file0 = new TFile("Sim_ab_initio_NMUONS_100000_PLANES_3.0x3.0_RADIO_12_0.root");
// TFile *file = new TFile("Sim_ab_initio_NMUONS_100000_PLANES_3.0x3.0_RADIO_12_.root");
// TFile *file = new TFile("Sim_ab_initio_NMUONS_50000_PLANES_1x1_RADIO_5_CCDSIZE_400x600_.root");
//TFile *file = new TFile("../build/Sim_NMUONS_100000.root");
TFile *file = new TFile("./root_files/muons_1M_vacuum_250x529_file_m.root");

// TFile *file = new TFile("Sim_ab_initio_Barra_NMUONS_200000_PLANES_150x150_RADIO_100.root");
// TFile *file = new TFile("./treesROOT_Barra/Sim_ab_initio_Barra_NMUONS_300000_PLANES_150x150_RADIO_450_0.root");


// TTree *tree0 = (TTree*) file0->Get("tree");
TTree *tree = (TTree*) file->Get("B02Evts");
// tree->Print();

// int NB = 90;
int NB = 90;
double tlow = 0;
double thi = TMath::Pi()/2.0;
TH1F *theta_all = new TH1F("theta_all", "Distribuci#acute{o}n angular #theta de todos los muones simulados", NB, tlow, thi);
theta_all->GetXaxis()->SetTitle("#acute{A}ngulo (Rad)");
theta_all->SetLineColor(2);
theta_all->SetStats(0);

TH1F *theta_in = new TH1F("theta_in", "Distribuci#acute{o}n angular #theta de los muones que impactaron la CCD", NB, tlow, thi);
theta_in->GetXaxis()->SetTitle("#acute{A}ngulo (Rad)");
theta_in->SetLineColor(2);
theta_in->SetStats(0);

// TH1F *theta_incut = new TH1F("theta_incut", "", NB, tlow, thi);

// Fill histograms //
// tree0->Draw("thet>>theta_all");
tree->Draw("thetaPri>>theta_all");

// tree0->Draw("thet>>theta_in", "l>0 ");
tree->Draw("thetaPri>>theta_in", "nHitBar>0");
// tree->Draw("thet>>theta_incut", "thet>22*TMath::Pi()/180");

// Define fuctions //
TF1 *func1 = new TF1("func1", "[0]*sin(x)*(cos(x))^2", 0.01, 85*TMath::Pi()/180);
func1->SetLineColor(4);
func1->SetParameter(0, 5000);

// ====================== Función para la CCD =========================== //
// TF1 *func2 = new TF1("func2", "[0]*((12.169116/1)*sin(x)*(cos(x))^3+(4.06464/(2*TMath::Pi()))*(sin(x))^2*(cos(x))^2)", 0.01, 85*TMath::Pi()/180);
// TF1 *func2 = new TF1("func2", "[0]*((3.042279/1)*sin(x)*(cos(x))^3+(0.25404/(1*TMath::Pi()))*(sin(x))^2*(cos(x))^2)", 0.01, 85*TMath::Pi()/180);
// TF1 *func2 = new TF1("func2", "[0]*((3.042279/1)*sin(x)*(cos(x))^3+(2*(0.1389825 + 0.1150575)/(1*TMath::Pi()))*(sin(x))^2*(cos(x))^2)", 0.01, 85*TMath::Pi()/180);
// double Ah = 0.54;
// double Al = 0.06525;
// double Ac = 0.0435;

double Ah = 0.2975625;
double Al = 0.05752875;
double Ac = 0.0271875;

TF1 *func2 = new TF1("func2", Form("[0]*((%f)*sin(x)*(cos(x))^3+(2*(%f + %f)/(1*TMath::Pi()))*(sin(x))^2*(cos(x))^2)", Ah, Al, Ac), 0.01, 85*TMath::Pi()/180);
func2->SetLineColor(4);
// ====================== Funcion para la barra ========================  //
// TF1 *func2 = new TF1("func2", "[0]*((10)*sin(x)*(cos(x))^3+(22/(TMath::Pi()))*(sin(x))^2*(cos(x))^2)", 0, 90*TMath::Pi()/180);
// TF1 *func2 = new TF1("func2", "[0]*((10)*sin(x)*(cos(x))^3+(22/4)*(sin(x))^2*(cos(x))^2)", 0, 90*TMath::Pi()/180);
// TF1 *func2 = new TF1("func2", "[0]*((1.25)*sin(x)*(cos(x))^3+(11/(TMath::Pi()))*(sin(x))^2*(cos(x))^2)", 0, 90*TMath::Pi()/180);

func2->SetParameter(0, 2000);

// Fit functions //
theta_all->Fit("func1", "R");
theta_in->Fit("func2", "R");

double Prob1 = func1->GetProb();
double chi1 = func1->GetChisquare();
// double NFD1 = func1->GetN

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
leg->AddEntry(func1, "A_{0}sin#theta cos^{2}#theta", "l");
leg->AddEntry(theta_all, "Datos Simulados", "f");
leg->Draw();

canv->cd(2);
theta_in->Draw();
func2->Draw("same");

leg = new TLegend(0.7, 0.8, 0.9, 0.9);
// leg->SetTextAlign(11);
leg->SetHeader("", "C");
leg->SetFillStyle(0);
leg->AddEntry(func1, "A_{1} [(0.2975625)sin#theta cos^{3}#theta", "l");
leg->AddEntry((TObject*)0, "", "");
leg->AddEntry((TObject*)0, "+ #left(#frac{0.1694325}{#pi}#right)sin^{2}#theta cos^{2}#theta]", " ");
leg->AddEntry((TObject*)0, "", "");
leg->AddEntry(theta_all, "Datos Simulados", "f");
leg->Draw();




canv->Print("Dis_thet.pdf");

}
