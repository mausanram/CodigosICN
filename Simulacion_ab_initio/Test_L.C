void Test_L(){
TChain *tree = new TChain("tree");
// tree->Add("Sim_ab_initio_NMUONS_1000000_PLANES_1.5_RADIO_8_CCDSIZE_400X600_C.root");
// tree->Add("Sim_ab_initio_NMUONS_2000000_PLANES_1.5_RADIO_8_CCDSIZE_400X600_C_0.root");
// tree->Add("Sim_ab_initio_NMUONS_1000000_PLANES_1.5_RADIO_8_CCDSIZE_400X600_C_1.root");
// tree->Add("Sim_ab_initio_NMUONS_1000000_PLANES_1.5_RADIO_8_CCDSIZE_400X600_C_2.root");
// tree->Add("Sim_ab_initio_NMUONS_1000000_PLANES_1.5_RADIO_8_CCDSIZE_400X600_C_3.root");
// tree->Add("Sim_ab_initio_NMUONS_1000000_PLANES_1.5_RADIO_8_CCDSIZE_400X600_C_4.root");
// tree->Add("Sim_ab_initio_NMUONS_1000000_PLANES_1.5_RADIO_8_CCDSIZE_400X600_C_5.root");
    
// tree->Add("Sim_ab_initio_NMUONS_1000000_PLANES_1.5_RADIO_8_CCDSIZE_250X529_C_0.root");
// tree->Add("Sim_ab_initio_CONNIE_NMUONS_1000000_PLANES_1.7_RADIO_7_CCDSIZE_420X1022_C_0.root");

tree->Add("Sim_ab_initio_NMUONS_1000000_PLANES_1.5x1.5_RADIO_8_CCDSIZE_250x529_SIGMA_LV_1.0_.root");


int NB = 165;
// int NB = 120;
double tlow = 0;

double thi = 0.4;
// double thi = 30;

TH1F *L = new TH1F("L", "Distribucion de Longitudes", NB, tlow, thi);
L->GetXaxis()->SetTitle("Distancia (cm)");
L->SetStats(0);


TH1F *Lcut = new TH1F("Lcut", "", NB, tlow, thi);

// Fill histograms //
tree->Draw("l>>L", "l>0");
tree->Draw("l>>Lcut", "l>0 & thet>22*TMath::Pi()/180");
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

TLine *line = new TLine(0.0725,0,0.0725,12000);
// TLine *line = new TLine(0.068,0,0.068,30000);
line->SetLineStyle(2);
line->SetLineWidth(2);


// Create Canvas //
TCanvas *canv = new TCanvas("canv","", 2*700, 600);
canv->Divide(2,1);
canv->cd(1);
L->Draw();
line->Draw("same");
// Lcut->Draw("same");
// func1->Draw("same");

TLegend *leg = new TLegend(0.7, 0.8, 0.9, 0.9);
// leg->SetTextAlign(11);
leg->SetFillStyle(0);
leg->AddEntry(L, "Datos Simulados", "l");
leg->AddEntry(line, "Grosor de la CCD: 0.0725 cm", "l");
// leg->AddEntry(line, "Grosor de la CCD: 0.068 cm", "l");
leg->Draw();

canv->cd(2);
Lcut->Draw();
// func2->Draw("same");
canv->Print("Dis_L.pdf");

}
