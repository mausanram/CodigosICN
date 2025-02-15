void Test_edep(){
    // TFile *file = new TFile("Sim_ab_initio_NMUONS_10000_PLANES_1.5x1.5_RADIO_8_CCDSIZE_400x600_SIGMA_LV_1.0_DELTAPDG_.root");
    // TFile *file = new TFile("Sim_ab_initio_NMUONS_10000_PLANES_1.5x1.5_RADIO_8_CCDSIZE_400x600_SIGMA_LV_1.0_DELTALEO_.root");
    // TFile *file = new TFile("Sim_ab_initio_NMUONS_300000_PLANES_1.5x1.5_RADIO_8_CCDSIZE_400x600_SIGMA_LV_1.0_.root");
    // TFile *file = new TFile("Sim_ab_initio_NMUONS_1000000_PLANES_1.5x1.5_RADIO_8_CCDSIZE_400x600_SIGMA_LV_1.0_.root");

    TFile *file = new TFile("Sim_ab_initio_NMUONS_1000000_PLANES_1.5_RADIO_8_CCDSIZE_400X600_C.root");
    // TFile *file = new TFile("Sim_ab_initio_Barra_NMUONS_1000000_PLANES_150x150_RADIO_450.root");


    TTree *tree = (TTree*) file->Get("tree");

    // TFile *file0 = new TFile("Edep_allclusters_NSAMP324_MeV.root");
    // TTree *tree0 = (TTree*) file->Get("tree");

    int NB = 120;
    double tlow = 0;
    double thi = 1; // PAra la CCD
    TH1F *edep = new TH1F("edep", "", NB, tlow, thi);
    edep->GetXaxis()->SetTitle("Energy (MeV)");
    // edep->SetGrid(1);


    TH1F *edep_cut = new TH1F("edep_cut", "", NB, tlow, thi);

    // TH1F *edep0 = new TH1F("edep0", "", NB, tlow, thi);

    // Fill histograms //
    tree->Draw("edep>>edep", "edep > 0");
    tree->Draw("edep>>edep_cut", "thet>22*TMath::Pi()/180 & edep>0");

    // tree0->Draw("edep>>edep0");

    // // Define fuctions //
    // TF1 *func1 = new TF1("func1", "[0]*sin(x)*(cos(x))^2", 0, 90);
    // func1->SetParameter(0, 5000);

    // TF1 *func2 = new TF1("func2", "[0]*sin(x)*(cos(x))^3+[1]*(sin(x))^2*(cos(x))^2", 0, 90);
    // func2->SetParameter(0, 2000);
    // func2->SetParameter(1, 1000);

    // Fit functions //
    // theta_all->Fit("func1", "R");
    // theta_in->Fit("func2", "R");

    // Create Canvas //
    TCanvas *canv = new TCanvas("canv","", 2*700, 600);
    canv->Divide(2,1);
    canv->cd(1);
    edep->Draw();
    // edep_cut->Draw("same");
    // func1->Draw("same");

    canv->cd(2);
    edep_cut->Draw();
    // edep0->Draw("same");
    // func2->Draw("same");
    canv->Print("Dis_edep.pdf");

}
