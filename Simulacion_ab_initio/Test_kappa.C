void Test_kappa(){
    // TFile *file = new TFile("Sim_ab_initio_NMUONS_2000000_PLANES_1.5_RADIO_8_CCDSIZE_400X600_C.root");
    TFile *file = new TFile("Sim_ab_initio_NMUONS_1000000_PLANES_1.5_RADIO_8_CCDSIZE_250X529_SIGMA_1.0_C_0.root");


    TTree *tree = (TTree*) file->Get("tree");

    // TFile *file0 = new TFile("Edep_allclusters_NSAMP324_MeV.root");
    // TTree *tree0 = (TTree*) file->Get("tree");

    int NB = 120;
    double tlow = 0;
    double thi = 1; // PAra la CCD

    TH1F *edep_Landau = new TH1F("edep_Landau", "", NB, tlow, thi);
    edep_Landau->GetXaxis()->SetTitle("Energy (MeV)");

    TH1F *edep_Vav = new TH1F("edep_Vav", "", NB, tlow, thi + 3);
    edep_Vav->GetXaxis()->SetTitle("Energy (MeV)");
    // edep->SetGrid(1);

    // Fill histograms //
    tree->Draw("edep>>edep_Landau", "edep>0 & kappa>0 & kappa < 0.01");
    tree->Draw("edep>>edep_Vav", "edep>0 & kappa>0.01 & kappa < 10");

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
    edep_Landau->Draw();
    // edep_cut->Draw("same");
    // func1->Draw("same");

    canv->cd(2);
    edep_Vav->Draw();
    // edep0->Draw("same");
    // func2->Draw("same");
    canv->Print("Dis_kappa.pdf");

}
