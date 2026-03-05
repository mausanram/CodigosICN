void Test_edep_Chain(){

    TChain *chain = new TChain("B02Evts");
    chain->Add("./root_files/muons_1M_vacuum_250x529_file_SDLog_Cu_*.root");


    // Sección de Energía depositada //
    int NB = 150;
    double tlow = 0;
    double thi = 1000;

    TH1F *edep = new TH1F("edep", "", NB, tlow, thi); // Histogram simulation all
    edep->GetXaxis()->SetTitle("Energy (KeV)");
    // edep->Scale(5);

    TH1F *edep_cut = new TH1F("edep_cut", "Sigma: 1", NB, tlow, thi); // Histogram simulation cut 
    edep_cut->GetXaxis()->SetTitle("Energy (KeV)");
    edep_cut->SetLineColor(1);

    TH1F *edep_exp = new TH1F("edep_exp", "", NB, tlow, thi); // Histogram simulation cut 
    edep_exp->GetXaxis()->SetTitle("Energy (KeV)");
    // edep_exp->Scale(5);

    // chain->Draw("EevtBar*1000>>edep", "EevtBar > 0");
    chain->Draw("EevtBar*1000>>edep", "nHitBar>0 && LengthMuLAr>0");
    chain->Draw("EevtBar*1000>>edep_cut", "nHitBar>0 && LengthMuLAr>0 && thetaPri>22*TMath::Pi()/180");
    // chain->Draw("edep>>edep_exp", "l<0");


    //scale hint1 to the pad coordinates
    double rightmax = 0.00168*edep->GetMaximum();
    // double scale = gPad->GetUymax()/rightmax;
    // hint1->SetLineColor(kRed);
    // edep->Scale(scale);
    edep->SetLineColor(2);
    // edep->SetLineStyle(0);
    // edep->SetFillStyle(0);
    // hint1->Draw("same");


    rightmax = 0.00168*edep_cut->GetMaximum();
    // scale = gPad->GetUymax()/rightmax;
    // hint1->SetLineColor(kRed);
    // edep_cut->Scale(scale);
    // hint1->Draw("same");

    TCanvas *canv = new TCanvas("canv","", 2*700, 600);
    canv->Divide(2,1);
    canv->cd(1);
    edep->Draw("hist");
    // edep_exp->Draw("hist same");

    canv->cd(2);
    edep_cut->Draw("hist");
    // edep_exp->Draw("same");


    canv->Print("Dis_edep_chain_SIGMA_1.pdf");

}
