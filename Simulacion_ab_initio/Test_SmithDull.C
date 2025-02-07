
double Smith_Dull(double *lx , double *lpar){

    double E_mu = lx[0];
    double theta = TMath::Pi() * lpar[0]/180;

    //## ------------- Constantes físicas -------------- ##
    double k = 8.0 / 3.0;
    double b = 0.771;
    double lambda_pi = 120;    // g/cm^2
    double y_0 = 1000;      // g/cm^2
    double r = 0.76;
    double b_mu = 0.8;   

    // # if units == 0:
    double a = 2.5;      //## MeV cm^2/g
    double m_mu = 105.7;    //## MeV/c^2 
    double m_pi = 139.6;    //## MeV/c^2

    double tau_mu_0 = 2.2 * pow(10,-6);   //## s
    double tau_0 = 2.6 * pow(10, -8);     //## s
    double rho_0 = 0.00129; //## g/cm^3
    double c = TMath::C()* 100; //## cm/s

    //### ---------------------- Parámetros ---------------------- ###
    double E_pi = (1 / r) * (E_mu + a * y_0 * ((1/TMath::Cos(theta)) - 0.1));
    double B_mu = (b_mu * m_mu * y_0)/(tau_mu_0 * rho_0 * c);
    // double B_mu = (b_mu * m_mu * y_0 * c)/(tau_mu_0 * rho_0);
    double P_mu = pow((0.1 * TMath::Cos(theta)) * (1 - (a * (y_0 *(1/TMath::Cos(theta)) - 100))/( r * E_pi)), ((B_mu)/((r * E_pi + 100 * a) * TMath::Cos(theta))));
    double j_pi = (m_pi * y_0)/(c * tau_0 * rho_0);
    // double j_pi = (m_pi * y_0 * c)/(tau_0 * rho_0);

    //## Intensidad diferencial
    //# C_1 = E_pi ** (-k) * P_mu * lambda_pi * b * j_pi
    double C_1 = pow(E_pi, -k)* P_mu * lambda_pi * b * j_pi;
    double C_2 = E_pi * TMath::Cos(theta);
    double C_3 = b * j_pi;

    //# return (C_1 * np.sin(theta)) / (C_2 + C_3)
    double En_k = (C_1 ) / (C_2 + C_3);
    return En_k;
    }
/// ============================================================= ///


void Test_SmithDull(){
    cout << "Introduce el ángulo Theta del muon: (en grados): ";	// ---------------------------------------- //
	double p;						// Esta sección es para pedir que se ingrese el momento de los muones a mano 
   	std::cin >> p;	

    double thet = p * TMath::Pi()/180;


    double nbins = 100;
	double hlow = 1;
	double hhi = pow(10, 5);


	TF1 *f = new TF1("f", Smith_Dull, 1, pow(10, 5), 1);
	f->SetParameter(0,thet );

    TH1F *hist = new TH1F("hist", "", nbins, hlow, hhi);
	hist->FillRandom("f", 100000);
    hist->Scale(0.0000000000001);


    int NB = 100;
    double tlow = 1;
    double thi = pow(10,5);
    double xbins[NB+1];
    double low_bin;

    double Delta_log = (log10(thi) - log10(tlow))/(NB -1);
    for (int i = 0; i<NB;i++){
        low_bin = pow(10, log10(tlow) + i * Delta_log);
        cout<<low_bin <<endl;
        xbins[i] = low_bin;
    };

    xbins[NB] = thi;

    TH1F *hist2 = new TH1F("hist", "", NB, xbins);
	hist2->FillRandom("f", 100000);
    hist2->Scale(0.000000000002);

	// f->SetTitle("Landau-Vavilov distribution (for 0.0725cm of Si);#font[12]{Energy} (MeV);Probability");
	// f->SetTitle("Landau-Vavilov distribution (0.0725cm of Si);Energy (MeV);Probability");

    TCanvas *cnv = new TCanvas("cnv", "", 2 *   900, 700);
	cnv->SetGrid();
    cnv->Divide(2,1);

    cnv->cd(1);
	f->Draw();
    hist->Draw("h0 same");

    gPad->SetLogy(1);
    gPad->SetLogx(1);

    cnv->cd(2);
    f->Draw();
    hist2->Draw("h0 same");
    // hist->Draw();
    // f->Draw("same");

    gPad->SetLogy(1);
    gPad->SetLogx(1);

    // ==================================================== %%

    

    // hist->Draw();
    // f->Draw("same");

    gPad->SetLogy(1);
    gPad->SetLogx(1);

    cnv->Print("Dis_Smith_Duller.pdf");
}