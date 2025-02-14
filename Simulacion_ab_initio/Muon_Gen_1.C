
# include <chrono>

double dis_thet(double *lx, double *lpar){ //## Distribucion dis_angular
    double theta = lx[0];
    return (pow(TMath::Cos(theta), 2)) * (TMath::Sin(theta));
};

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

double Smith_Dull_Log(double *lx , double *lpar){

    double E_mu = pow(10, lx[0]);
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
    double log_ek = En_k;

    return log_ek;
    // return log10(log_ek);
    // return pow(10, log_ek);
    // return pow(10, log10(log_ek));
    }
    
double intersection_CCD(bool *flags_CCD, double *list_z, double medida_z, double Random_th){
    //delta_L = None
    double delta_L = 0;
    double h;

    //## Caras 1 y 2
    if (flags_CCD[0] & flags_CCD[1]){
        delta_L = (2 * medida_z) / TMath::Cos(Random_th);  //## cm
    }

    //## Caras 1 y 3
    if (flags_CCD[0] & flags_CCD[2]){
        h = medida_z - list_z[0];
        delta_L = h /  TMath::Cos(Random_th);   //## cm
    }
    //## Caras 1 y 4
    if (flags_CCD[0] & flags_CCD[3]){
        h = medida_z - list_z[1];
        delta_L =  (h /   TMath::Cos(Random_th));  //## cm
    }
    // ## Caras 1 y 5
    if (flags_CCD[0] & flags_CCD[4]){
        h = medida_z - list_z[2];
        delta_L =  (h /   TMath::Cos(Random_th));  //## cm
    }
    //## Caras 1 y 6
    if (flags_CCD[0] & flags_CCD[5]){
        h = medida_z - list_z[3];
        delta_L =  (h /   TMath::Cos(Random_th));  //## cm
    }

    // ## Caras 3 y 2
    if (flags_CCD[2] and flags_CCD[1]){
        h = medida_z + list_z[0]; 
        delta_L =  (h / TMath::Cos(Random_th)); //## cm
    }

    // ## Caras 3 y 4
    if ((list_z[0] > list_z[1]) & flags_CCD[2] & flags_CCD[3]){
        h = list_z[0] - list_z[1];
        delta_L =  (h /  TMath::Cos(Random_th));  //## cm
    }

    // ## Caras 3 y 5
    if ((list_z[0] > list_z[2]) & flags_CCD[2] & flags_CCD[4]){
        h = list_z[0] - list_z[2];
        delta_L =  (h /  TMath::Cos(Random_th));  //## cm
    }

    // ## Caras 3 y 6
    if ((list_z[0] > list_z[3]) & flags_CCD[2] & flags_CCD[5]){
        h = list_z[0] - list_z[3];
        delta_L =  (h/  TMath::Cos(Random_th));  //## cm
    }

    // ## Caras 4 y 2
    if (flags_CCD[3] and flags_CCD[1]){
        h = medida_z + list_z[1];
        delta_L =  (h /  TMath::Cos(Random_th));  //## cm
    }

    // ## Caras 4 y 3
    if ((list_z[1] > list_z[0]) & flags_CCD[3] & flags_CCD[2]){
        h = list_z[1] - list_z[0];
        delta_L =  (h /  TMath::Cos(Random_th));  //## cm
    }

    // ## Caras 4 y 5
    if ((list_z[1] > list_z[2]) & flags_CCD[3] & flags_CCD[4]){
        h = list_z[1] - list_z[2];
        delta_L =  (h /  TMath::Cos(Random_th));  //## cm
    }

    // ## Caras 4 y 6
    if ((list_z[1] > list_z[3]) & flags_CCD[3] & flags_CCD[5]){
        h = list_z[1] - list_z[3];
        delta_L =  (h /  TMath::Cos(Random_th));  //## cm
    }

    // ## Caras 5 y 2
    if (flags_CCD[4] & flags_CCD[1]){
        h = list_z[2] + medida_z;
        delta_L =  (h /  TMath::Cos(Random_th));  //## cm
    }

    // ## Caras 5 y 3
    if ((list_z[2] > list_z[0]) & flags_CCD[4] & flags_CCD[2]){
        h = list_z[2] - list_z[0];
        delta_L =  (h /  TMath::Cos(Random_th));  //## cm
    }

    // ## Caras 5 y 4
    if ((list_z[2] > list_z[1]) & flags_CCD[4] & flags_CCD[3]){
        h = list_z[2] - list_z[1];
        delta_L =  (h /  TMath::Cos(Random_th));  //## cm
    }

    // ## Caras 5 y 6
    if ((list_z[2] > list_z[3]) & flags_CCD[4] & flags_CCD[5]){
        h = list_z[2]  - list_z[3];
        delta_L =  (h /  TMath::Cos(Random_th));  //## cm
    }

    // ## Caras 6 y 2
    if (flags_CCD[5] and flags_CCD[1]){
        h = list_z[3] + medida_z;
        delta_L =  (h /TMath::Cos(Random_th)); //## cm
    }

    // ## Caras 6 y 3
    if ((list_z[3] > list_z[0]) & flags_CCD[5] & flags_CCD[2]){
        h = list_z[3] - list_z[0];
        delta_L =  (h /  TMath::Cos(Random_th));  //## cm
    }

    // ## Caras 6 y 4
    if ((list_z[3] > list_z[1]) & flags_CCD[5] & flags_CCD[3]){
        h = list_z[3] - list_z[1];
        delta_L =  (h /  TMath::Cos(Random_th));  //## cm
    }

    // ## Caras 6 y 5
    if ((list_z[3] > list_z[2]) & flags_CCD[5] & flags_CCD[4]){
        h = list_z[3] - list_z[2];
        delta_L =  (h /  TMath::Cos(Random_th));  //## cm
    }

    return delta_L ;
};

double LV (double *lx, double *lpar) {
	double Delta = lx[0];	// Energy loss in absorber
	double L = lpar[0];		// Thickness of absorber (Distance crossed by the particle)
	double p = lpar[1];		// Momentum (in MeV/c)


	double z = 1.0;		// Charge number of incident particle
	double ZA = 0.498487;	// Atomic number over Atomic mass of absorber (for Si)
	double c = TMath::C();	// Speed of light
	double me = 0.510998928;	// Electron mass in MeV/c^2
	double M = 105.65839;	// Muon mass in MeV/c^2
	double I = 0.000173;		// Mean excitation energy (for Si)

	double bg = p/M;
	double beta = bg/sqrt(1+(pow(bg,2)));	// Beta factor
	double gamma = 1/sqrt(1-(pow(beta,2)));	// Gamma factor
	double pi = TMath::Pi();
	double rho = 2.33;	// Density of material (for Si) gcm^-3

	double d;	// Variable for the Density effect

	double K = 0.1535;	// K coefficient = 2*pi*N*r^2*m*c^2 (in MeV mol^-1 cm^2)
	// double K = 0.307075;	// K coefficient = 4*pi*N*r^2*m*c^2 (in MeV mol^-1 cm^2)


	double a = 0.1492;	// Parameters (taken from W.R. Leo for SI)
	double k = 3.2546;		//
	double X0 = 0.2014;		//
	double X1 = 2.8715;		//
	double C = 4.4351;		//
	double d0 = 0.14;		//
	double X = log10(bg);

	if (X>=X1) {
		d = 2*log(10.0)*X-C;
		}
	else if (X0<=X && X<X1) {
		// cout<<"Entro aquí" << endl;
		d = 2*log(10.0)*X-C+a*(pow((X1-X),k));
		}
	else if (X<X0) {
		d = d0*(pow(10,(2*(X-X0))));
		// d = 0;
		}

	// double WM = 2*me*(pow((beta*gamma),2))/(1+(2*me/M)*(sqrt(1 + pow(beta*gamma, 2)))+pow((me/M),2));   // Maximum energy tranfer
	double WM = (2*me*pow(bg,2))/(1+(2*me*gamma/M) + pow((me/M),2));   // Maximum energy tranfer

	// double loge = log((1-(pow(beta,2)))*(pow(I,2))/(2*M*(pow(beta,2))))+(pow(beta,2)); // log epsilon variable
	double loge = log(((1-pow(beta,2))*(pow(I,2)))/(2*me*pow(beta,2)))+(pow(beta,2)); // log epsilon variable

	double EC = 0.577;	// Euler's constant

    double xi = K * rho * ZA * L * pow((1/beta),2);		// Xi variable 


	double DeltaAv = xi * (log((2*me*(pow(bg,2)))*WM/(pow(I,2)))-(2 * pow(beta,2))-(d));// Mean energy loss (Bethe-Bloch)
	// double DeltaAv = K*rho*L*(pow(z,2))*ZA*(1.0/(pow(beta,2)))*((1/2)* log(2*me*(pow(gamma,2))*(pow(beta,2))*WM/(pow(I,2)))-(pow(beta,2))-(d/2.0));// Mean energy loss (Bethe-Bloch)

	// ===================================================== //

	double Log_1 = log((2 * me * pow(bg,2))/(I));
	double Log_2 = log((xi)/(I));
	double SumLogs = Log_1 + Log_2;
	double SumLogsj = Log_1 + Log_2 + 0.2;
	double SumLogsjbet = Log_1 + Log_2 + 0.2 - (pow(beta,2));
	double SumTerms = Log_1 + Log_2 + 0.2 - (pow(beta,2)) - d;
	double MostP = xi * (Log_1 + Log_2 + 0.2 - pow(beta,2) - d);

	// cout << "Csi Value " << xi << " MeV" << endl;
	// cout << "1st Log: " << Log_1 << endl;
	// cout << "2nd Log: " << Log_2 << endl;
	// cout << "density corr: " << d << endl;
	// cout << "beta^2: " << pow(beta,2) << endl;
	// cout << "SumLogs: " << SumLogs << endl;
	// cout << "SumLogsj: " << SumLogsj << endl;
	// cout << "SumLogsjbet: " << SumLogsjbet << "\n " <<endl;
	// cout << "SumT: " << SumTerms << endl;
	// cout << "MostP: " << MostP << endl;
	// double xi = (K/2)*rho*ZA*L*pow((z/beta),2);		// Xi variable 

	// ==================================================== //

	double lambda = (Delta-xi*(log(xi)-loge+1-EC))/xi;	// Lambda parameter

  	// double Deltamp = xi * (log(xi) - loge + 0.198-d);		// Most probable energy loss
    double Deltamp = xi * (log((2 * me * pow(bg,2))/I) + log(xi/I) + 0.2 - pow(beta,2) - d);		//# Most probable energy loss from PDG
  	double lambdamp = (Deltamp - xi*(log(xi)-loge+1-EC))/xi;

	double kappa = xi/WM;		// Kappa ratio
	double beta2 = pow(beta,2);

	double sigma2 = (pow(xi,2))*(1-beta2/2)/kappa;		// Standard deviation for relativistic particles

	//  cout << lambda << " " << Delta <<"  "<< phi/csi << endl;

	//printf(Deltamp);
	// cout << "Most Probably Energy in KeV " << Deltamp * 1000 << endl;
	//td::cout << Deltamp * 1000 << std::endl;
	// std::cout << "DMP "<<  Deltamp * 1000 << std::endl;
	// std::cout << "MP "<<  lambdamp << std::endl;
	// std::cout << "Kappa: "<<  kappa << std::endl;
	// cout << "LambdaMP: " << lambdamp << " , Deltamp: " << Deltamp << ", DeltaAv: " << DeltaAv << endl;

	if (kappa<=0.01) {
		// std::cout << "lambMP "<<  lambdamp << std::endl;
		// std::cout << "DeltaMP "<<  Deltamp << std::endl;
		double phi = TMath::Landau(lambda, lambdamp, 1.0);
        // double phi = TMath::Landau(lambda, Deltamp, 1.0);
		return phi/xi;
		// return phi;
		}
	else if (0.01<kappa && kappa<10) {
		// std::cout << "2"<< std::endl;
		double vav = TMath::Vavilov(lambda-lambdamp, kappa, beta2);
        // double vav = TMath::Vavilov(Delta-Deltamp, kappa, beta2);
		return vav;
		}
	else {
//		double gauss = exp(((Delta-DeltaAv)**2)/(2*sigma2));
		double gauss =  TMath::Gaus(Delta, DeltaAv, sqrt(sigma2));
		return gauss;
		}
}

double Kappa(double Delta_l, double momentum){
    double z = 1.0;		// Charge number of incident particle
	double ZA = 0.498487;	// Atomic number over Atomic mass of absorber (for Si)
	double c = TMath::C();	// Speed of light
	double me = 0.510998928;	// Electron mass in MeV/c^2
	double M = 105.65839;	// Muon mass in MeV/c^2
	double I = 0.000173;		// Mean excitation energy (for Si)
    double bg = momentum/M;
	double beta = bg/sqrt(1+(pow(bg,2)));	// Beta factor
	double gamma = 1/sqrt(1-(pow(beta,2)));	// Gamma factor
	double pi = TMath::Pi();
	double rho = 2.33;	// Density of material (for Si) gcm^-3
    double K = 0.1535;	// K coefficient = 2*pi*N*r^2*m*c^2 (in MeV mol^-1 cm^2)

    double WM = 2 * me * (pow((beta * gamma),2))/(1+(2 * me * gamma / M) + pow((me / M),2));   // Maximum energy tranfer
	double xi = (K) * rho * ZA * Delta_l * pow((z/beta),2);		// Xi variable 

    double kappa = xi/WM;		// Kappa ratio
    return kappa;
}


/// ============================================================= ///

void Muon_Gen_1(){

    auto start = chrono::high_resolution_clock::now();

    int number_thet = 2000000;

    double Radio = 8; //cm
    double half_plane_size = 0.75; // cm 

    double sizex_pixels = 400; // cm
    double sizey_pixels = 600; // cm
    double pixel_size =  0.0015; //cm

    double medida_x = (sizex_pixels * pixel_size) / 2;   //# cm
    double medida_y = (sizey_pixels * pixel_size)  / 2;   //# cm
    double medida_z = 0.0725 / 2;  //# cm

    // double list_thet_in_CCD[number_thet];
    // double list_phi_in_CCD[number_thet];
    // double list_energy_pri_in_CCD[number_thet];
    // double list_delta_L[number_thet];
    // double list_edep[number_thet];
    // double list_kappa[number_thet];
    // double list_muonID[number_thet];

    TRandom3 rand;

    // int Len_th = sizeof(list_thet_in_CCD);
    // int Len_phi = sizeof(list_phi_in_CCD);
    // int Len_epri = sizeof(list_energy_pri_in_CCD);


    // cout<< "Len thet list: " <<Len_th << endl;
    // cout<< "Len phi list: " <<Len_phi << endl;
    // cout<< "Len epri list: " <<Len_epri << endl;

    TF1 *f_thet = new TF1("", dis_thet, 0, TMath::Pi()/2, 0);
    TF1 *f_SD = new TF1("", Smith_Dull, 1, pow(10,5), 1);
    TF1 *f_LandVav = new TF1("", LV, 0.001, 10, 2);

    int muon_inbucle = 0;

    // Sim_ab_initio_NMUONS_200000_PLANES_1.5x1.5_RADIO_8_CCDSIZE_400x600_SIGMA_LV_1.0_.root

    TFile *file = TFile::Open("Sim_ab_initio_NMUONS_2000000_PLANES_1.5_RADIO_8_CCDSIZE_400X600_C.root", "recreate");
    TTree *tree = new TTree("tree", "tree");

    double Rand_thet;
    double Rand_phi;
    double Rand_epri;
    double Delta_L;
    double Rand_edep;
    double kappa;
    double MuonID;

    double vec[3];
    double Point[3];
    double vec_thet[3];
    double vec_phi[3];
    double P_vector[3];
    double random_plane_point[3];
    double normal_Vec[3];
    bool list_flags[6];
    double list_z[4];

    tree->Branch("muonID", &MuonID);
    tree->Branch("thet", &Rand_thet);
    tree->Branch("phi", &Rand_phi);
    tree->Branch("epri", &Rand_epri);
    tree->Branch("l", &Delta_L);
    tree->Branch("edep", &Rand_edep);
    tree->Branch("kappa", &kappa);

    double max_long = sqrt(pow((medida_x * 2),2) + (pow((medida_y * 2),2) + (pow((medida_z * 2),2))));
    int muons_in_CCD = 0;

    for (int i = 0; i<number_thet;i++){
        gRandom->SetSeed(0);

        //### ================== Seleccion aleatoria de theta, phi y en_pri(Smith-Duller) =============== ###
        Rand_thet = f_thet->GetRandom(); // Rad
        Rand_phi = rand.Rndm() * (2 * TMath::Pi()); // Rad

        f_SD->SetParameter(0, Rand_thet);
        // Rand_epri = pow(10, f_SD->GetRandom()); // MeV
        Rand_epri = f_SD->GetRandom(); // MeV

        // tree->Fill();
        // =============================================================================================== //

        // ### ==================== Momento del muon ======================= ###
        double m_mu = 105.7;    //## MeV/c^2
        double En_tot = Rand_epri + m_mu;  
        double momentum = sqrt(pow(En_tot,2) - pow(m_mu, 2));
        // ### ============================================================= ###
        
        // ============ Vector de dirección sobre la esfera ========== //
        vec[0] = TMath::Sin(Rand_thet) * TMath::Cos(Rand_phi);
        vec[1] = TMath::Sin(Rand_thet) * TMath::Sin(Rand_phi);
        vec[2] = TMath::Cos(Rand_thet);

        Point[0] = Radio * vec[0]; 
        Point[1] = Radio * vec[1]; 
        Point[2] = Radio * vec[2];  //## Genera un punto sobre la esfera

        vec_thet[0] = TMath::Cos(Rand_thet) * TMath::Cos(Rand_phi); 
        vec_thet[1] = TMath::Cos(Rand_thet) * TMath::Sin(Rand_phi);
        vec_thet[2] = -1 * TMath::Sin(Rand_thet);

        vec_phi[0] = -1 * TMath::Sin(Rand_phi);
        vec_phi[1] = TMath::Cos(Rand_phi); 
        vec_phi[2] = 0;

        double random_a = -half_plane_size + rand.Rndm() * 2 * half_plane_size;  // ## Selecciona un valor uniforme para el parámetro a
        double random_b = -half_plane_size + rand.Rndm() * 2 * half_plane_size; // ##      ''      ''      ''      ''          ''    b
        
        P_vector[0] = random_a * vec_thet[0] + random_b * vec_phi[0]; 
        P_vector[1] = random_a * vec_thet[1] + random_b * vec_phi[1]; 
        P_vector[2] = random_a * vec_thet[2] + random_b * vec_phi[2];

        random_plane_point[0] = Point[0] + P_vector[0];
        random_plane_point[1] = Point[1] + P_vector[1];
        random_plane_point[2] = Point[2] + P_vector[2];

        normal_Vec[0] = -1 * TMath::Sin(Rand_thet) * TMath::Cos(Rand_phi);
        normal_Vec[1] = -1 * TMath::Sin(Rand_thet) * TMath::Sin(Rand_phi);
        normal_Vec[2] = -1 * TMath::Cos(Rand_thet); //## Es un vector normal unitario apuntando  hacia el centro de coordenadas


        // ### ==================== Intersecciones con cada cara =====================  ####
        bool flag_cara_1 = false; 
        bool flag_cara_2 = false; 
        bool flag_cara_3 = false; 
        bool flag_cara_4 = false; 
        bool flag_cara_5 = false;
        bool flag_cara_6 = false;  

        // #### Cara Superior ###
        double t_1 = (medida_z - random_plane_point[2]) / normal_Vec[2];
        double x_1 = random_plane_point[0] + normal_Vec[0] * t_1; 
        double y_1 = random_plane_point[1] + normal_Vec[1] * t_1; 

        // #### Cara Inferior ###
        double t_2 = (-medida_z - random_plane_point[2]) / normal_Vec[2];
        double x_2 = random_plane_point[0] + normal_Vec[0] * t_2;
        double y_2 = random_plane_point[1] + normal_Vec[1] * t_2;

        // ### Caras en X ###
        // ### Cara 3 ###
        double t_3 = (medida_x - random_plane_point[0]) / normal_Vec[0];
        double z_3 = random_plane_point[2] + normal_Vec[2] * t_3;
        double y_3 = random_plane_point[1] + normal_Vec[1] * t_3;

        // ### Cara 4 ###
        double t_4 = (-medida_x - random_plane_point[0]) / normal_Vec[0];
        double z_4 = random_plane_point[2] + normal_Vec[2] * t_4;
        double y_4 = random_plane_point[1] + normal_Vec[1] * t_4;

        // #### Caras en Y ###
        // ### Cara 3 ###
        double t_5 = (medida_y - random_plane_point[1]) / normal_Vec[1];
        double z_5 = random_plane_point[2] + normal_Vec[2] * t_5;
        double x_5 = random_plane_point[0] + normal_Vec[0] * t_5;

        // ### Cara 4 ###
        double t_6 = (-medida_y - random_plane_point[1]) / normal_Vec[1];
        double z_6 = random_plane_point[2] + normal_Vec[2] * t_6;
        double x_6 = random_plane_point[0] + normal_Vec[0] * t_6;

        bool flag_faces = true;
        int n_flags = 0;
        while (flag_faces){
            if ((-medida_x < x_1) & (x_1 < medida_x) & (-medida_y < y_1) & (y_1 < medida_y)){
                flag_cara_1 = true;
                n_flags += 1;
                if (n_flags == 2){
                    flag_faces = false;
                } }
            if ((-medida_x < x_2) & (x_2 < medida_x) & (-medida_y < y_2) & (y_2 < medida_y)){
                flag_cara_2 = true;
                n_flags += 1;
                if (n_flags == 2){
                    flag_faces = false;
                }}
            if ((-medida_y < y_3) & (y_3 < medida_y) & (-medida_z < z_3) & (z_3 < medida_z)){
                flag_cara_3 = true;
                n_flags += 1;
                if (n_flags == 2){
                    flag_faces = false;
                }}
            if ((-medida_y < y_4) & (y_4 < medida_y) & (-medida_z < z_4) & (z_4 < medida_z)){
                flag_cara_4 = true;
                n_flags += 1;
                if (n_flags == 2){
                    flag_faces = false;
                }}
            if ((-medida_x < x_5) & (x_5 < medida_x) & (-medida_z < z_5) & (z_5 < medida_z)){
                flag_cara_5 = true;
                n_flags += 1;
                if (n_flags == 2){
                    flag_faces = false;
                }}
            if ((-medida_x < x_6) & (x_6 < medida_x) & (-medida_z < z_6) & (z_6 < medida_z)){
                flag_cara_6 = true;
                n_flags += 1;
                if (n_flags == 2){
                    flag_faces = false;
                }}
            else{flag_faces = false;}
            };
        list_flags[0] = flag_cara_1;
        list_flags[1] = flag_cara_2;
        list_flags[2] = flag_cara_3;
        list_flags[3] = flag_cara_4;
        list_flags[4] = flag_cara_5;
        list_flags[5] = flag_cara_6;

        list_z[0] = z_3;
        list_z[1] = z_4;
        list_z[2] = z_5;
        list_z[3] = z_6;


        Delta_L = intersection_CCD(list_flags, list_z, medida_z, Rand_thet);

        if ((0 < Delta_L) & (Delta_L < max_long)){
            
            f_LandVav->SetParameters(Delta_L, momentum);
            Rand_edep = f_LandVav->GetRandom();   // MeV
            kappa = Kappa(Delta_L, momentum);
            MuonID = i;

            muons_in_CCD += 1;
            cout << "Muon " << i << "/" << number_thet << "\r" << endl;
            // cout << Delta_L << endl;

            // if (i == number_thet-1){
            //     cout << "Muon " << i+1 << "/" << number_thet << "\r" << endl;
            // }
        }
        else{
            Delta_L = 0;
            Rand_edep = 0;
            kappa = -1;
            MuonID = i;

            cout << "Muon " << i << "/" << number_thet << "\r" << endl;
        }
        tree->Fill();
    }; // end of loop

    // Len_th = sizeof(list_thet_in_CCD);
    // Len_phi = sizeof(list_phi_in_CCD);
    // Len_epri = sizeof(list_energy_pri_in_CCD);

    auto stop = chrono::high_resolution_clock::now();

    using namespace std::chrono;
    auto duration = duration_cast<microseconds>(stop - start);

    cout<< "Muons in CCD: " << muons_in_CCD << endl;
    cout<< "Execution time: " << duration.count() / 1000000 << " seconds" << endl;


    // cout<< "Len thet list: " <<Len_th << endl;
    // cout<< "Len phi list: " <<Len_phi << endl;
    // cout<< "Len epri list: " <<Len_epri << endl;

    tree->Write();
    file->Close();
}