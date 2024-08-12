// Macro to draw the Landau-Vavilov(-Gauss) formula.
// Here is rewritten the entire Bethe-Bloch formula because is needed in the Landau-Vavilov distribution.
// All the steps were taken from the Leo's book.
// Hare are needed 2 parameters: the distance traveled by the particle inside the bar and its momentum.
#include <iostream> 
#include <random>
#include <ctime>
#include <chrono>

int seed_rand(){
	srand(1003);
	// cout <<  << endl;
	return 0;
} 

// double srand(time(NULL));

double LV (double *lx, double *lpar) {
	double Delta = lx[0];	// Energy loss in absorber
	double L = lpar[0];		// Thickness of absorber (Distance crossed by the particle)

	double p = lpar[1];		// Momentum (in MeV/c)
	double K = 0.307075;	// K coefficient = 4*pi*N*r^2*m*c^2 (in MeV mol^-1 cm^2)
	int z = -1;		// Charge number of incident particle
	double ZA = 0.498487;	// Atomic number over Atomic mass of absorber (for Si)
	double c = TMath::C();	// Speed of light
	double me = 0.510998928;	// Electron mass in MeV/c^2
	double M = 105.65839;	// Muon mass in MeV/c^2
	double I = 0.000173;		// Mean excitation energy (for Si)
	double bg = p/M;
	double beta = bg/sqrt(1+(pow(bg,2)));	// Beta factor
	double gamma = 1/sqrt(1-(pow(beta,2)));	// Gamma factor
	double pi = TMath::Pi();
	double rho = 2.33;	// Density of material (for Si)

	double d;	// Variable for the Density effect

	double a = 0.1492;	// Parameters (taken from W.R. Leo for SI)
	double k = 3.25;		//
	double X0 = 0.2014;		//
	double X1 = 2.87;		//
	double C = -4.44;		//
	// double d0 = 0.0;		//
	double X = log10(bg);
		if (X>=X1) {
			d = 2*log(10.0)*X-C;
			}
		else if (X0<=X && X<X1) {
			d = 2*log(10.0)*X-C+a*(pow((X1-X),k));
			}
		else if (X<X0) {
			// d = d0*(pow(10,(2*(X-X0))));
			d = 0;
			}

	double WM = 2*me*(pow((beta*gamma),2))/(1+(2*me*gamma/M)+pow((me/M),2));   // Maximum energy tranfer

	double loge = log((1-(pow(beta,2)))*(pow(I,2))/(2*M*(pow(beta,2))))+(pow(beta,2)); // log epsilon variable

	double EC = 0.577;	// Euler's constant

	double DeltaAv = K*rho*L*(pow(z,2))*ZA*(1.0/(pow(beta,2)))*((1.0/2.0)*log(2*me*(pow(beta,2))*(pow(gamma,2))*WM/(pow(I,2)))-(pow(beta,2))-(d/2.0));// Mean energy loss (Bethe-Bloch)

	double xi = (K/2)*rho*ZA*L*pow((z/beta),2);		// Xi variable 

	double lambda = (Delta-xi*(log(xi)-loge+1-EC))/xi;	// Lambda parameter

  	double Deltamp = xi*(log(xi/exp(loge))+0.198-d);		// Most probable energy loss
  	double lambdamp = (Deltamp-xi*(log(xi)-loge+1-EC))/xi;

	double kappa = xi/WM;		// Kappa ratio
	double beta2 = pow(beta,2);

	double sigma2 = (pow(xi,2))*(1-beta2/2)/kappa;		// Standard deviation for relativistic particles

	//  cout << lambda << " " << Delta <<"  "<< phi/csi << endl;

	//printf(Deltamp);
	// cout << "Most Probably Energy in KeV " << Deltamp * 1000 << endl;
	//td::cout << Deltamp * 1000 << std::endl;

	if (kappa<=0.01) {
		double phi = TMath::Landau(lambda, lambdamp, 1.0);
		return phi/xi;
		}
	else if (0.01<kappa && kappa<10) {
		double vav = TMath::Vavilov(Delta-Deltamp, kappa, beta2);
		return vav;
		}
	else {
//		double gauss = exp(((Delta-DeltaAv)**2)/(2*sigma2));
		double gauss = TMath::Gaus(Delta, DeltaAv, sqrt(sigma2));
		return gauss;
		}

}

void LandauVavilov_Mau() {
	// gRandom->SetSeed(0);	// Cambia la semilla aleatoria para el GetRandom 

	double s = 0.0725;	// Distance of CCD (in cm)
	// double p = 600; // Momentum parameter (in MeV)
	
	// double En_Smith;
	// char En_Smith_char[100] =  getenv("EN_SMITH");
	float p = atof(getenv("EN_SMITH"));	// Momentum parameter (in MeV)

	// cout << p endl;

	// cout << "Introduce un entero: ";	// ---------------------------------------- //
	// double p;						// Esta sección es para pedir que se ingrese el momento de los muones a mano 
   	// std::cin >> p;					// ---------------------------------------- //

	TCanvas *cnv = new TCanvas("cnv", "", 900, 700);
	cnv->SetGrid();

	TLatex lat;

	TF1 *f = new TF1("f", LV, 0, 100, 2);
	f->SetNpx(1000);


	f->SetParameter(0, s);
	f->SetParameter(1, p);
	f->SetRange(0, 0.7);
	//f->SetTitle("Landau-Vavilov distribution (for 0.0725cm of Si);#font[12]{Energy} (MeV);Probability");
	f->SetTitle("Landau-Vavilov distribution (0.0725cm of Si);Energy (MeV);Probability");
	f->Draw();

	double  If = f->Integral(0,2.0);

	double loE= 0.0;
	double hiE= 0.7;
	int NB = 100;
	double bWidth = (hiE-loE)/NB;

	TH1F *h = new TH1F("h", "", NB, loE, hiE);
	double Edep = 0;

	int SetSeed(0);
	Edep = f->GetRandom(0, 0.7); 
	std::cout << "Edep = "<< Edep * 1000 << " KeV" <<std::endl;
	char buff[100];

	// sprintf(buff,"%f",Edep);

	// // printf("%s", buff);
	// char env[] = "LD_LIBRARY_PATH";
	// putenv(env);

	// setenv("EDEP", buff, 1);
	// setenv("EDEP", Edep, 1);

	
	// cout << "The value of EDEP is " << getenv("PATH")<<" MeV"<<endl;
	// cout << getenv("MYENV") << endl;

	// diagnostico	
	
	// int Nevt = 1000;
	// for (int i = 0 ; i < Nevt; i++){
	// 	Edep = f->GetRandom();
	// 	h->Fill(Edep);
	// }
	// double sf = (If/(bWidth*h->Integral()));
	// h->Scale(sf);
	// h->SetTitle("Landau-Vavilov distribution (0.0725cm of Si);Energy (MeV);Probability");
	// h->Draw();
	// f->Draw("same");
	

	// std::cout << "Most Probably Energy" << std::endl;
	// std::cout << Deltamp * 1000 << std::endl;

}
