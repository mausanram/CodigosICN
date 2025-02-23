void Test_edep(){

TFile *file = new TFile("./root_files/muons_2M_vacuum_file.root");
// TFile *file = new TFile("./root_files/muons_1M_vacuum_file.root");
TTree *tree = (TTree*) file->Get("B02Evts");
//TTree *tree = (TTree*) file->Get("B02Hits");

TFile *file0 = new TFile("Edep_allclusters_NSAMP324_MeV.root"); // INFO ALL_CLUSTERS
// TFile *file0 = new TFile("Edep_NSAMP324_KeV.root");	// INFO MUONS ONLY
TTree *tree0 = (TTree*) file0->Get("tree");

TFile *file1 = new TFile("../../../../Simulacion_ab_initio/Sim_ab_initio_NMUONS_100000_PLANES_1.5_RADIO_8_CCDSIZE_400X600_C.root");
TTree *tree1 = (TTree*) file1->Get("tree");

TFile *file2 = new TFile("Edep_CONNIE_NSAMP400_MeV.root");
TTree *tree2 = (TTree*) file2->Get("tree");

TFile *file3 = new TFile("../../../../Simulacion_ab_initio/Sim_ab_initio_NMUONS_100000_PLANES_1.5_RADIO_8_CCDSIZE_400X600_C.root");
TTree *tree3 = (TTree*) file3->Get("tree");

TLatex lat;

int NB = 200;
double tlow = 0;
double thi = 1000;

TH1F *edep = new TH1F("edep", "Energy Spectrum", NB, tlow, thi);
edep->GetXaxis()->SetTitle("Energy (MeV)");
edep->SetLineStyle(1);
edep->SetLineColor(3);
edep->SetStats(0);

//TH1F *edep0 = new TH1F("edep0", "Lab. Det. ", NB, tlow, thi);
TH1F *edep0 = new TH1F("edep0", "Energy Spectrum (Sim. PP - All Clusters (ICN))", NB, tlow, thi);
edep0->GetXaxis()->SetTitle("Energy (MeV)");
edep0->SetLineStyle(1);
edep0->SetLineColor(1);
edep0->SetStats(0);

TH1F *edep1 = new TH1F("edep1", "Simulation no Birks", NB, tlow, thi);
edep1->SetStats(0);
edep1->SetLineStyle(1);
edep1->SetLineColor(2);

TH1F *edep2 = new TH1F("edep2", "CONNIE 2021-2022", NB, tlow, thi);
edep2->SetStats(0);
edep2->SetLineStyle(2);
edep2->SetLineColor(1);

TH1F *edep3 = new TH1F("edep3", "Simulation PP", NB, tlow, thi);
edep3->SetStats(0);
edep3->SetLineStyle(1);
edep3->SetLineColor(4);

TH1F *edep4 = new TH1F("edep4", "Simulation no Birks scale", NB, tlow, thi);
edep4->SetStats(0);
edep4->SetLineStyle(2);
edep4->SetLineColor(1);

// === Template for fit ==
int NBmu = 200;
double maxEd = 700.;

double muonsbw = maxEd/NBmu;

TH1F *muons = new TH1F("muons","",NBmu,0,maxEd);
muons->GetXaxis()->SetTitle("Energy (KeV)");
muons->GetYaxis()->SetTitle("events");

// TH1F *muons = new TH1F("muons", "Simulation no Birks", NBmu, 0, maxEd);
tree->Draw(" EevtBar>>muons", "EevtBar>0"); // GEANT4 INFO (NO BIRKS)

// ============= Fill histograms =========== //
tree->Draw("WevtBar>>edep", "WevtBar>0"); // GEANT4 INFO (BIRKS)
tree->Draw(" EevtBar*1000>>edep1", "EevtBar>0"); // GEANT4 INFO (NO BIRKS)
//tree->Draw("Ehitbar>>edep_cut", "thet>22*TMath::Pi()/180 & edep>0");
// ========================================= //

tree0->Draw("edep*1000>>edep0", "edep>0"); // EXPERIMENTAL INFO
// tree0->Draw("edep>>edep0", "edep>0"); // EXPERIMENTAL INFO
tree2->Draw("edep>>edep2", "edep>0"); // EXPERIMENTAL INFO CONNIE
tree3->Draw("edep*1000>>edep3", "edep>0"); // SIM_AB_INITIO INFO

tree->Draw("EevtBar*0.8*1000>>edep4", "EevtBar>0"); // GEANT4 INFO (NO BIRKS)


tree->Draw("EevtBar*1000*0.8>>muons", "EevtBar>0");
// tree0->Draw("edep>>muons", "edep>0");

double cont = edep->Integral();
double cont0 = edep0->Integral();
double cont1 = edep1->Integral();
double cont3 = edep3->Integral();
double cont4 = edep4->Integral();

cout << "Integral }gean4: " <<  cont1 <<endl;

// ========== Scale histograms =========== //
//edep->Scale(0.65, "");
edep->Scale(1);
//edep->SetLineColor(2);

edep0->Scale(1.5);
// edep0->Scale(50);
//edep0->SetLineColor(4);

edep3->Scale(cont1/cont3);
//edep1->SetLineColor(1);

edep2->Scale(0.86);
//edep2->SetLineColor(7);

// edep1->Scale(cont3/cont1);

// ======================================= //



// ============ Create Canvas ============== //
TCanvas *canv = new TCanvas("canv","Edep", 2*800, 600);
canv->Divide(1,1);
canv->cd(1);
// canv->Grid();
// edep->Draw("hist"); 		// edepG4 with birks
// edep4->Draw("hist");
// edep0->Draw("hist same"); 	// ICN data 
edep1->Draw("hist "); 	// edepG4 no Birks
// edep2->Draw("he0 same"); 	// CONNIE data
edep3->Draw("hist same"); 	// edepPP

TLegend *leg = new TLegend(0.5, 0.7, 0.9, 0.9);
// leg->AddEntry(edep, "SimG4-Birks: 0.09 cm/MeV", "lep");
leg->AddEntry(edep1, "Sim-GEANT4", "LEP");
leg->AddEntry(edep0, "All Clusters (ICN-NSAMP324)", "LEP");
//leg->AddEntry(edep2, "Datos CONNIE-NSAMP400", "LEP");
leg->AddEntry(edep3, "Sim-PP", "LEP");
leg->AddEntry(edep4, "Sim-GEANT4 (0.77)", "LEP");
leg->Draw();



TFile *fout = new TFile("ccdhisto.root", "recreate");
edep0->Write();
fout->Close();


//----------------------------
//-- muons template data file 
//----------------------------
  
  ofstream myfile;
  myfile.open ("muons.dat");
  for (int i=0;i<NBmu; i++) 
  myfile << muons->GetBinContent(i+1)/(muons->Integral(1,NBmu)*muonsbw) << "\n";
  myfile.close();

  cout << "Integral muons: " << muons->Integral() << " entries." << endl;
  cout << "Integral muons (+ovflw): " << muons->Integral(0,NBmu+1) << " entries." << endl;

  //- Find muon peak
  int peakBin = 50;
  for (int i=50;i<NBmu;i++){
   if (muons->GetBinContent(i+1) > muons->GetBinContent(peakBin))
     peakBin = i+1;
  } //for i
  printf("Peak: %d  \n", peakBin);
  printf("Peak: %3.3f MeV \n", peakBin*maxEd/NBmu);
  double Epeak = peakBin*maxEd/NBmu;




//canv->cd(2);
//edep_cut->Draw();
// func2->Draw("same");
canv->Print("Dis_edep_GEANT-PP.pdf");

}
