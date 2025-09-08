void TestDiffMod(){
TChain *tree = new TChain("tree");
tree->Add("tree_DiffusionMod_Ext1_Vert.root");

float sprd_val;
float deep_val;

tree->SetBranchAddress("sprd", &sprd_val);
tree->SetBranchAddress("deep", &deep_val);

int n_entries = tree->GetEntries();
double arr_sprd[n_entries];
double arr_deep[n_entries];

// std::cout<<tree->GetEntries()<<endl;

// Fill histograms //
// tree->Draw("sprd:deep");

for (int i = 0; i < tree->GetEntries(); i++){
    tree->GetEntry(i);
    arr_sprd[i] = sprd_val;
    arr_deep[i] = deep_val;
};

TCanvas *c1 = new TCanvas("c1");
TGraph *diff_data = new TGraph(n_entries, arr_deep, arr_sprd);
diff_data->SetMarkerStyle(2);
diff_data->GetXaxis()->SetRangeUser(0, 730);
diff_data->GetYaxis()->SetRangeUser(0, 1.3);
diff_data->GetXaxis()->SetTitle("Deep (\mu m)");
diff_data->GetYaxis()->SetTitle("Spread (px)");
diff_data->SetTitle("Diffusion Model (Ext1-Vert)");

TF1 *diff_curve = new TF1("diff_curve", "sqrt([0]* log(1 - [1]*x))/15", 0, 725);
diff_curve->SetParameters(-300, 0.0004);
diff_data->Fit(diff_curve);

TGraph *frame = new TGraph(10, 0, 725);

// Create Canvas //
// TCanvas *canv = new TCanvas("canv","", 2*700, 600);
diff_data->Draw("AP");
diff_curve->Draw("same");

// TLegend *leg = new TLegend(0.7, 0.8, 0.9, 0.9);
// // leg->SetTextAlign(11);
// leg->SetFillStyle(0);
// leg->AddEntry(diff_curve, "Diff. Mod, fit: ", "l");
// leg->AddEntry(line, "Grosor de la CCD: 0.0725 cm", "l");
// // leg->AddEntry(line, "Grosor de la CCD: 0.068 cm", "l");
// leg->Draw();

std::cout<<"Chi2: " << diff_curve->GetChisquare() << endl;
std::cout<<"NDF: " << diff_curve->GetNDF() << endl;
std::cout<<"Prob: " << diff_curve->GetProb() << endl;

}
