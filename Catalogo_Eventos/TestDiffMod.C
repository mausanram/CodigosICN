#include <string>

void TestDiffMod(){
TChain *tree = new TChain("tree");
// tree->Add("tree_DiffusionMod_Ext1_Vert.root");
// tree->Add("tree_DiffusionMod_Ext1_Horz.root");
tree->Add("tree_DiffusionMod_Ext2_Vert.root");
// tree->Add("tree_DiffusionMod_Ext2_Horz.root");

float sprd_val;
float deep_val;

tree->SetBranchAddress("sprd", &sprd_val);
tree->SetBranchAddress("deep", &deep_val);

int nentries = tree->GetEntries();
double arr_sprd[1000];
double arr_deep[1000];

// std::cout<<tree->GetEntries()<<endl;

// Fill histograms //
// tree->Draw("sprd:deep");

for (int i = 0; i < tree->GetEntries(); i++){
    tree->GetEntry(i);
    arr_sprd[i] = sprd_val;
    arr_deep[i] = deep_val;
};

TCanvas *c1 = new TCanvas("c1");
TGraph *diff_data = new TGraph(nentries, arr_deep, arr_sprd);
diff_data->SetMarkerStyle(2);
diff_data->GetXaxis()->SetRangeUser(0, 730);
diff_data->GetYaxis()->SetRangeUser(0, 1.3);
diff_data->GetXaxis()->SetTitle("Profundidad (#mu m)");
diff_data->GetYaxis()->SetTitle("Anchura (px)");
// diff_data->SetTitle("Modelo de Difusi#acute{o}n(Extensi#acute{o}n 1)");
diff_data->SetTitle("Modelo de Difusi#acute{o}n(Extensi#acute{o}n 2)");

TF1 *diff_curve = new TF1("diff_curve", "sqrt([0]* log(1 - [1]*x))/15", 0, 725);
diff_curve->SetParameters(-300, 0.0004);
diff_data->Fit(diff_curve);

TGraph *frame = new TGraph(10, 0, 725);

// Create Canvas //
// TCanvas *canv = new TCanvas("canv","", 2*700, 600);
diff_data->Draw("AP");
diff_curve->Draw("same");

TLegend *leg = new TLegend(0.7, 0.8, 0.9, 0.9);
// leg->SetTextAlign(11);
// leg->SetFillStyle(0);
// leg->AddEntry(diff_curve, "#sqrt{#alpha log(1 - #beta z)} fit: #alpha = -306.329 #mu m^{2}, #beta = 0.0007903 #mu m^{-1}", "l"); // Estension 1
leg->AddEntry(diff_curve, "#sqrt{#alpha log(1 - #beta z)} fit: #alpha = -332.143 #mu m^{2}, #beta = 0.0007341 #mu m^{-1}", "l"); // Estension 2
// leg->AddEntry(line, "Grosor de la CCD: 0.068 cm", "l");
leg->Draw();

std::cout<<"Chi2: " << diff_curve->GetChisquare() << endl;
std::cout<<"NDF: " << diff_curve->GetNDF() << endl;
std::cout<<"Prob: " << diff_curve->GetProb() << endl;


TCanvas *canv = new TCanvas("canv","", 2*700, 600);
TF1 *diff_curve_0 = new TF1("diff_curve_0", "sqrt([0]* log(1 - [1]*x))/15", 0, 725);
diff_curve_0->SetLineColor(1);
diff_curve_0->GetXaxis()->SetTitle("Deep (mu m)");
diff_curve_0->GetYaxis()->SetTitle("Spread (px)");
diff_curve_0->SetTitle("Diffusion Model");

TF1 *diff_curve_1 = new TF1("diff_curve_1", "sqrt([0]* log(1 - [1]*x))/15", 0, 725);
diff_curve_1->SetLineColor(2);

TF1 *diff_curve_2 = new TF1("diff_curve_2", "sqrt([0]* log(1 - [1]*x))/15", 0, 725);
diff_curve_2->SetLineColor(4);

TF1 *diff_curve_3 = new TF1("diff_curve_3", "sqrt([0]* log(1 - [1]*x))/15", 0, 725);
diff_curve_3->SetLineColor(6);


diff_curve_0->SetParameters(-306.329, 0.0007903);
diff_curve_1->SetParameters(-844.967, 0.0003287);
diff_curve_2->SetParameters(-332.143, 0.0007341);
diff_curve_3->SetParameters(-487.071, 0.0004444);
diff_curve_0->Draw();
diff_curve_1->Draw("same");
diff_curve_2->Draw("same");
diff_curve_3->Draw("same");

// TLegend *leg = new TLegend(0.7, 0.8, 0.9, 0.9);
// // leg->SetTextAlign(11);
// // leg->SetFillStyle(0);
// leg->AddEntry(diff_curve_0, "", "l");
// leg->AddEntry(diff_curve_1, "Ext1-Horz", "l");
// leg->AddEntry(diff_curve_2, "Ext2-Vert", "l");
// leg->AddEntry(diff_curve_3, "Ext2-Horz", "l");
// leg->Draw();



}
