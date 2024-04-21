#|/bin/bash
root -l << EOF
new TTree("tree", "tree");
tree->ReadFile("$1");
tree->Draw("gamma:mass_w-124.77:chi2_over_ndf", "ptbin==$2", "goff");
int n = tree->GetSelectedRows();
new TCanvas("c");
auto *g = new TGraph(tree->GetSelectedRows(), tree->GetV1(), tree->GetV2());
g->Draw("AP");
for(int i=0; i<n; i++){
   auto *l = new TLatex(tree->GetV1()[i]+0.1, tree->GetV2()[i], Form("%.2f", tree->GetV3()[i]));
   l->SetTextSize(0.01);
   l->Draw();
}
TF1 *f = new TF1("f", "-sqrt(x)*[0]+[1]", 0, 60);
g->Fit(f, "Q");
cout << $2 << "," << f->GetParameter(0) << "," << f->GetParError(0) << "," << f->GetParameter(1) << "," << f->GetParError(1) << endl;
c->SaveAs("$3.root");
c->SaveAs("$3.pdf");
EOF
