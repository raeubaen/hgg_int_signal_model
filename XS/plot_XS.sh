root << EOF
new TTree("tree", "tree");
tree->ReadFile("$1");
tree->Draw("gamma:xsec*1000:abs(xsec*0.3*1000):0*xsec", "", "goff");
new TCanvas("c");
auto *g = new TGraphErrors(tree->GetSelectedRows(), tree->GetV1(), tree->GetV2(), tree->GetV4(), tree->GetV3());
g->Draw("AP");
g->GetXaxis()->SetTitle("#Gamma/#Gamma_{SM}");
g->GetYaxis()->SetTitle("XS [fb]");
c->SaveAs("XS_$2.root");
EOF
