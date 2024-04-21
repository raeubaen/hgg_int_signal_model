root << EOF
new TTree("tree", "tree");
tree->ReadFile("${1}.csv");
tree->Draw("lowetaedge:relerror", "", "goff");
new TH1F("h_${2}", "h", 19, tree->GetV1());
for(int i=1; i<tree->GetSelectedRows(); i++){ h_${2}->SetBinContent(i, tree->GetV2()[i-1]); h_${2}->SetBinError(i, 1e-6); }
h_${2}->SaveAs("${1}.root");
EOF
