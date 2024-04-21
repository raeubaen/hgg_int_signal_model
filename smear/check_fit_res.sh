root $1 << EOF
new TTree("temp", "temp");
temp->ReadFile("reso_eta_inclusive_2018.csv");
temp->Draw("lowetaedge:relerror", "", "goff");
new TTree("temp2", "temp2");
temp2->ReadFile("fitted_reso_inc.csv");
temp2->Draw("sigma:errsigma:0*errsigma", "", "goff");
auto *g = new TGraph(temp->GetSelectedRows(), temp->GetV1(), temp->GetV2());
auto *ge = new TGraphErrors(temp->GetSelectedRows(), temp->GetV1(), temp2->GetV1(), temp->GetV3(), temp2->GetV2());
new TCanvas("c");
g->Draw("AP");
ge->Draw("P");
c->SaveAs("check.root");
EOF
