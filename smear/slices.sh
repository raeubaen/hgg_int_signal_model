root $1 << EOF
new TTree("temp", "temp");
temp->ReadFile("reso_eta_inclusive_2018.csv");
temp->Draw("lowetaedge:relerror", "", "goff");
for(int i=0; i<19; i++) cout << temp->GetV1()[i] << endl;
new TH2F("h", "h", 19, temp->GetV1(), 50, -0.1, 0.1);
tree->Draw("(lead_en-lead_en_gen)/lead_en_gen:abs(lead_eta)>>h", "lead_has_conv==1", "zcol");
h->FitSlicesY(nullptr, 0, -1, 0, "", nullptr);
h->SaveAs("prova.root");
EOF
