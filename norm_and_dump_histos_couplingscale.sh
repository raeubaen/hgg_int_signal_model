#!/usr/bin/bash
#1 file smearato - 2 tipo - 4 normalizzazione - 3 gamma
#va chiamato da un python
root $1 << EOF
new TTree("binstree", "binstree");
binstree->ReadFile("higgsptbins.csv");
binstree->Draw("lowhptedge:lowhptedge", "", "goff");
double *lowhptedges = binstree->GetV1();

cout << "nbins: " << binstree->GetSelectedRows()-1 << endl;

for(int i=0; i<binstree->GetSelectedRows()-1; i++){
  auto *h = new TH1F(Form("mass_hptbin%i", i), "mass", 20, 120, 130);
  if (lowhptedges[i]==-1) tree->Draw(Form("mass>>mass_hptbin%i", i), "(weight/1500 * $4) * (lead_pt > 40 && sublead_pt > 30) * sqrt($3)");
  else tree->Draw(Form("mass>>mass_hptbin%i", i), Form("weight/1500 * $4 * (higgs_pt > %f && higgs_pt <= %f) * (lead_pt > 40 && sublead_pt > 30) * sqrt($3)", lowhptedges[i], lowhptedges[i+1]));
  h->SaveAs(Form("histos/mass_couplingscale_${2}_${3}_hptbin%i.root", i));
  delete h;
}
EOF
