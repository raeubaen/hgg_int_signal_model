#!/usr/bin/bash
#1 file smearato - 2 tipo - 4 normalizzazione - 3 gamma
#va chiamato da un python
root $1 << EOF
auto *h = new TH1F("mass", "mass", 35, 105, 140);
tree->Draw("mass>>mass", "(weight/1500 * $4) * (lead_pt > 40 && sublead_pt > 30)");
h->SaveAs("histos/mass_${2}_${3}_th.root");

EOF
