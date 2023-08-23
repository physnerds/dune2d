void test_fluxes(){
  std::string rnom = "test_nominal.root";
  std::string rfid = "test_fidVol.root";

  auto  nom = new TFile(rnom.c_str(),"READ");
  auto fid = new TFile(rfid.c_str(),"READ");
  std::string r_center = "numu/1D/hFluxCenter";
  std::string r_fid = "numu/1D/hFluxRndDetXYZ";
  
  auto rnom_center = (TH1D*)nom->Get(r_center.c_str());
  auto rnom_fid = (TH1D*)nom->Get(r_fid.c_str());
  auto rfid_center = (TH1D*)fid->Get(r_center.c_str());
  auto rfid_fid = (TH1D*)fid->Get(r_fid.c_str());
  
  auto c = new TCanvas();
  c->cd();
  rnom_center->SetLineColor(1);
  rfid_center->SetLineColor(2);
  rnom_center->SetTitle("Default");
  rfid_center->SetTitle("Fid. ND Center");
  rfid_fid->SetTitle("Fid. ND Volume");
  rnom_center->GetXaxis()->SetRangeUser(0,20);
  rnom_center->Draw("hist");
  rfid_center->Draw("histSame");
  c->BuildLegend();
  c->Print("DUNEND_Comparisons_CentralFluxes.png");
  
  auto c2 = new TCanvas();
  c2->cd();
  rnom_center->SetLineColor(1);
  rfid_fid->SetLineColor(2);
  rnom_center->GetXaxis()->SetRangeUser(0,20);
  rnom_center->Draw("hist");
  rfid_fid->Draw("histSame");
  c2->BuildLegend();
  c2->Print("DUNEND_Comparisons_FidFluxes.png");
}
