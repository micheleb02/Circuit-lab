#include <fstream>

#include "TCanvas.h"
#include "TF1.h"
#include "TH1F.h"

Double_t Gauss(Double_t *x, Double_t *par) {
  return par[0] * exp(-(x[0] - par[1]) * (x[0] - par[1]) / (2 * par[2]));
}

void NotchGauss() {
  TH1F *h = new TH1F("h", "Gauss", 50, 2.052, 2.0528);
  TF1 *f = new TF1("f", Gauss, 2.052, 2.0528, 3);
  f->SetParameter(1, 2.051);
  TCanvas *c = new TCanvas("c", "MyCanvas", 200, 200, 1000, 600);
  std::ifstream in;
  in.open("GaussA.dat");
  Float_t x, y;
  while (1) {
    in >> x >> y;
    if (!in.good()) {
      break;
    }
    h->Fill(y);
  }
  in.close();
  h->Fit(f, "S0");
  h->Draw("histo");
  f->Draw("same");
  //std::cout << "RMS: " << h->GetRMS() << " +/- " << h->GetRMSError() << '\n';
}