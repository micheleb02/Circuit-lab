#include <fstream>

#include "TCanvas.h"
#include "TF1.h"
#include "TH1F.h"
#include "TStyle.h"

void NotchGaussV2() {  // fit sulla distribuzione di ampiezze
  gStyle->SetOptFit(1111);
  TH1F *h = new TH1F("h", "Fitting A(#omega)", 20, 2.4185,
                     2.42);                      // controllare gli estremi
  TF1 *f = new TF1("f", "gaus", 2.4185, 2.42);  // controllare gli estremi
  TCanvas *c = new TCanvas("c", "MyCanvas", 200, 200, 1000, 600);
  std::ifstream in;
  in.open("A_R.dat");  // nome del file
  Double_t x;
  Double_t y;
  while (1) {
    in >> x >> y;
    if (!in.good()) {
      break;
    }
    if (x == 5000.11) {
      h->Fill(y);
    }
  }
  in.close();
  h->Fit("f");
  h->Draw("hist");
  f->Draw("same");
}