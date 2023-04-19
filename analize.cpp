#include "TCanvas.h" 
#include "TF1.h"
#include "TGraph.h"

/*Double_t Amplitude(Double_t *x, Double_t *par) {
  // par[0] = V0, par[1] = R, par[2] = L, par[3] = C
  Float_t xx = x[0];
  return (par[0] * par[1]) /
         sqrt(par[1] * par[1] +
              (par[2] * par[2] * xx * xx) / ((par[2] * par[3] * xx * xx - 1) *
                                             (par[2] * par[3] * xx * xx - 1)));
}*/

Double_t Amplitude(Double_t *x, Double_t *par) {
  // par[0] = V0, par[1] = R, par[2] = L, par[3] = C
  Float_t xx = x[0];
  return par[0] *
         abs(par[1] *
             (par[2] * par[3] * (xx / (2 * M_PI)) * (xx / (2 * M_PI)) - 1)) /
         sqrt(
             par[1] * par[1] *
                 (par[2] * par[3] * (xx / (2 * M_PI)) * (xx / (2 * M_PI)) - 1) *
                 (par[2] * par[3] * (xx / (2 * M_PI)) * (xx / (2 * M_PI)) - 1) +
             par[2] * par[2] * (xx / (2 * M_PI)) * (xx / (2 * M_PI)));
}

void NotchGraph(Double_t V0, Double_t R, Double_t L, Double_t C) {
  TCanvas *c = new TCanvas("c", "myCanvas", 200, 200, 1000, 600);
  TGraph *g = new TGraph("NotchData.txt", "%lg%lg");
  TF1 *f = new TF1("f", Amplitude, 2400.03, 22700.2, 4);
  f->SetParameters(V0, R, L, C);
  g->Fit(f, "S0");
  g->Draw("AP");
  f->Draw("same");
}

/*void NotchFunc(Double_t V0, Double_t R, Double_t L, Double_t C) {
  TCanvas *c = new TCanvas("c", "myCanvas", 200, 200, 1000, 600);
  TF1 *f = new TF1("f", Amplitude, 0, 1, 4);
  f->SetParameters(V0, R, L, C);
  f->Draw();
}*/

