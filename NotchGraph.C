#include "TCanvas.h"
#include "TF1.h"
#include "TGraph.h"
#include "TStyle.h"

Double_t AmplitudeR(Double_t *x, Double_t *par) {
  // par[0] = V0, par[1] = R, par[2] = L, par[3] = C
  Float_t xx = x[0];
  return par[0] *
         abs(par[1] *
             (par[2] * par[3] * (xx * (2 * M_PI)) * (xx * (2 * M_PI)) - 1)) /
         sqrt(
             par[1] * par[1] *
                 (par[2] * par[3] * (xx * (2 * M_PI)) * (xx * (2 * M_PI)) - 1) *
                 (par[2] * par[3] * (xx * (2 * M_PI)) * (xx * (2 * M_PI)) - 1) +
             par[2] * par[2] * (xx * (2 * M_PI)) * (xx * (2 * M_PI))) /* +
     0.545479*/
      ;
}

Double_t AmplitudeLC(Double_t *x, Double_t *par) {
  // par[0] = V0, par[1] = R, par[2] = L, par[3] = C
  Float_t xx = x[0];
  return par[0] * par[2] * (xx * (2 * M_PI)) /
         sqrt(
             par[1] * par[1] *
                 (par[2] * par[3] * (xx * (2 * M_PI)) * (xx * (2 * M_PI)) - 1) *
                 (par[2] * par[3] * (xx * (2 * M_PI)) * (xx * (2 * M_PI)) - 1) +
             par[2] * par[2] * (xx * (2 * M_PI)) * (xx * (2 * M_PI)));
}

Double_t GetF0Err(Double_t *x, Double_t *par) {
  // par[0] = L, par[1] = C, par[2] = DL, par[3] = DC
  return (sqrt((par[2] * par[2]) / (par[0] * par[0]) +
               (par[3] * par[3]) / (par[1] * par[1]))) /
             (2 * M_PI * sqrt(par[0] * par[1])) +
         x[0];
}

Double_t GetQErr(Double_t *x, Double_t *par) {
  // par[0] = R, par[1] = L, par[2] = C, par[3] = DR, par[4] = DL, par[5] = DC
  return sqrt(par[2] * par[3] * par[3] / par[1] +
              par[0] * par[0] * par[5] * par[5] / (4 * par[1] * par[2]) +
              par[0] * par[0] * par[2] * par[4] * par[4] /
                  (4 * par[1] * par[1] * par[1])) +
         x[0];
}

void NotchGraphR(Double_t V0, Double_t R, Double_t L, Double_t C) {
  gStyle->SetOptFit(1111);
  TCanvas *c = new TCanvas("c", "myCanvas", 200, 200, 1000, 600);
  TGraph *g = new TGraph("NotchData.txt", "%lg%lg");
  TF1 *f = new TF1("f", AmplitudeR, 2400.04, 23000.1, 4);
  g->SetMarkerStyle(25);
  f->SetParameters(V0, R, L, C);
  g->Fit(f, "S0");
  g->Draw("AP");
  f->Draw("same");
  std::cout << "w0 :" << f->GetMinimumX() << '\n';
}

void NotchGraphLC(Double_t V0, Double_t R, Double_t L, Double_t C) {
  gStyle->SetOptFit(1111);
  TCanvas *c = new TCanvas("c", "myCanvas", 200, 200, 1000, 600);
  TGraph *g = new TGraph("MagliaData.txt", "%lg%lg");
  TF1 *f = new TF1("f", AmplitudeLC, 2400.03, 22700.2, 4);
  g->SetMarkerStyle(25);
  f->SetParameters(V0, R, L, C);
  g->Fit(f, "S0");
  g->Draw("AP");
  f->Draw("same");
  std::cout << "w0 :" << f->GetMaximumX() << '\n';
}

void GetF0AndQ_R(Double_t V0, Double_t R, Double_t L, Double_t C) {
  TGraph *g = new TGraph("NotchData.txt", "%lg%lg");
  TF1 *fit = new TF1("fit", AmplitudeR, 2400.04, 23000.1, 4);
  TF1 *F0err = new TF1("err", GetF0Err, 0, 1, 4);
  TF1 *Qerr = new TF1("Qerr", GetQErr, 0, 1, 6);
  fit->SetParameters(V0, R, L, C);
  g->Fit(fit, "S0");
  F0err->SetParameters(fit->GetParameter(2), fit->GetParameter(3),
                       fit->GetParError(2), fit->GetParError(3));
  Qerr->SetParameters(fit->GetParameter(1), fit->GetParameter(2),
                      fit->GetParameter(3), fit->GetParError(1),
                      fit->GetParError(2), fit->GetParError(3));
  std::cout << "\n";
  std::cout << "*****************************" << '\n';
  std::cout << "f0 = "
            << 1 / (2 * M_PI *
                    sqrt(fit->GetParameter(2) * fit->GetParameter(3)))
            << " +/- " << F0err->Eval(0) << '\n';
  std::cout << "Q = "
            << fit->GetParameter(1) *
                   sqrt(fit->GetParameter(3) / fit->GetParameter(2))
            << " +/- " << Qerr->Eval(0) << '\n';
}

void GetF0AndQ_LC(Double_t V0, Double_t R, Double_t L, Double_t C) {
  TGraph *g = new TGraph("MagliaData.txt", "%lg%lg");
  TF1 *fit = new TF1("fit", AmplitudeLC, 2400.03, 22700.2, 4);
  TF1 *F0err = new TF1("err", GetF0Err, 0, 1, 4);
  TF1 *Qerr = new TF1("Qerr", GetQErr, 0, 1, 6);
  fit->SetParameters(V0, R, L, C);
  g->Fit(fit, "S0");
  F0err->SetParameters(fit->GetParameter(2), fit->GetParameter(3),
                       fit->GetParError(2), fit->GetParError(3));
  Qerr->SetParameters(fit->GetParameter(1), fit->GetParameter(2),
                      fit->GetParameter(3), fit->GetParError(1),
                      fit->GetParError(2), fit->GetParError(3));
  std::cout << "\n";
  std::cout << "*****************************" << '\n';
  std::cout << "f0 = "
            << 1 / (2 * M_PI *
                    sqrt(fit->GetParameter(2) * fit->GetParameter(3)))
            << " +/- " << F0err->Eval(0) << '\n';
  std::cout << "Q = "
            << fit->GetParameter(1) *
                   sqrt(fit->GetParameter(3) / fit->GetParameter(2))
            << " +/- " << Qerr->Eval(0) << '\n';
}