void Parabola_Fit_Linear (float sigma=0.1, float a=0.05) {
// Parabola_Fit_Linear(0.05,0.05)
// Subir hasta a=0.10
// 
const int npoints=11;
double xx[npoints] = {-5,-4,-3,-2,-1,0,+1,+2,+3,+4,+5};
double ex[npoints] = {0};    // sin errores en X
double ey[npoints] = {0};    // errores Y 
double yyL[npoints] = {0};
double yyP[npoints] = {0};

double b=0.3;

TRandom3 r(0);
for (int i=0; i<npoints; i++) {
   double x = xx[i];
   ey[i] = sigma;
   yyP[i] = r.Gaus (b*x+a*x*x,ey[i]);  // data de H1: parabola
}

auto grL = new TGraphErrors(npoints,xx,yyL,ex,ey);
auto grP = new TGraphErrors(npoints,xx,yyP,ex,ey);
auto grQ = new TGraphErrors(npoints,xx,yyP,ex,ey);

gStyle->SetOptStat(0);
gStyle->SetOptFit(11);

auto c1 = new TCanvas("c1","Lineal vc Parabola",700,1000);
c1->Divide(1,2);
c1->cd(1);

grP->SetMarkerColor(1);
grP->SetMarkerStyle(20);
grP->SetTitle(";X;Y");
grP->Draw("apm");
gPad->SetGridx();
gPad->SetGridy();

double xmin = grL->GetXaxis()->GetXmin();
double xmax = grL->GetXaxis()->GetXmax();
auto fLinear = new TF1("fLinear","        [b]*x+[c]",xmin,xmax);
auto fParabo = new TF1("fParabo","[A]*x*x+[B]*x+[C]",xmin,xmax);
fLinear->SetParameters(b,1);
fParabo->SetParameters(0.03,b,1);
fLinear->SetLineColor(kGreen+2);
fParabo->SetLineColor(kBlue);
grP->Fit(fLinear);

c1->cd(2);
grQ->Fit(fParabo);
grQ->SetMarkerColor(1);
grQ->SetMarkerStyle(20);
grQ->SetTitle(";X;Y");
grQ->Draw("apm");
gPad->SetGridx();
gPad->SetGridy();

return;
}