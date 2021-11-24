// ToH: hipotesis simples
// H0: los datos vienen de una dada recta.
// H1: los datos vienen de una dada parabola. 

// STD es un test de chi cuadrado de los datos contra cada hipotesis

void TH_STD(double a=0.03, int Nexp=20000) {

// "a" es el coeficiente cuadrático de la parábola
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
   ey[i] = 0.5;
   yyL[i] = r.Gaus (      b*x,ey[i]);  // data de H0: lineal
   yyP[i] = r.Gaus (a*x*x+b*x,ey[i]);  // data de H1: parabola
}

auto grL = new TGraphErrors(npoints,xx,yyL,ex,ey);
//auto grP = new TGraphErrors(npoints,xx,yyP,ex,ey);

gStyle->SetOptStat(0);
auto c1 = new TCanvas("c1","Lineal vs Parabola",700,1000); 
c1->Divide(1,2);
c1->cd(1);

grL->SetMarkerColor(1);
grL->SetMarkerStyle(20);
grL->SetTitle(";X;Y");
grL->Draw("apm");

gPad->SetGrid();

double xmin = grL->GetXaxis()->GetXmin();
double xmax = grL->GetXaxis()->GetXmax();
auto fLinear = new TF1("fLinear","[0]*x"        ,xmin,xmax);
auto fParabo = new TF1("fParabo","[0]*x+[1]*x*x",xmin,xmax);
fLinear->SetParameter(0,b); // fijo en b la ordenada al origen de la recta
fParabo->SetParameter(0,b); // fijo en b la ordenada al origen de la parabola
fParabo->SetParameter(1,a); // fijo en b el coeficiente cuadrático de la parábola
fLinear->SetLineColor(kGreen+2);
fParabo->SetLineColor(kBlue);
fLinear->DrawCopy("same");
fParabo->DrawCopy("same");

c1->Print("Simple_Hipotesis.pdf["); // open file
c1->Print("Simple_Hipotesis.pdf"); // print canvas
gPad->Update();
//gPad->WaitPrimitive();

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Lo anterior sirvió para mostrar una realización particular de los 11 los puntos y dos modelos
// Ahora se inicia el bucle de Nexp experimentos
// en cada uno de ellos se generan 11 puntos alrededor de la recta 
// y otros 11 alrededor de la parábola
 
auto hTestSTD_H0 = new TH1F("hTestSTD_H0","",200,0,50); // histograma donde voy a guardar los chi2 asumiendo H0 verdadera
auto hTestSTD_H1 = new TH1F("hTestSTD_H1","",200,0,50); // histograma donde voy a guardar los chi2 asumiendo H1 verdadera

for (int iexp=0; iexp<Nexp; iexp++){
  for (int i=0; i<npoints; i++) {
      double x = xx[i];
      ey[i] = 0.5;
      yyL[i] = r.Gaus (      b*x,ey[i]);  // data de H0: lineal
      yyP[i] = r.Gaus (a*x*x+b*x,ey[i]);  // data de H1: parabola
   }

   auto grL = new TGraphErrors(npoints,xx,yyL,ex,ey);
   auto grP = new TGraphErrors(npoints,xx,yyP,ex,ey);

   float chi2_H0_Linear = grL->Chisquare(fLinear);  // Calculo chi2 respecto de una recta los puntos generados a partir de la recta
   float chi2_H1_Linear = grP->Chisquare(fLinear);  // Calculo chi2 respecto de recta a los puntos generados a partir de la parabola

   hTestSTD_H0 -> Fill(chi2_H0_Linear); // voy llenando un histograma con los chi2 que resultan cuando H0 es cierta
   hTestSTD_H1 -> Fill(chi2_H1_Linear); // voy llenando un histograma con los chi2 que resultan cuando H1 es cierta

}

hTestSTD_H0->SetLineWidth(2);
hTestSTD_H0->SetLineColor(kBlue);

hTestSTD_H1->SetLineWidth(2);
hTestSTD_H1->SetLineColor(2);
   
/// Calcular Quantiles para STD /////////////////////////////////

TH1F* h;
double STDerrT1, STDerrT2;
//double LLRerrT1, LLRerrT2;
double xq[1] = {0.95};  // Obtener el quantil 0.95
double yq[1];           // array donde guardarlo

/// Calcular Quantiles para STD ////////////////////////////////

h = hTestSTD_H0;
// Calculo el error tipo 1
h->GetQuantiles(1,yq,xq); // Busco el cuantil en la distribución obtenida bajo H0 verdadera
float QuantileSTD = yq[0];
// Integral a izquierda de la distribución cuando H0 es verdadera
STDerrT1 = h->Integral(h->FindBin(   0),h->FindBin(QuantileSTD),""); 

h = hTestSTD_H1;
// Calculo el error tipo 2
// La integral a izquierda de la distribución cuando H1 es verdadera (usando el mismo cuantil que antes!)
STDerrT2 = h->Integral(h->FindBin(   0),h->FindBin(QuantileSTD),""); 

printf("\n STD   Corte:%+7.2f   CL:%5.2f    Error T2:%6.3f \n",
       QuantileSTD,STDerrT1/Nexp,STDerrT2/Nexp);

////////// Dibujar STD //////////////////////////////////////////

auto  XaxisSTD = hTestSTD_H0->GetXaxis();
int   QuantilBinSTD = XaxisSTD->FindBin(QuantileSTD);
float QuantilModSTD = XaxisSTD->GetBinUpEdge(QuantilBinSTD);

// Crear histogramas de los errores Tipo1 y Tipo2

auto RejRegionSTD_H0 = (TH1F*)hTestSTD_H0->Clone("RejRegionH0");
RejRegionSTD_H0->GetXaxis()->SetRangeUser(QuantilModSTD,50);
RejRegionSTD_H0->SetFillStyle(3335);
RejRegionSTD_H0->SetFillColor(4);

auto RejRegionSTD_H1 = (TH1F*)hTestSTD_H1->Clone("RejRegionH1");
RejRegionSTD_H1->GetXaxis()->SetRangeUser(0,QuantilModSTD);
RejRegionSTD_H1->SetFillStyle(3353);
RejRegionSTD_H1->SetFillColor(2);

// Dibujar los histogramas de STD

c1->cd(2);
gPad->SetGrid();

hTestSTD_H0->Draw();
c1->Print("Simple_Hipotesis.pdf"); // print canvas
hTestSTD_H1->Draw("same");
c1->Print("Simple_Hipotesis.pdf"); // print canvas
RejRegionSTD_H0->Draw("same");
c1->Print("Simple_Hipotesis.pdf"); // print canvas
RejRegionSTD_H1->Draw("same");
c1->Print("Simple_Hipotesis.pdf"); // print canvas

// Dibujar la leyenda

gStyle->SetLegendBorderSize(0);
gStyle->SetLegendFillColor(0);
gStyle->SetLegendFont(42);
gStyle->SetLegendTextSize(0.05);

auto legend1 = new TLegend(0.7,0.75,0.9,0.89);
legend1->AddEntry(RejRegionSTD_H0,"Test STD|H0","f");
legend1->AddEntry(RejRegionSTD_H1,"Test STD|H1","f");
legend1->Draw();
c1->Print("Simple_Hipotesis.pdf"); // 'print canvas

// Guarda en el pdf y lo cierra.
c1->Print("Simple_Hipotesis.pdf");  // print canvas
c1->Print("Simple_Hipotesis.pdf]"); // close file

// Guarda histogramas para hacer las curvas ROC
TFile *f = new TFile("Histos.root","recreate");
hTestSTD_H0->Write();
hTestSTD_H1->Write();

f->Close();

}
