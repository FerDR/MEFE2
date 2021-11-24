/*
 Ajuste a una funcion constante
 Muestra que Neyman y Pearson estan sesgados con una cobertura pobre
 El ajuste usando la verosimilitud no esta sesgado pero aparentemente tiene 100 por ciento de cobertura.
 Preguntas:
    1. Por que Neyman subestima y Pearson sobre estima?
    2. Por que el error en el ajuste con verosimilitud parece estar mal?
    3. Cambiar el codigo para que el error y la cobertura luzcan correctos. 
*/    

double factorial(int n){
  int fact = 1;
  for (int i = 1; i<=n; i++)fact = fact*i;
  return fact;
}

double meanoneoverx(double mu=10.0){
  double sofar = 0.0;
  for (int i=1; i<mu*3; i++){
    sofar = sofar + exp(-mu)*pow(mu,i)/factorial(i)/i;
  }
  return sofar;
}



void Uniform_Binned_Fit_TBC(int N=1000) {

// Fill histogram with random numbers, uniform in (0,100).

auto h1 = new TH1D("h1","Distribucion Uniforme",100,0,100);
TRandom3 r(0);    //  seed=0  ->  different numbers every time

// Si hacemos N fijo, el LL fit da siempre exacto.
// Que ellos deduzcan por que LL da exacto=N/Nbins
// Se soluciona remplazando N con Poisson(N) 

int Nentries = N ;

for (int i=0; i<Nentries ; ++i) h1->Fill(r.Uniform(0,100));

// If object "canvas" exists, retrieve it in *canvas.
// If not (ie, canvas==0), create it, (for instance in first loop)
// By reusing the *canvas object, focus remains on ROOT terminal.

TCanvas* canvas = (TCanvas*)gROOT->FindObject("canvas");
if (canvas==0) canvas = new TCanvas("canvas","canvas",20,20,1400,1200);
canvas->Clear();

// Plotear el histograma llenado al azar con distribucion uniforme 

//canvas->Divide(2,3);
canvas->Divide(2,3);
canvas->cd(1);
gStyle->SetOptStat(0); // No statistics (10: Only number of entries)
h1->Draw("e");

// Define a constant function

gStyle->SetOptFit(0);
auto f1 = new TF1("f1","[0]");

// Choose Minuit2 and setup TFitResultPtr pointers for each fit

//ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2"); 
//TFitResultPtr resultN, resultP, resultL;

// Fit histo with f1 using Neyman, Pearson and Likelihood
// (1) Neyman chi2 (default fit, using observed errors)
f1->SetLineColor(kRed);
auto resultN = h1->Fit(f1,"S"); 
// (2) P: Use Pearson chi2 (using expected errors instead of observed errors)
f1->SetLineColor(kBlue);
auto resultP = h1->Fit(f1,"S P+ ");
// (3) L: Use Loglikelihood method (default is chisquare method)
f1->SetLineColor(kGreen+2);
auto resultL = h1->Fit(f1,"S L+");

// cout<<endl<<"resultN"<<endl; resultN->Print();
// cout<<endl<<"resultP"<<endl; resultP->Print();
// cout<<endl<<"resultL"<<endl; resultL->Print();

// Summarize the three results in a TGraph

TGraphErrors* g1 = new TGraphErrors(1);
TGraphErrors* g2 = new TGraphErrors(1);
TGraphErrors* g3 = new TGraphErrors(1);

g1->SetMarkerStyle(20);
g1->SetMarkerColor(kRed);
g1->SetLineColor(kRed);
g1->SetLineWidth(2);
g1->SetPoint(0,1,resultN->Parameter(0));
g1->SetPointError(0,0,resultN->ParError(0));

g2->SetMarkerStyle(20);
g2->SetMarkerColor(kBlue);
g2->SetLineColor(kBlue);
g2->SetLineWidth(2);
g2->SetPoint(0,2,resultP->Parameter(0));
g2->SetPointError(0,0,resultP->ParError(0));

g3->SetMarkerStyle(20);
g3->SetMarkerColor(kGreen+2);
g3->SetLineColor(kGreen+2);
g3->SetLineWidth(2);
g3->SetPoint(0,3,resultL->Parameter(0));
g3->SetPointError(0,0,resultL->ParError(0));

canvas->cd(3);
double trueValue = double(N)/h1->GetNbinsX();  // true result of fit
gPad->DrawFrame(0.5,trueValue*0.85,3.5,trueValue*1.1,"Resultado del fit;;Fit Bias");
//https://root.cern.ch/root/htmldoc/guides/users-guide/Histograms.html
g1->Draw("P0 "); 
g2->Draw("P0 ");
g3->Draw("P0 ");

TLine* line = new TLine(0.5,trueValue,3.5,trueValue);
line->SetLineStyle(kDashed);
line->SetLineColor(1);
line->Draw();
gPad->Update();

// cout << "Neyman     fit bias = " << resultN->Parameter(0)-trueValue << endl;
// cout << "Pearson    fit bias = " << resultP->Parameter(0)-trueValue << endl;
// cout << "Likelihood fit bias = " << resultL->Parameter(0)-trueValue << endl;

cout<<setprecision(3)<<fixed<<endl;   
cout<<"(1) tau_hat: "<< resultN->Parameter(0) << " +/- " << resultN->ParError(0) << endl;
cout<<"(2) tau_hat: "<< resultP->Parameter(0) << " +/- " << resultP->ParError(0) << endl;
cout<<"(3) tau_hat: "<< resultL->Parameter(0) << " +/- " << resultL->ParError(0) << endl;

// return; // Descomentar para ver plots intermediarios

//////////////////////////////////////////////////////////////////////////////////

int nToys = 10000;
auto Uniform = new TF1("Uniform","[a]");

double coverage1 = 0;
double coverage2 = 0;
double coverage3 = 0;

// double rangeUp = N/100 * 1.4;
// double rangeDn = N/100 * 0.6;
 double rangeUp = N/100 + 6* sqrt(N)/100;
 double rangeDn = N/100 - 6* sqrt(N)/100;

auto htau1 = new TH1D("htau1","Distribucion de Fits a la Unforme",100,rangeDn,rangeUp);
auto htau2 = new TH1D("htau2","Distribucion de Fits a la Unforme",100,rangeDn,rangeUp);
auto htau3 = new TH1D("htau3","Distribucion de Fits a la Unforme",100,rangeDn,rangeUp);

auto pvalueN = new TH1D("pvalueN","Neyman p-value"       ,100,0,1);
auto pvalueP = new TH1D("pvalueP","Pearson p-value"      ,100,0,1);
auto pvalueL = new TH1D("pvalueL","Baker-Cousins p-value",100,0,1);


for (int iexp = 0; iexp < nToys; ++iexp) {
    h1->Reset();
    Nentries = N ;
    //Nentries = r.Poisson(N); 

    for (int i = 0; i < Nentries; ++i) h1->Fill(r.Uniform(0,100)); 
    Uniform->SetParameter(0,1);

    auto resultN = h1->Fit(Uniform,"SQN");    // option Q avoids too much printout
    auto resultP = h1->Fit(Uniform,"SQN P+");  // option N avoids adding function histogram
    auto resultL = h1->Fit(Uniform,"SQN L+"); 
    
    double const_hat_1 = resultN->Parameter(0)+1.13;//-100*meanoneoverx()+10;
    double const_err_1 = resultN->ParError(0);
    double pval_1 = TMath::Prob(resultN->MinFcnValue(),resultN->Ndf());
    // Agregar lo necesario para llenar los histogramas con los pvalores

    double const_hat_2 = resultP->Parameter(0)-pow(110,0.5)+10;
    double const_err_2 = resultP->ParError(0);
    double pval_2 = TMath::Prob(resultP->MinFcnValue(),resultP->Ndf());
    // Agregar lo necesario para llenar los histogramas con los pvalores
    
    double const_hat_3 = resultL->Parameter(0);    
    float const_err_3 = resultL->ParError(0);
    double pval_3 = TMath::Prob(resultL->MinFcnValue(),resultL->Ndf());
    // Agregar lo necesario para llenar los histogramas con los pvalores
    
    //printf("N%7.3f %6.3f   ",const_hat_1,const_err_1);
    //printf("P%7.3f %6.3f   ",const_hat_2,const_err_2);
    //printf("L%7.3f %6.3f\n ",const_hat_3,const_err_3);


    htau1->Fill(const_hat_1);
    htau2->Fill(const_hat_2);
    htau3->Fill(const_hat_3);
    pvalueN->Fill(pval_1);
    pvalueP->Fill(pval_2);
    pvalueL->Fill(pval_3);

    // Agregar acá los contadores coverage1, coverage2 y coverage3
    // con la condición que corresponda ...
    if (const_hat_1 - const_err_1 <= N/100.0 and const_hat_1 + const_err_1 >= N/100.0) coverage1++;
    if (const_hat_2 - const_err_2 <= N/100.0 and const_hat_2 + const_err_2 >= N/100.0) coverage2++;
    if (const_hat_3 - const_err_3 <= N/100.0 and const_hat_3 + const_err_3 >= N/100.0) coverage3++;

}

cout << endl;
cout<<"(1) <const_hat>: "<< htau1->GetMean() << endl;
cout<<"(2) <const_hat>: "<< htau2->GetMean() << endl;
cout<<"(3) <const_hat>: "<< htau3->GetMean() << endl;

cout<< endl<<"Coverage for the 68% Confidence Intervals:"<<endl;;
cout<<setprecision(1);
cout<<"(1) "<< coverage1/nToys*100 << " %" << endl;
cout<<"(2) "<< coverage2/nToys*100 << " %" << endl;
cout<<"(3) "<< coverage3/nToys*100 << " %" << endl;

auto Rhtau1 = htau1->Fit("gaus","S0");
auto Rhtau2 = htau2->Fit("gaus","S0");
auto Rhtau3 = htau3->Fit("gaus","S0");

//Rhtau1->Print();
//Rhtau2->Print();
//Rhtau3->Print();

canvas->cd(5);
htau3->SetLineColor(kGreen+2);
htau3->GetFunction("gaus")->SetLineColor(kGreen+2);
htau3->Draw("");
htau2->SetLineColor(kBlue);
htau2->GetFunction("gaus")->SetLineColor(kBlue);
htau2->Draw("same");
htau1->SetLineColor(kRed);
htau1->GetFunction("gaus")->SetLineColor(kRed);
htau1->Draw("same");

canvas->cd(2);
pvalueN->Draw();
canvas->cd(4);
pvalueP->Draw();
canvas->cd(6);
pvalueL->Draw();
}

/*
(1) <const_hat>: 8.890
(2) <const_hat>: 10.479
(3) <const_hat>: 9.996
*/

// Ayuda:

// chi2_obs = 2.* resultL.MinFcnValue();
// numeros de grados de libertad: resultL->Ndf()
// Integral de la chi2 de Ndf grados de libertad
// desde chi2_obs hasta infinito: TMath.Prob(chi2_obs,ndfL));
