/*
 Test Wilks para Gamma:  LLR = -2*[LL(alfa,beta)-LL(alfaHat,betaHat)] es chisq_2 ?
 No hay formula cerrada para alfaHat, hay que obtener (alfaHat,betaHat) numericamente
*/

#include "TMinuit.h"  // needed for TTree::UnbinnedFit	

void Wilks_Gamma_2D(int n=10, double beta=0.5, double alfa=3.0, int nToys = 10000) 
{

gStyle->SetOptStat(0);      // only histo name and Nentries
gStyle->SetOptFit(0);       // display fit results
gStyle->SetPadTickX(1);      // right ticks also
gStyle->SetPadTickY(1);      // upper ticks also
gStyle->SetFuncWidth(3);     // thicker function lines
gStyle->SetHistLineWidth(2); // thickrr histo lines

TRandom3 r(0);          // Initialize random generator
vector<double> x(n);    // Will store the n random variables
cout.precision(4);      // Print using four decimals precision
cout<<fixed;            // ...exactly four decimals

// Log Likelihood LL(alfa,beta), x se pasa por [&] 

auto LL = [&](double alfa, double beta){ 
   double sum = 0; 
   for (auto x_i : x) sum += log(ROOT::Math::gamma_pdf (x_i,alfa,1./beta));
   return sum;
};

auto fGamma = new TF1("fGamma","ROOT::Math::gamma_pdf(x,[0],1./[1])",0.001,4*alfa/beta);
fGamma->SetParameter(0,alfa);
fGamma->SetParameter(1,beta);

// Statistics and histos to fill for each experiment
auto hbetaHat = new TH1D("hbetaHat","Estimador de beta",200,0,4*beta);
auto halfaHat = new TH1D("halfaHat","Estimador de alfa",200,0,4*alfa);
auto hLL0     = new TH1D("hLL0","-2log[L(x|alfa,beta)] y -2log[L(x|alfaHat,betaHat)] (azul/rojo)",100,0,n*10);
auto hLL1     = new TH1D("hLL1","-2log[L(x|alfa,beta)] y -2log[L(x|alfaHat,betaHat)] (azul/rojo)",100,0,n*10);
auto hLLR     = new TH1D("hLLR","LLR: -2*log[L(alfa,beta)/L(alfaHat,betaHat)]",200,0,10);
auto hPvalue  = new TH1D("hPvalue","Distribucion of p-value",100,0,1);

auto c1 = new TCanvas("c1");

/////////////////////////////////////
//                                 //
//      Loop over experiments      //
//                                 //
/////////////////////////////////////

for (int iexp = 0; iexp < nToys; ++iexp) {

    // Generate n gamma random numbers
    fGamma->SetParameter(0,alfa);
    fGamma->SetParameter(1,beta);
    for (int i = 0; i < n; ++i) x[i] = fGamma->GetRandom();

    // Save them on TTree
    double rndm; 
    TTree T("T","tree");
    T.Branch("rndm",&rndm,"rndm/D");
    for (int i = 0; i < n; ++i) { rndm=x[i]; T.Fill();}

    // Fit Tree to fGamma
    T.UnbinnedFit("fGamma","rndm","","Q");       // E:Minos Q:quiet
    double alfaHat = fGamma->GetParameter(0);
    double betaHat = fGamma->GetParameter(1);

    // Wilks: distribucion de -2log[LL(alfa,beta)/LL(alfaHat,betaHat)] es chi2(ndf=2) ?

    double LL0 = LL(alfa,beta);
    double LL1 = LL(alfaHat,betaHat);
    double LLR = (-2*LL0)-(-2*LL1);
    double pvalue = ROOT::Math::chisquared_cdf(LLR,2);
    
    hLLR->Fill(LLR);
    hLL0->Fill(-2*LL0);
    hLL1->Fill(-2*LL1);
    hPvalue->Fill(pvalue);
    halfaHat->Fill(alfaHat);
    hbetaHat->Fill(betaHat);
}                                                                                   

halfaHat->SetLineColor(kBlue);
if (n>5) halfaHat->Fit("gaus");
halfaHat->Draw();
gPad->Update();
gPad->WaitPrimitive();

hbetaHat->SetLineColor(kBlue);
if (n>5) hbetaHat->Fit("gaus");
hbetaHat->Draw();
gPad->Update();
gPad->WaitPrimitive();

hLL0->SetLineColor(kBlue);
hLL0->Draw("E");
hLL1->SetLineColor(kRed);
hLL1->Draw("E same");
gPad->Update();
gPad->WaitPrimitive();

gPad->SetLogy(0);
hLLR->Scale(1./nToys);
hLLR->Draw("E");
gPad->Update();
gPad->WaitPrimitive();

auto fchi2 = new TF1("chi2","(10./200)*ROOT::Math::chisquared_pdf(x,2)");
fchi2->SetRange(0,10);
fchi2->Draw("same");
gPad->Update();
gPad->WaitPrimitive();

gPad->SetLogy(0);
hPvalue->SetMinimum(0);
hPvalue->Draw("E");
hPvalue->Fit("pol0");
gPad->Update();


}




//double errorup,errordn,errparabolic,corr;
//gMinuit->mnerrs(0,errorup,errordn,errparabolic,corr);
//cout<<"TTree::UnbinnedFit : "<<alfaHat<<" "<<betaHat<<endl;
//    <<"+/-"<<fGamma->GetParError(0)
//	<<" ["<<errordn<<", +"<<errorup<<"]  "<<errparabolic<<endl;
