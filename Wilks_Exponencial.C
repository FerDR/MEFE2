void Wilks_Exponencial(int n=5, double tau=2, int nToys = 10000) 
{
// int n=5;
// double tau=2;

//Set style defaults
gStyle->SetOptStat(0);      // only histo name and Nentries
gStyle->SetOptFit(0);     // display fit results
gStyle->SetPadTickX(1);      // right ticks also
gStyle->SetPadTickY(1);      // upper ticks also
gStyle->SetFuncWidth(3);     // thicker function lines
gStyle->SetHistLineWidth(2); // thickrr histo lines

TRandom3 r(0);          // Initialize random generator
vector<double> x(n);    // Will store the n random variables
cout.precision(4);      // Print using four decimals precision
cout<<fixed;            // ...exactly four decimals

// Log Likelihood function for exponential: LL(tau)
// tau = val[0], la variable independiente
// par = vacia, no hay parametros
// Los valores de x se acceden por [&] (por ee necesito una lambda function

auto LL = [&](double *val, double *par){ 
   double sum = 0; 
   double tau = val[0]; 
   for (auto x_i : x)  sum += log(1./tau * exp(-x_i/tau));
   return sum; 
};

// Y la cargo en un TF1. last 3 arguments xmin,xmax,npar(=0)
auto fLL = new TF1("LL",LL,0.001,8*tau,0);
fLL->SetTitle("Log-Likelihood LL(tau|x) for exponential case;tau;LL");

// To Draw logL for each experiment
// auto fLLDraw = new TF1("log L(tau|xi)",LL,0.001,8*tau,0);
// fLLDraw->SetTitle(";Parameter #tau;log L(#tau|x)");

// Exponential function. Used by TTree::UnbinnedFit
// auto fExpo = new TF1("fExpo", "(1/[0])*exp(-x/[0])",0.001,8*tau);
// fExpo->SetParameter(0,tau);

// Transform input CL => number of sigmas => Delta(logL)
// double nsigmas = ROOT::Math::gaussian_quantile((1+CL)/2,1);
// double delta = nsigmas*nsigmas/2;  

// Statistics and histos to fill for each experiment
auto htau =      new TH1D("htau1"    ,"Distribution of estimated tau",200,0,3*tau);
auto hLLR =      new TH1D("hLLR"     ,"Distribution of -2log(LLR)",200,0,8);
auto hLLtau =    new TH1D("hLLtau"   ,"Distribution of -2log(LLtau)",100,0,n*6);
auto hLLtauhat = new TH1D("hLLtauhat","Distribution of -2log(LLtauhat)",100,0,n*6);
auto hPvalue =   new TH1D("hPvalue"  ,"Distribution of pvalue de -2log(LLR)",100,0,1);

auto c1 = new TCanvas("c1");
c1->Draw();

/////////////////////////////////////
//                                 //
//     Loop over experiments       //
//                                 //
/////////////////////////////////////

for (int iexp = 0; iexp < nToys; ++iexp) {

    // Generate n exponential randon numbers
    for (int i = 0; i < n; ++i) x[i] = r.Exp(tau);

    // Check Wilks: distribution of -2log[LL(tau)/LL(tauhat)] is chi2(ndf=1) ?
    double x_sum = 0;
    for (auto x_i : x)  x_sum += x_i;
    double tauhat = x_sum/n;
    double LLtau = -2*fLL->Eval(tau);
    double LLtauhat = -2*fLL->Eval(tauhat);
    double LLR = -2*fLL->Eval(tau) + 2*fLL->Eval(tauhat);
    double pvalue = 1-ROOT::Math::chisquared_cdf(LLR,1);
    hLLR->Fill(LLR);
    hLLtau->Fill(LLtau);
    hLLtauhat->Fill(LLtauhat);
    hPvalue->Fill(pvalue);
    htau->Fill(tauhat);
}                                                                                   

// Distribucion de hat_tau
htau->SetLineColor(kBlue);
if (n>5) htau->Fit("gaus");
htau->Draw();
gPad->Update();
gPad->WaitPrimitive();

// Distribucion de -2lnL(x|hat_tau)
hLLtauhat->SetLineColor(kRed);
hLLtauhat->Draw("E");
gPad->Update();
gPad->WaitPrimitive();

// Distribucion de -2lnL(x|tau)
hLLtau->Draw("E same ");
gPad->Update();
gPad->WaitPrimitive();

// Distribucion de -2logLLR: [-2logL(x|tau)] - [-2logL(x|hat_tau)]. Es una chi2_1?
gPad->SetLogy(1);
hLLR->Scale(1./nToys);
hLLR->Draw("E");
auto fchi2 = new TF1("chi2","(8./200)*ROOT::Math::chisquared_pdf(x,1)");
fchi2->SetRange(0,10);
fchi2->Draw("same");
gPad->Update();
gPad->WaitPrimitive();

// Distribucion de p-values -2logLLR suponiendo una chi2_1
gPad->SetLogy(0);
hPvalue->SetMinimum(0);
hPvalue->Draw("E");
hPvalue->Fit("pol0");
gPad->Update();

}



// Opcional en la efinicion de la funcion lambda:
//for (auto x_i : x) sum += log(ROOT::Math::exponential_pdf (x_i,1/tau));
