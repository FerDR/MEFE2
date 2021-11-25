/*
 Test Wilks para Gamma, con alfa fijo  LLR = -2*[LL(beta)-LL(betaHat)] es chisq_1 ?
 la formula para beta_hat analitica es betahat = alfa/x_mean
*/

void Wilks_Gamma_1D(int n=5, double beta=0.5, double alpha=3.0, int nToys = 10000) 
{

//Set style defaults
gStyle->SetOptStat(10);      // only histo name and Nentries
gStyle->SetOptFit(1111);     // display fit results
gStyle->SetPadTickX(1);      // right ticks also
gStyle->SetPadTickY(1);      // upper ticks also
gStyle->SetFuncWidth(3);     // thicker function lines
gStyle->SetHistLineWidth(2); // thickrr histo lines

TRandom3 r(0);          // Initialize random generator
vector<double> x(n);    // Will store the n random variables
cout.precision(4);      // Print using four decimals precision
cout<<fixed;            // ...exactly four decimals

// x[n] is fixed and passed as [&].
// val[0] is beta, the independent variable, There are no parameters.
auto LL = [&](double *val, double *par){ 
   double sum = 0; 
   double betaP = val[0];
   for (auto x_i : x) sum += log(ROOT::Math::gamma_pdf (x_i,alpha,1./betaP));
   return sum;
};

// TF1 LL: last 3 arguments xmin,xmax,npar
// Here xmin=betamin, xmax=betamax, npar=0
auto fLL = new TF1("LL",LL,0.001,8*beta,0);
fLL->SetTitle("Log-Likelihood LL(x|#alpha,#beta) for gamma pdf;beta;LL");

// Gamma function, to generate random numbers
 auto fGamma = new TF1("fGamma","ROOT::Math::gamma_pdf(x,[0],1./[1])",0.001,4*alpha/beta);
 fGamma->SetParameter(0,alpha);
 fGamma->SetParameter(1,beta);

// Statistics and histos to fill for each experiment
auto hbetahat =   new TH1D("hbeta1"    ,"Distribucion de beta estimator",200,0,4*beta);
auto hLLR =       new TH1D("hLLR"      ,"Distribucion de LLR",200,0,10);
auto hLLbeta =    new TH1D("hLLbeta"   ,"Distribucion de -2logL(beta) y -2logL(betahat) (azul/rojo)",100,0,n*10);
auto hLLbetahat = new TH1D("hLLbetahat","Distribucion de -2logL(beta) y -2logL(betahat) (azul/rojo)",100,0,n*10);
auto hPvalue =    new TH1D("hPvalue"   ,"Distribucion de p-value",100,0,1);

auto c1 = new TCanvas("c1");
c1->Draw();

/////////////////////////////////////
//                                 //
//     Loop over experiments       //
//                                 //
/////////////////////////////////////

for (int iexp = 0; iexp < nToys; ++iexp) {

    // Generate n gamma(alfa,beta) random numbers
    for (int i = 0; i < n; ++i) x[i] = fGamma->GetRandom();

    // Wilks: distribution of -2log[LL(beta)/LL(betahat)] is chi2(ndf=1) ?
    double x_sum = 0;
    for (auto x_i : x)  x_sum += x_i;
    double betahat = alpha * n/x_sum;
    double LLbeta = -2*fLL->Eval(beta);
    double LLbetahat = -2*fLL->Eval(betahat);
    double LLR = LLbeta - LLbetahat ;
    double pvalue = 1 - ROOT::Math::chisquared_cdf(LLR,1);
    hLLR->Fill(LLR);
    hLLbeta->Fill(LLbeta);
    hLLbetahat->Fill(LLbetahat);
    hPvalue->Fill(pvalue);
    hbetahat->Fill(betahat);
}                                                                                   

hbetahat->SetLineColor(kBlue);
if (n>5) hbetahat->Fit("gaus");
hbetahat->Draw();
gPad->Update();
gPad->WaitPrimitive();

hLLbeta->Draw("E");

hLLbetahat->SetLineColor(kRed);
hLLbetahat->Draw("E same");
gPad->Update();
gPad->WaitPrimitive();

gPad->SetLogy(1);
hLLR->Scale(1./nToys);
hLLR->Draw("E");
gPad->Update();
gPad->WaitPrimitive();

auto fchi2 = new TF1("chi2","(10./200)*ROOT::Math::chisquared_pdf(x,1)");
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


/*
La distribucion gamma es:
f(x|alfa,beta) = [beta^alfa/Gamma(alfa)] x^(alfa-1) e^(-beta*x)

donde la funcion Gamma(x) en root se calcula como ROOT::Math::tgamma(x)

Tomemos como alfa=3 y beta=0.5

Sea X una VA con distribucion f(x|alfa,beta), y defino una nueva variable Y=-2log[f(X|alfa,beta)],
la distribucion de Y es extranha, no baja de 4 y diverge en 4 viniendo de la derecha, hLLbeta

Esto se entiende analizando al funcion -2log[f(X|alfa,beta):

TF1 *f = new TF1("f","-2*([a]*log([b])-log(ROOT::Math::tgamma([a]))+([a]-1)*log(x)-[b]*x)",0.1,10);
     f->SetParameters(3,0.5);
     f->Draw();

donde se ve que a bajo x domina ([a]-1)*log(x) y a alto x domina -[b]*x
*/


// Transform input CL => number of sigmas => Delta(logL)
// double nsigmas = ROOT::Math::gaussian_quantile((1+CL)/2,1);
// double delta = nsigmas*nsigmas/2;  
//
// To Draw logL for each experiment
// auto fLLDraw = new TF1("log L(beta|xi)",LL,0.001,8*beta,0);
// fLLDraw->SetTitle(";Parameter #beta;log L(#beta|x)");

