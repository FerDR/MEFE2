/* Ajuste no bineado a n números aleatorios con distribución exponencial usando la verosimilitud
 * El estimador es siempre el de maxima verosimilitud pero y se comparan dos formas de estimar la incerteza 
 * dejando que un paquete de ROOT use la verosimilitud o hacerlo a mano.
 *
 *   (1) TTree::UnbinnedFit Likelihood  // FALTA AGREGAR OPCION PARA MAS DE 1 SIGMA
 *   (2) LL(ci) = LLmax-1/2*n^2 (Likelihood)
 *
*/

#include "TMinuit.h"  // needed for TTree::UnbinnedFit	

void Exponential_Unbinned_Fit_LL(int n=5, int nexperiments = 1000, double CL=0.9544) // 0.683 0.9544 0.9973
{

double tau=4; // exponential distribution parameter

// Cosmetic stuff //////////////////////////////////////////////////////
gStyle->SetOptStat(11);      // only histo name and Nentries
gStyle->SetOptFit(1111);     // display fit results
gStyle->SetPadTickX(1);      // right ticks also
gStyle->SetPadTickY(1);      // upper ticks also
gStyle->SetFuncWidth(3);     // thicker function lines
gStyle->SetHistLineWidth(2); // thickrr histo lines
////////////////////////////////////////////////////////////////////////

TRandom3 r(0);          // Initialize random generator
vector<double> x(n);    // Will store the n random variables
cout.precision(4);      // Print using four decimals precision
cout<<fixed;            // ...exactly four decimals

// Unbinned Log Likelihood LL(tau): 
// x[n] is fixed and passed as reference.
// p[0] is the independent variable. There are no parameters.
auto LL = [&](double * p, double *){ 
   double sum = 0; 
   double tau = p[0]; 
   for (auto x_i : x)  sum += log(1./tau * exp(-x_i/tau));
   return sum; 
};

// TF1 LL: ultimos tres argumentos xmin,xmax,npar
// Aca se define LL(tau), entonces: xmin=taumin, xmax=taumax, npar=0
auto fLL = new TF1("LL",LL,0.001,8*tau,0);
fLL->SetTitle("Log-Likelihood LL(tau|x) for exponential case;tau;LL");

// Funcion exponencial. Usada por TTree::UnbinnedFit
auto fExpo = new TF1("fExpo", "(1/[0])*exp(-x/[0])",0.001,8*tau);

// Transforma el Confidence Level en número de sigmas => Delta(logL)
double nsigmas = ROOT::Math::gaussian_quantile((1+CL)/2,1);
double delta = nsigmas*nsigmas/2;  

/////////////////////////////////////////////
//                                         //
//   Comienza el bucle sobre experimentos  //
//                                         //
/////////////////////////////////////////////

double coverage1 = 0;
double coverage2 = 0;

for (int iexp = 0; iexp < nexperiments; ++iexp) {

   // Generate n exponential randon numbers ///////////////////////////
   for (int i = 0; i < n; ++i) x[i] = r.Exp(tau);

   // (0) TTree::UnbinnedFit  //////////////////////////////////////////

   TTree T("T","tree");
   double rndm; 
   T.Branch("rndm",&rndm,"rndm/D");

   for (int i = 0; i < n; ++i) { rndm=x[i]; T.Fill();}
      
   fExpo->SetParameter(0,4);
   ROOT::Math::MinimizerOptions::SetDefaultErrorDef(delta);
   T.UnbinnedFit("fExpo","rndm","","EQ");       // E:Minos Q:quiet ////
       
   // Ahora le pedimos que calcule los errores por el método de la verosimilitud (errorup,errordn)
   // y también le pedimos el error usando la aproximación parabólica (errparabolic)
   double errorup,errordn,errparabolic,corr;
   gMinuit->mnerrs(0,errorup,errordn,errparabolic,corr);  
   gMinuit->SetPrintLevel(1);

   // (1) Intervalo por el metodo de la verosimilitud  /////////////////
   double tau1   = fExpo->GetParameter(0); // Tau estimado a partir de MLM;
   double ci1_dn = tau1+errordn;  // Limite inferior del intervalo
   double ci1_up = tau1+errorup;  // Limite superior del intervalo
       	
   // (2) LL(ci) = LLmax-delta (Likelihood) ///////////////////////////
   double tau2   = fLL->GetMaximumX();
   double LL_max = fLL->Eval(tau2); 
   double ci2_dn = fLL->GetX(LL_max-delta,0.0001,tau2);
   double ci2_up = fLL->GetX(LL_max-delta,tau2,5*tau2);

	cout.precision(2);

   cout<<"Resultado utilizando la verosimilitud:    "<<tau1<<" ["<<ci1_dn<<", +"<<ci1_up<<"]"<<endl<<endl;
   cout<<"Método de la verosimilitud hecho a mano:  "<<tau2<<" ["<<ci2_dn<<", +"<<ci2_up<<"]"<<endl<<endl;                                                
      
   if (tau > ci1_dn && tau < ci1_up) coverage1++; // Contamos cuantas veces el intervalo incluye el valor real
   if (tau > ci2_dn && tau < ci2_up) coverage2++;
    
}                                                                                   
                                                                                    
   cout<< endl<<"Cobertura del intervalo de confianza:"<<endl<<endl;;

   cout<<"Resultado utilizando la verosimilitud: :"<< coverage1/nexperiments*100 << " %" << endl;
   cout<<"Método de la verosimilitud hecho a mano:"<< coverage2/nexperiments*100 << " %" << endl<<endl;

}