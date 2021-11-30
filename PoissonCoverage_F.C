// Este script calcula el intervalo frecuentista para el mu de una variable aleatoria de Poisson.

void PoissonCoverage_F(double mu_min = 0.0, double mu_max = 6.5, double CL = 0.6827,  int nscan_points= 1000, float CoverageMin = 0.0) {
	
// Transformamos el CL en un número de sigmas de la Gaussiana N(0, 1)
// esto es de utilidad para el cálculo de intervalos usando la varianza y la verosimilitud (ambos pendientes de implementar)

double nsigmas = ROOT::Math::gaussian_quantile((1+CL)/2,1);
double delta = nsigmas*nsigmas/2; 

// (1) Intervalo usando la varianza:  nobs +- sqrt(nobs) /////////////////////
auto IsInside1 = [&](double nobs, double mu){
    return (nobs-nsigmas*pow(nobs,0.5)<=mu && mu<=nobs+nsigmas*pow(nobs,0.5));
};

// (2) Intervalo usando la verosimilitud: LL(ci_limits) = LLmax-1/2 ///////////
auto IsInside2 = [&](double nobs, double mu){
    auto LL = new TF1("LL","[0]*log(x)-x",0,100.01);
    LL->SetParameter(0,nobs);
    double LLMax = LL->GetMaximum();
    double muMin = (nobs==0)? 0 : LL->GetX(LLMax-delta,1e-4,nobs);
    double muMax = (nobs==0)? delta : LL->GetX(LLMax-delta,nobs+1e-4,100);
    return (muMin <= mu && mu <= muMax);
};

// (3) Intervalo frecuentista (analítico) ///////////////////////////////////
auto IsInside3 = [&](double nobs, double mu) { 
    double alpha = 1.- CL;
    double xmin = (nobs >0) ? 0.5*ROOT::Math::chisquared_quantile( alpha/2, 2*nobs) : 0;
    double xmax = 0.5* ROOT::Math::chisquared_quantile_c( alpha/2, 2*(nobs+1));
    return (xmin <= mu && mu <= xmax); 
};


// (4) Intervalo bayesiano ////////////////////////////
auto IsInside4 = [&](double nobs, double mu) {
    auto posterior = new TF1("posterior","[0]*TMath::PoissonI([1],x)",0,1000);
    posterior->SetParameter(1,nobs);
    posterior->SetParameter(0,1);
    double norm = posterior->Integral(0,1000);
    posterior->SetParameter(0,1./norm);
    double xmin;
    double xmax;
    double muMax = posterior->GetMaximumX(0,1000);
    //Solucion exacta para el caso n=0 para asegurarme que el minimo del intervalo no es el maximo de la funcion
    if(nobs==0){
      xmin=0;
      xmax=-log(1-CL);
    }
    else{
      double thisx;
      double thisy;
      double thisintegral;
      double prevx;
      double prevy;
      double previntegral;
      for(int i = 0; i < 100; i++){
	thisx = muMax/100*(i+1);
	thisy = posterior->Eval(thisx);
	thisintegral = posterior->Integral(thisx,posterior->GetX(thisy,muMax,1000));
	if(thisintegral>CL){
	  prevx = thisx;
	  prevy = thisy;
	  previntegral = thisintegral;
	  continue;
	}
	else{
	  xmin = (abs(thisintegral-CL) < abs(previntegral-CL)) ? prevx : thisx;
	  xmax = (abs(thisintegral-CL) < abs(previntegral-CL)) ? posterior->GetX(prevy,muMax,1000) : posterior->GetX(thisy,muMax,1000);
	  break;
        }
      }
    }
    return (xmin <= mu && mu <= xmax);
};


////////////////////////////////////////////////////////////////////////
// Cobertura vs mu a ser guardadas en tres objetos de la clase TGraphs ////
TGraph * g1 = new TGraph(nscan_points);
TGraph * g2 = new TGraph(nscan_points);
TGraph * g3 = new TGraph(nscan_points);
TGraph * g4 = new TGraph(nscan_points);

// Bucle sobre mu, desde mu_min a mu_max.
// Para cada mu calcula la cobertura del intervalo y lo guarda en el TGraphs.
for (int i = 0; i < nscan_points; ++i) {
   double mu = mu_min + (mu_max-mu_min)/(nscan_points-1) * i;
   //cout<<i<<" "<<mu<<" "<<mu_min<<" "<<mu_max<<endl;
   
   int Nmax = int(mu+10*sqrt(mu))+1; // esperanza mas 10 veces sigma

   // Inicializar las variables que faltan
   double probInside1 = 0;
   double probInside2 = 0;
   double probInside3 = 0;
   double probInside4 = 0;

   // Barro desde n observado hasta la esperanza mas 10 veces sigma
   for (int nobs = 0; nobs < Nmax; ++nobs) {

      // Completar con lo correspondiente a los otros dos intervalos
      bool inside1 = IsInside1(nobs,mu);
      bool inside2 = IsInside2(nobs,mu);
      bool inside3 = IsInside3(nobs,mu);
      bool inside4 = IsInside4(nobs,mu);

      // Calculo la probabilidad de obtener nobs de una Poisson con media mu
      double prob = ROOT::Math::poisson_pdf(nobs,mu); 
      // Si nobs esta dentro del intervalo, agrega la probabilidd de ese caso particular, i.e. Poisson(nobs,mu).
      if (inside1) probInside1 += prob;
      if (inside2) probInside2 += prob;
      if (inside3) probInside3 += prob;
      if (inside4) probInside4 += prob;
   }
   
   g1->SetPoint(i,mu,probInside1);
   g2->SetPoint(i,mu,probInside2);
   g3->SetPoint(i,mu,probInside3);
   g4->SetPoint(i,mu,probInside4);
   
}

////////////////////////////////////////////////////////////////////////
// Create TCanvas & TFrame, and a TLine at CL for reference ////////////
TCanvas * canvas = new TCanvas("canvas","canvas",700,500);

gPad->DrawFrame(mu_min,CoverageMin,mu_max,1,"Comparacion de los Intervalos de Confianza para Poisson;\\text{Parametro de Poisson }\\mu;Cobertura");

auto l = new TLine(mu_min,CL,mu_max,CL);
l->SetLineStyle(kDashed);
l->Draw(); // Draw a reference line at y-axis = CL

// The 3 TGraphs can be plotted together (I added a little offset between them),
// but first better look at each of them one at a time. Uncomment as necessary:
 
// Naive Poisson interval (Red) ////////////////////////////////////////
  g1->SetLineWidth(2);
  g1->SetLineColor(kRed);
  g1->Draw("L");
  gPad->Update();
  gPad->WaitPrimitive();

// LogLikelihood interval (Blue line) //////////////////////////////////
  g2->SetLineWidth(2);
  g2->SetLineColor(kBlue);
  g2->Draw("L");
  gPad->Update();
  gPad->WaitPrimitive();
  
// Analytic frequentist interval (Black line) //////////////////////////
  g3->SetLineWidth(2);
  g3->SetLineColor(kBlack);
  g3->Draw("L");
  gPad->Update();
  gPad->WaitPrimitive();

  g4->SetLineWidth(2);
  g4->SetLineColor(kGreen);
  g4->Draw("L");

  auto leg = new TLegend();
  leg->AddEntry(g1,"Varianza");
  leg->AddEntry(g2,"LogLikelihood");
  leg->AddEntry(g3,"Frecuentista");
  leg->AddEntry(g4,"Bayes");
  leg->Draw("same");
  gPad->Update();

}
