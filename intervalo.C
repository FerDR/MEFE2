double promedio(vector<double> datos,int n=10,double tau = 1){
  double suma=0;
  for(int i=0;i<n;i++)suma = suma+datos[i];
  suma = suma/n;
  return suma;
}

double mierda(vector<double> datos,int n=10,double tau = 1){//sumas de las raices n+1-esimas del dato n
  //perdon pero la que sugirieron con 10 datos agarra valores ~10^13, se complica
  double total=0;
  for(int i=0;i<n;i++)total = total+pow(datos[i],1.0/(i+1));
  return total;
}

double logLike(vector<double> datos,int n=10,double tau=1){
  double ll = 0;
  for(int i=0;i<n;i++){
    ll += -datos[i]/tau-log(tau);
  }
  return ll;
}

double wilks(vector<double> datos,int n=10,double tau = 1){
  double prom = promedio(datos,n,tau);
  return -2*logLike(datos,n,tau)+2*logLike(datos,n,prom);
}

void cinturon(double (*func)(vector<double>,int,double),string name,TGraph*& gUp, TGraph*& gDown, int n=10,int N=100000,bool debug=false){
  string pdfName = string("intervalo_");
  pdfName += name;
  pdfName += string(".pdf");
  string pdfBegin = pdfName;
  pdfBegin += string("(");
  string pdfEnd = pdfName;
  pdfEnd += string(")");
  TRandom3 r(0);
  vector<double> x(n);
  vector<double> stat(N);
  gStyle->SetOptStat(0);      // only histo name and Nentries
  gStyle->SetPadTickX(1);      // right ticks also
  gStyle->SetPadTickY(1);      // upper ticks also
  gStyle->SetHistLineWidth(2); // thickrr histo lines
  

  auto c1 = new TCanvas("c1","Cinturon de confianza",700,1000);
  c1->Print(pdfBegin.c_str(),"pdf");
  double Highs[100];
  double Lows[100];
  int iTau = 0;
  double Taus[100];
  for(double tau=0.1;tau<=10;tau=tau+0.1){
    Taus[iTau] = tau;
    for(int toy = 0; toy<N; toy++){  
      for(int i = 0; i < n; i++){
        x[i] = r.Exp(tau);
      }
      stat[toy] = func(x,n,tau);

    }
    //Mas lento pero asi tengo bin edges dinamicos
    std::sort(stat.begin(),stat.end());
    if(debug && N<20){
      for(int i=0;i<N;i++)cout<<stat[i]<<endl;
    }
    //El histograma esta SOLO para plotear, el intervalo se saca exacto
    auto h = new TH1F("h",name.c_str(),2000,stat[0]-0.1,stat[N-1]+0.1);
    for(int iStat = 0; iStat<N;iStat++)h->Fill(stat[iStat]);
    
    //Find best interval
    double lowEdge;
    double highEdge;
    int distance = N*68/100;
    double thisDist = 1e6;
    for(int iStat=distance;iStat<N-1;iStat++){
      double lowValue = stat[iStat-distance];
      double highValue = stat[iStat];
      if((highValue-lowValue)<thisDist){
	thisDist = highValue-lowValue;
	lowEdge = lowValue;
	highEdge = highValue;
      }
    }
    if(debug || r.Uniform(1)<0.05){
      auto Interval = (TH1F*)h->Clone("Interval");
      Interval->GetXaxis()->SetRangeUser(lowEdge,highEdge);
      Interval->SetFillStyle(1001);
      Interval->SetFillColor(2);
      string hTitle = string("Tau = ");
      hTitle += to_string(tau);
      hTitle += ";";
      hTitle += name;
      hTitle += ";Entries";
      h->SetTitle(hTitle.c_str());
      h->Draw("hist");
      Interval->Draw("same");
      c1->Print(pdfName.c_str(),"pdf");
      c1->Clear();
    }
    Lows[iTau] = lowEdge;
    Lows[iTau] = lowEdge;
    Highs[iTau] = highEdge;
    iTau++;
  }
  auto gHighs = new TGraph(100,Highs,Taus);
  auto gLows = new TGraph(100,Lows,Taus);
  auto axis = gHighs->GetXaxis();
  axis->SetLimits(Lows[0]-0.5,Highs[99]+0.5);
  string gTitle = string("Cinturon para estadistico ");
  gTitle += name;
  gTitle += string(";");
  gTitle += name;
  gTitle += string(";#tau");
  gHighs->SetTitle(gTitle.c_str());
  gHighs->Draw();
  gLows->Draw("same");
  c1->Print(pdfEnd.c_str(),"pdf");
  c1->Clear();
  gUp = gHighs;
  gDown = gLows;
  //TODO
  //Chequear cobertura
  //Hacer lo de -2LnLog(Like/OptLike)
}

void coverage(double (*func)(vector<double>,int,double),string name,TGraph* gUp,TGraph* gDown,int N=10000,int n=10){
  TRandom3 r(0);
  double taus[100];
  for (int i =0;i<100;i++)taus[i] = gUp->GetY()[i];
  vector<double> x(n);
  double coverages[100];
  double erx[100];
  double errores[100];
  double tau;
  double stat[N];
  double abajo;
  double arriba;
  for(int i=0;i<100;i++){
    coverages[i] = 0;
    tau = taus[i];
    for(int toy = 0; toy<N; toy++){  
      for(int i = 0; i < n; i++){
        x[i] = r.Exp(tau);
      }
      stat[toy] = func(x,n,tau);
      if(name!=string("wilks")){  
        abajo = min(gDown->Eval(stat[toy]),gUp->Eval(stat[toy]));
        arriba = max(gDown->Eval(stat[toy]),gUp->Eval(stat[toy]));
        if(abajo<tau && tau<arriba)coverages[i]++;
      }
      //Calculo ad-hoc despues de ver el intervalo porque si no la extrapolacion de TGraph.Eval() lo arruina
      else if(stat[toy]<=1)coverages[i]++;
    }
    erx[i] = 0;
    errores[i] = pow(coverages[i],0.5);
    errores[i] = errores[i]/N;
    coverages[i] = coverages[i]/N;
  }
  auto Cov_vs_tau = new TGraphErrors(100,taus,coverages,erx,errores);
  auto c2 = new TCanvas("c2","Cobertura vs tau",700,1000);
  Cov_vs_tau->Draw();
  string gName = string("Cobertura_");
  gName+=name;
  Cov_vs_tau->SetTitle(gName.c_str());
  gName+= string(".pdf");
  c2->Print(gName.c_str(),"pdf");
}

void intervalo(){
  cout.precision(4);
  cout<<fixed;
  gROOT->SetBatch(kTRUE);
  auto gUp = new TGraph(100);
  auto gDown = new TGraph(100);
  //cinturon(&promedio,string("promedio"),gUp,gDown);
  //coverage(&promedio,string("promedio"),gUp,gDown);
  //cinturon(&mierda,string("mierda"),gUp,gDown);
  //coverage(&mierda,string("mierda"),gUp,gDown);
  cinturon(&wilks,string("wilks"),gUp,gDown);
  coverage(&wilks,string("wilks"),gUp,gDown);
}
