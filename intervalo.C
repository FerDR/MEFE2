double promedio(vector<double> datos,int n=10){
  double suma=0;
  for(int i=0;i<n;i++)suma = suma+datos[i];
  suma = suma/n;
  return suma;
}

double mierda(vector<double> datos,int n=10){//sumas de las raices n+1-esimas del dato n
  //perdon pero la que sugirieron con 10 datos agarra valores ~10^13, se complica
  double total=0;
  for(int i=0;i<n;i++)total = total+pow(datos[i],1.0/(i+1));
  return total;
}

double wilks(vector<double> datos,int n=10){
  double Like = 0;
  return Like;
}

void cinturon(double (*func)(vector<double>,int),string name, int n=10,int N=100000,bool debug=false){
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
  cout.precision(4);
  cout<<fixed;

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
      stat[toy] = func(x,n);

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
  //TODO
  //Chequear cobertura
  //Hacer lo de -2LnLog(Like/OptLike)
}

void intervalo(){
  cinturon(&promedio,string("promedio"));
  cinturon(&mierda,string("mierda"));
}
