double promedio(vector<double> datos,int n=10){
  double suma=0;
  for(int i=0;i<n;i++)suma = suma+datos[i];
  suma = suma/n;
  return suma;
}

double mierda(vector<double> datos,int n=10){
  double total=0;
  for(int i=0;i<n;i++)total = total+pow(datos[i],1.0/(i+1));
  return total;
}


void cinturon(double (*func)(vector<double>,int),string name, int n=10,int N=10000,bool debug=true){
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
  //auto hMierda = new TH1F("hMierda","",2000,8,40);
  double Highs[100];
  double Lows[100];
  //double MierdaHighs[100];
  //double MierdaLows[100];
  int iTau = 0;
  double Taus[100];
  for(double tau=0.1;tau<=10;tau=tau+0.1){
    Taus[iTau] = tau;
    //hPromedio->Reset();
    //hMierda->Reset();
    
    for(int toy = 0; toy<N; toy++){  
      double promedio = 0;
      //double mierda = 0;
      for(int i = 0; i < n; i++){
        x[i] = r.Exp(tau);
	if(debug&&N<20)cout<<x[i]<<endl;
	//promedio = promedio + x[i];
        //mierda = mierda + pow(x[i],1/(i+1));
      }
      //promedio = promedio/n;
      stat[toy] = func(x,n);
      if(debug&&N<20)cout<<stat[toy]<<endl;

    }
    //Mas lento pero asi tengo bin edges dinamicos
    std::sort(stat.begin(),stat.end());
    auto h = new TH1F("h","",2000,stat[0]-0.1,stat[N-1]+0.1);
    for(int iStat = 0; iStat<N;iStat++)h->Fill(stat[iStat]);
    //hMierda->Fill(mierda);
    if(debug){
      c1->Clear();
      //c1->Divide(1,2);
      //c1->cd(1);
      h->Draw();
      //c1->cd(2);
      //hMierda->Draw();
      c1->Print(pdfName.c_str(),"pdf");
    }
    //Find best interval
    auto c = h->GetCumulative();
    //auto cMierda = hMierda->GetCumulative();
    c->Scale(1/h->Integral());
    //cMierda->Scale(1/hMierda->Integral());
    if(debug){
      c1->Clear();
      //c1->Divide(1,2);
      //c1->cd(1);
      c->Draw();
      //c1->cd(2);
      //cMierda->Draw();
      c1->Print(pdfName.c_str(),"pdf");
    }
    int intDist = 1e6;
    double lowEdge;
    double highEdge;
    //double lowEdgeMier;
    //double highEdgeMier;
    double binContent;
    double mintofind;
    int tempDist;
    int lowBin;
    for (int iBin = 1; iBin<c->GetNbinsX();iBin++){
      binContent = c->GetBinContent(iBin);
      if (binContent<0.68)continue;
      mintofind = binContent-0.68;
      for (int jBin = 1; jBin<iBin;jBin++){
	if (c->GetBinContent(jBin)>mintofind){
	  if (abs(c->GetBinContent(jBin)-mintofind)<abs(c->GetBinContent(jBin-1))){
	    tempDist = iBin-jBin;
	    lowBin = jBin;
	  }
	  else {
	    tempDist = iBin-jBin+1;
	    lowBin = jBin-1;
	  }

	break;
	}
      }
      if (tempDist<intDist){
	intDist = tempDist;
        lowEdge = c->GetBinCenter(lowBin);
        highEdge = c->GetBinCenter(iBin);
      }
    }
    /*intDist = 1e6;
    //if(debug)std::cout<<"Promedio: El borde bajo es "<<lowEdgeProm<<" y el alto es "<<highEdgeProm<<std::endl;
    for (int iBin = 1; iBin<cMierda->GetNbinsX();iBin++){
      binContent = cMierda->GetBinContent(iBin);
      if (binContent<0.68)continue;
      mintofind = binContent-0.68;
      for (int jBin = 1; jBin<iBin;jBin++){
	if (cMierda->GetBinContent(jBin)>mintofind){
	  if (abs(cMierda->GetBinContent(jBin)-mintofind)<abs(cMierda->GetBinContent(jBin-1))){
	    tempDist = iBin-jBin;
	    lowBin = jBin;
	  }
	  else {
	    tempDist = iBin-jBin+1;
	    lowBin = jBin-1;
	  }
	break;
	}
      }
      if (tempDist<intDist){
	intDist = tempDist;
        lowEdgeMier = cMierda->GetBinCenter(lowBin);
        highEdgeMier = cMierda->GetBinCenter(iBin);
      }
    }*/
    //if(debug)std::cout<<"Mierda: El borde bajo es "<<lowEdgeMier<<" y el alto es "<<highEdgeMier<<std::endl;
    Lows[iTau] = lowEdge;
    Highs[iTau] = highEdge;
    //MierdaLows[iTau] = lowEdgeMier;
    //MierdaHighs[iTau] = highEdgeMier;
    iTau++;
    //if(debug)break;
  }
  auto gHighs = new TGraph(100,Highs,Taus);
  auto gLows = new TGraph(100,Lows,Taus);
  //auto gMierdaHighs = new TGraph(100,MierdaHighs,Taus);
  //auto gMierdaLows = new TGraph(100,MierdaLows,Taus);
  auto axis = gHighs->GetXaxis();
  //auto axisMier = gMierdaHighs->GetXaxis();
  axis->SetLimits(Lows[0]-0.5,Highs[99]+0.5);
  //axisMier->SetLimits(MierdaLows[0]-0.5,MierdaHighs[99]+0.5);
  c1->Clear();
  //c1->Divide(1,2);
  //c1->cd(1);
  gHighs->SetTitle(";Estadistico;#Tau");
  gHighs->Draw();
  gLows->Draw("same");
  //c1->cd(2);
  //gMierdaHighs->SetTitle(";Mierda;#Tau");
  //gMierdaHighs->Draw();
  //gMierdaLows->Draw("same");
  c1->Print(pdfEnd.c_str(),"pdf");
 //TODO
 //Chequear cobertura
 //Hacer lo de -2LnLog(Like/OptLike)

}

void intervalo(){
  cinturon(&mierda,string("mierda"));
}
