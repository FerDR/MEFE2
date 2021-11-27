void intervalo(int n=10,int N=100000,bool debug=false){
  TRandom3 r(0);
  vector<double> x(n);
  gStyle->SetOptStat(0);      // only histo name and Nentries
  gStyle->SetPadTickX(1);      // right ticks also
  gStyle->SetPadTickY(1);      // upper ticks also
  gStyle->SetHistLineWidth(2); // thickrr histo lines
  cout.precision(4);
  cout<<fixed;

  auto c1 = new TCanvas("c1","Cinturon de confianza",700,1000);
  c1->Print("intervalo.pdf(","pdf");
  auto hPromedio = new TH1F("hPromedio","",2000,0,30);
  auto hMierda = new TH1F("hMierda","",2000,8,40);
  double PromedioHighs[100];
  double PromedioLows[100];
  double MierdaHighs[100];
  double MierdaLows[100];
  int iTau = 0;
  double Taus[100];
  for(double tau=0.1;tau<=10;tau=tau+0.1){
    Taus[iTau] = tau;
    hPromedio->Reset();
    hMierda->Reset();
    for(int toy = 0; toy<N; toy++){  
      double promedio = 0;
      double mierda = 0;
      for(int i = 0; i < n; i++){
        x[i] = r.Exp(tau);
        promedio = promedio + x[i];
        mierda = mierda + pow(x[i],1/(i+1));
      }
      promedio = promedio/n;
      hPromedio->Fill(promedio);
      hMierda->Fill(mierda);
    }
    if(debug){
      c1->Clear();
      c1->Divide(1,2);
      c1->cd(1);
      hPromedio->Draw();
      c1->cd(2);
      hMierda->Draw();
      c1->Print("intervalo.pdf","pdf");
    }
    //Find best interval
    auto cPromedio = hPromedio->GetCumulative();
    auto cMierda = hMierda->GetCumulative();
    cPromedio->Scale(1/hPromedio->Integral());
    cMierda->Scale(1/hMierda->Integral());
    if(debug){
      c1->Clear();
      c1->Divide(1,2);
      c1->cd(1);
      cPromedio->Draw();
      c1->cd(2);
      cMierda->Draw();
      c1->Print("intervalo.pdf","pdf");
    }
    int intDist = 1e6;
    double lowEdgeProm;
    double highEdgeProm;
    double lowEdgeMier;
    double highEdgeMier;
    double binContent;
    double mintofind;
    int tempDist;
    int lowBin;
    for (int iBin = 1; iBin<cPromedio->GetNbinsX();iBin++){
      binContent = cPromedio->GetBinContent(iBin);
      if (binContent<0.68)continue;
      mintofind = binContent-0.68;
      for (int jBin = 1; jBin<iBin;jBin++){
	if (cPromedio->GetBinContent(jBin)>mintofind){
	  if (abs(cPromedio->GetBinContent(jBin)-mintofind)<abs(cPromedio->GetBinContent(jBin-1))){
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
        lowEdgeProm = cPromedio->GetBinCenter(lowBin);
        highEdgeProm = cPromedio->GetBinCenter(iBin);
      }
    }
    intDist = 1e6;
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
    }
    //if(debug)std::cout<<"Mierda: El borde bajo es "<<lowEdgeMier<<" y el alto es "<<highEdgeMier<<std::endl;
    PromedioLows[iTau] = lowEdgeProm;
    PromedioHighs[iTau] = highEdgeProm;
    MierdaLows[iTau] = lowEdgeMier;
    MierdaHighs[iTau] = highEdgeMier;
    iTau++;
    //if(debug)break;
  }
  auto gPromedioHighs = new TGraph(100,PromedioHighs,Taus);
  auto gPromedioLows = new TGraph(100,PromedioLows,Taus);
  auto gMierdaHighs = new TGraph(100,MierdaHighs,Taus);
  auto gMierdaLows = new TGraph(100,MierdaLows,Taus);
  auto axisProm = gPromedioHighs->GetXaxis();
  auto axisMier = gMierdaHighs->GetXaxis();
  axisProm->SetLimits(PromedioLows[0]-0.5,PromedioHighs[99]+0.5);
  axisMier->SetLimits(MierdaLows[0]-0.5,MierdaHighs[99]+0.5);
  c1->Clear();
  c1->Divide(1,2);
  c1->cd(1);
  gPromedioHighs->SetTitle(";Promedio;#Tau");
  gPromedioHighs->Draw();
  gPromedioLows->Draw("same");
  c1->cd(2);
  gMierdaHighs->SetTitle(";Mierda;#Tau");
  gMierdaHighs->Draw();
  gMierdaLows->Draw("same");
  c1->Print("intervalo.pdf)","pdf");
 //TODO
 //Chequear cobertura
 //Hacer lo de -2LnLog(Like/OptLike)

}
