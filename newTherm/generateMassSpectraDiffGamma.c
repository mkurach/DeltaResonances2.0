//#include "ParticleType.h"
//#include "ParticleDB.h"
//#include "Parser.h"
#include "moje.h"
#include <string>
#include <fstream>
#include <sstream>
#include <vector>

#define _N_PARTICLES_ 15


ParticleDB *particlesDB;
DecayTable* decaysTables[_N_PARTICLES_];
ParticleType *particles[_N_PARTICLES_];
 Int_t pid[_N_PARTICLES_] = {12112, 1214, 22112, 32214, 2122, 32212, 2216, 12216,5218, 2128,32124,22124, 9299,2218,12214};//12214
 //Int_t pid[_N_PARTICLES_] = {12112};
 int colors[]={1, 600, 629, 414, 802, 880, 819, 922,433,618};



Double_t momentum(Double_t m, Double_t m1, Double_t m2) {
    if(m != 0 && m >= (m1 + m2))
        return TMath::Sqrt(TMath::Abs((m*m - (m1+m2)*(m1+m2))*(m*m - (m1-m2)*(m1-m2))/(4.0*m*m)));
    else 
        return 0;
}


Double_t GammaK2Stable(Double_t m,ParticleType* mother,ParticleType* kid1, ParticleType* kid2,Double_t branchingRatio) {
  //Double_t czynnik = mother->GetMass() / m;
  Double_t czynnik = 1.0;
  Double_t Ltr = TMath::Abs(TMath::Abs(mother->GetSpin())-(TMath::Abs(kid1->GetSpin()) + TMath::Abs(kid2->GetSpin()))); 
  Double_t pierwszyNawias;
  if(momentum(mother->GetMass(),kid1->GetMass(),kid2->GetMass()) == 0)
    pierwszyNawias = 0;
  else  
   pierwszyNawias = TMath::Power( momentum(m,kid1->GetMass(),kid2->GetMass()) / momentum(mother->GetMass(),kid1->GetMass(),kid2->GetMass()) , 2.0*Ltr+1);

  Double_t beta = TMath::Power(mother->GetMass() - kid1->GetMass() - kid2->GetMass(),2.0) + TMath::Power(mother->GetGamma(),2.0)/4.0;
  Double_t drugiNawias;
  if((beta*beta + TMath::Power(momentum(m,kid1->GetMass(),kid2->GetMass()),2.0)) == 0)
    drugiNawias = 0;
  else
    drugiNawias = (beta*beta + TMath::Power(momentum(mother->GetMass(),kid1->GetMass(),kid2->GetMass()),2.0))/(beta*beta + TMath::Power(momentum(m,kid1->GetMass(),kid2->GetMass()),2.0));
  Double_t gammaK = branchingRatio * mother->GetGamma();

  if (mother->GetMass() > (kid1->GetMass() + kid2->GetMass() ))
    return czynnik * pierwszyNawias * drugiNawias *gammaK;
  else  
    return 0;

}



Double_t BWwithGamma(Double_t *x, Double_t *par) { // 0 - normalizacja, 1 - index czastki

    Int_t index =(Int_t)par[1];
   // DecayChannel *channel;
    ParticleType* part1;
     ParticleType* part2;
     Double_t gammaTot = 0;
    for(int i = 0; i <= decaysTables[index]->GetChannelCount(); i++) { //pętla po kanałach rozpadów
      //channel = decaysTables[index]->GetDecayChannel(i);
      part1 = particlesDB->GetParticleType(decaysTables[index]->GetDecayChannel(i)->GetParticle1());
      part2 = particlesDB->GetParticleType(decaysTables[index]->GetDecayChannel(i)->GetParticle2());
      gammaTot += GammaK2Stable(x[0],particles[index],part1,part2,decaysTables[index]->GetDecayChannel(i)->GetBranchingRatio());
    }
    Double_t poleMass = particles[index]->GetMass();
    if (((poleMass*poleMass-x[0]*x[0])*(poleMass*poleMass-x[0]*x[0]) + x[0]*x[0]*gammaTot*gammaTot) == 0 )
      return 0;
    else 
      return par[0]*x[0]*x[0]*gammaTot/((poleMass*poleMass-x[0]*x[0])*(poleMass*poleMass-x[0]*x[0]) + x[0]*x[0]*gammaTot*gammaTot);

}


void generateMassSpectra() {






    //std::vector<const DecayChannel*> decaysChannels[_N_PARTICLES_];



    Double_t mMin[_N_PARTICLES_];
    Double_t mMax[_N_PARTICLES_];
    Double_t gamma, mass;
   
    TF1* fun[_N_PARTICLES_]; 
    TCanvas *can[_N_PARTICLES_];
    TH1D *hist[_N_PARTICLES_];

    Double_t massSum, massTresh[_N_PARTICLES_];

    for(int i = 0; i < _N_PARTICLES_; i++) {

            particles[i] = particlesDB->FindByPID(pid[i]);
      mMin[i] = particles[i]->GetMass() - 2.0*particles[i]->GetGamma();
      mMax[i] = particles[i]->GetMass() + 12.0*particles[i]->GetGamma();

  cout<<particles[i]->GetName()<<endl;
    //  cout<<particle[i]->GetMass()<<endl;
      decaysTables[i] = particles[i]->GetTable();
      fun[i] = new TF1(particles[i]->GetName(),BWwithGamma,mMin[i],mMax[i],2);
      fun[i]->SetParameter(0,1);
      fun[i]->SetParameter(1,i);
      

      can[i]= new TCanvas(Form("can%s",particles[i]->GetName()),Form("can%s",particles[i]->GetName()),1000,1000);
     can[i]->cd();
      fun[i]->Draw();
      hist[i] = (TH1D*)fun[i]->GetHistogram();
      //hist[i]->Set
        //cout<<"min wartosc: "<<mMin[i]<<endl;
      //cout<<"gamma: "<<particle[i]->GetGamma()<<endl;
      //cout<<"szerokosz binu: "<<hist[i]->GetBinWidth(10)<<endl;
      //hist[i]->Draw();
      //for(int j= 0; j < 10; j++) 
        //  hist[i]->SetBinContent(j+1,0);
    
      //for(int j= 0; j < hist[i]->GetNbinsX(); j++) 
         // cout<<hist[i]->GetBinContent(j+1)<<endl;
    
      //cout<<hist[i]->GetNbinsX()<<endl;
      //cout<<hist[i]->GetBinContent(10)<<endl;
      can[i]->SaveAs(Form("outputDiffGamma/%sMassDiffGamma.png",particles[i]->GetName()));
   

    }

        //PRINT HISTOGRAMS
    /*ofstream fileTxt;
    fileTxt.open("histogramsDiffGamma.txt");


    fileTxt<<"TH1F *getBreitWigner(int pdg) {"<<endl;
    fileTxt<<"\tstatic std::map<int, std::vector<Double_t>> histograms;"<<endl;
    fileTxt<<"\tstatic std::map<int, Double_t> mMin;"<<endl;
    fileTxt<<"\tstatic std::map<int, Double_t> mMax;"<<endl;
    for(int i=0; i < _N_PARTICLES_; i++) {
        
        
        fileTxt<<"\thistograms["<<pid[i]<<"] = {";
        for(int j = 1; j < hist[i]->GetNbinsX()+1; j++) {
            if (j != hist[i]->GetNbinsX())
                fileTxt<<hist[i]->GetBinContent(j)<<",";
            else    
                fileTxt<<hist[i]->GetBinContent(j)<<"};"<<endl;
        }
        fileTxt<<"\tmMin["<<pid[i]<<"] = "<<mMin[i]<<";"<<endl;
        fileTxt<<"\tmMax["<<pid[i]<<"] = "<<mMax[i]<<";"<<endl;

    }
    fileTxt<<"\tstatic TH1F *hist = new TH1F(\"hist\",\"hist\",1000,mMin[pdg],mMax[pdg]);"<<endl;
    fileTxt<<"\tfor(int i = 1; i < hist->GetNbinsX(); i++)"<<endl;
    fileTxt<<"\t\thist->SetBinContent(i,histograms[pdg][i-1]);"<<endl;
    fileTxt<<endl<<"\treturn hist;"<<endl<<"}"<<endl<<endl;
    fileTxt.close();
*/

}

//Double_t m,ParticleType* mother,ParticleType* kid1, ParticleType* kid2,Double_t branchingRatio
Double_t GammaKForDrawing(Double_t *x, Double_t *par) {//0 - normalizacja, 1 - index czastki matki, 2 - kanal rozpadu

  ParticleType* mother = particles[(int)par[1]];
  ParticleType* kid1 = particlesDB->GetParticleType(decaysTables[(int)par[1]]->GetDecayChannel((int)par[2])->GetParticle1());
  ParticleType* kid2 = particlesDB->GetParticleType(decaysTables[(int)par[1]]->GetDecayChannel((int)par[2])->GetParticle2());
  Double_t branchingRatio = decaysTables[(int)par[1]]->GetDecayChannel((int)par[2])->GetBranchingRatio();
  //cout<<"Branching ratio "<<branchingRatio<<endl;
  //Double_t czynnik = mother->GetMass() / x[0];
  Double_t czynnik = 1.0;
  Double_t Ltr = TMath::Abs(TMath::Abs(mother->GetSpin())-(TMath::Abs(kid1->GetSpin()) + TMath::Abs(kid2->GetSpin()))); 
  Double_t pierwszyNawias;
  if(momentum(mother->GetMass(),kid1->GetMass(),kid2->GetMass()) == 0)
    pierwszyNawias = 0;
  else  
   pierwszyNawias = TMath::Power( momentum(x[0],kid1->GetMass(),kid2->GetMass()) / momentum(mother->GetMass(),kid1->GetMass(),kid2->GetMass()) , 2.0*Ltr+1);

  Double_t beta = TMath::Power(mother->GetMass() - kid1->GetMass() - kid2->GetMass(),2.0) + TMath::Power(mother->GetGamma(),2.0)/4.0;
  Double_t drugiNawias;
  if((beta*beta + TMath::Power(momentum(x[0],kid1->GetMass(),kid2->GetMass()),2.0)) == 0)
    drugiNawias = 0;
  else
    drugiNawias = (beta*beta + TMath::Power(momentum(mother->GetMass(),kid1->GetMass(),kid2->GetMass()),2.0))/(beta*beta + TMath::Power(momentum(x[0],kid1->GetMass(),kid2->GetMass()),2.0));
  Double_t gammaK = branchingRatio * mother->GetGamma();

  if (mother->GetMass() > (kid1->GetMass() + kid2->GetMass() ))
    return par[0]* czynnik * pierwszyNawias * drugiNawias *gammaK;
  else  
    return 0;

}

Double_t BreitWignerKForDrawing(Double_t *x, Double_t *par) {//0 - normalizacja, 1 - index czastki matki, 2 - kanal rozpadu

  ParticleType* mother = particles[(int)par[1]];
  ParticleType* kid1 = particlesDB->GetParticleType(decaysTables[(int)par[1]]->GetDecayChannel((int)par[2])->GetParticle1());
  ParticleType* kid2 = particlesDB->GetParticleType(decaysTables[(int)par[1]]->GetDecayChannel((int)par[2])->GetParticle2());
  Double_t branchingRatio = decaysTables[(int)par[1]]->GetDecayChannel((int)par[2])->GetBranchingRatio();
  //cout<<"Branching ratio "<<branchingRatio<<endl;
  //Double_t czynnik = mother->GetMass() / x[0];
  Double_t czynnik = 1.0;
  Double_t Ltr = TMath::Abs(TMath::Abs(mother->GetSpin())-(TMath::Abs(kid1->GetSpin()) + TMath::Abs(kid2->GetSpin()))); 
  Double_t pierwszyNawias;
  if(momentum(mother->GetMass(),kid1->GetMass(),kid2->GetMass()) == 0)
    pierwszyNawias = 0;
  else  
   pierwszyNawias = TMath::Power( momentum(x[0],kid1->GetMass(),kid2->GetMass()) / momentum(mother->GetMass(),kid1->GetMass(),kid2->GetMass()) , 2.0*Ltr+1);

  Double_t beta = TMath::Power(mother->GetMass() - kid1->GetMass() - kid2->GetMass(),2.0) + TMath::Power(mother->GetGamma(),2.0)/4.0;
  Double_t drugiNawias;
  if((beta*beta + TMath::Power(momentum(x[0],kid1->GetMass(),kid2->GetMass()),2.0)) == 0)
    drugiNawias = 0;
  else
    drugiNawias = (beta*beta + TMath::Power(momentum(mother->GetMass(),kid1->GetMass(),kid2->GetMass()),2.0))/(beta*beta + TMath::Power(momentum(x[0],kid1->GetMass(),kid2->GetMass()),2.0));
  Double_t gammaK = branchingRatio * mother->GetGamma();
  Double_t gammaTot = 0;
  if (mother->GetMass() > (kid1->GetMass() + kid2->GetMass() ))
    gammaTot =  czynnik * pierwszyNawias * drugiNawias *gammaK;


  Double_t poleMass = mother->GetMass();
    if (((poleMass*poleMass-x[0]*x[0])*(poleMass*poleMass-x[0]*x[0]) + x[0]*x[0]*gammaTot*gammaTot) == 0 )
      return 0;
    else 
      return par[0]*x[0]*x[0]*gammaTot/((poleMass*poleMass-x[0]*x[0])*(poleMass*poleMass-x[0]*x[0]) + x[0]*x[0]*gammaTot*gammaTot);



}

Double_t GammaTot(Double_t *x, Double_t *par) { // 0 - normalizacja, 1 - index czastki

    Int_t index =(Int_t)par[1];
   // DecayChannel *channel;
    ParticleType* part1;
     ParticleType* part2;
     Double_t gammaTot = 0;
    for(int i = 0; i < decaysTables[index]->GetChannelCount(); i++) { //pętla po kanałach rozpadów
      //channel = decaysTables[index]->GetDecayChannel(i);
      part1 = particlesDB->GetParticleType(decaysTables[index]->GetDecayChannel(i)->GetParticle1());
      part2 = particlesDB->GetParticleType(decaysTables[index]->GetDecayChannel(i)->GetParticle2());
      gammaTot += GammaK2Stable(x[0],particles[index],part1,part2,decaysTables[index]->GetDecayChannel(i)->GetBranchingRatio());
    }
    return gammaTot;

}

void drawGamma() {
  std::vector<TF1*> gammaChannels[_N_PARTICLES_];
  std::vector<TF1*> breitWignerChannels[_N_PARTICLES_];
  TF1* BWTotal[_N_PARTICLES_];
  TF1* GammaTotal[_N_PARTICLES_];

  TCanvas *can[_N_PARTICLES_];
  TCanvas *canBW[_N_PARTICLES_];
  TFile *fileOut = new TFile("outputDiffGamma/diffGamma.root","RECREATE");
  TFile *fileOutBW = new TFile("outputDiffGamma/diffGammaBW.root","RECREATE");
  TString name1,name2;
  
    Double_t mMin[_N_PARTICLES_];
    Double_t mMax[_N_PARTICLES_];
    Double_t mTresh;
    TH1D *hist[_N_PARTICLES_];


  for(int i = 0; i < _N_PARTICLES_; i++) {
      particles[i] = particlesDB->FindByPID(pid[i]);
      mMin[i] = particles[i]->GetMass() - 2.0*particles[i]->GetGamma();
      mMax[i] = particles[i]->GetMass() + 12.0*particles[i]->GetGamma();
      decaysTables[i] = particles[i]->GetTable();
      can[i]= new TCanvas(Form("can%s",particles[i]->GetName()),Form("can%s",particles[i]->GetName()),1000,1000);
      canBW[i]= new TCanvas(Form("can%sBW",particles[i]->GetName()),Form("can%sBW",particles[i]->GetName()),1000,1000);
      //cout<<particle[i]->GetGamma()<<endl;
      cout<<particles[i]->GetName()<<endl;
      cout<<"Liczba kanalow "<<decaysTables[i]->GetChannelCount()<<endl;

        /*BWTotal[i] = new TF1(particle[i]->GetName(),BWwithGamma,mMin[i],mMax[i],2);
        BWTotal[i]->SetParameter(0,1);
        BWTotal[i]->SetParameter(1,i);
        canBW[i]->cd();
        BWTotal[i]->Draw();*/

        GammaTotal[i] = new TF1(particles[i]->GetName(),GammaTot,mMin[i],mMax[i],2);
        GammaTotal[i]->SetParameter(0,1);
        GammaTotal[i]->SetParameter(1,i);
        GammaTotal[i]->SetMinimum(0);
        can[i]->cd();
        //GammaTotal[i]->Draw();

    for(int j = 0; j <= decaysTables[i]->GetChannelCount(); j++) {
      name1 = particlesDB->GetParticleType(decaysTables[i]->GetDecayChannel(j)->GetParticle1())->GetName();
      name2 = particlesDB->GetParticleType(decaysTables[i]->GetDecayChannel(j)->GetParticle2())->GetName();
      cout<<Form("%s+%s",name1.Data(),name2.Data())<<endl;
      cout<<decaysTables[i]->GetDecayChannel(j)->GetBranchingRatio()<<endl;
      mTresh = particlesDB->GetParticleType(decaysTables[i]->GetDecayChannel(j)->GetParticle1())->GetMass() + particlesDB->GetParticleType(decaysTables[i]->GetDecayChannel(j)->GetParticle2())->GetMass();
      cout<<"Mtresh "<<mTresh<<endl;
      if(mTresh > mMin[i])
        mMin[i] = mTresh;
      
        gammaChannels[i].push_back(new TF1(Form("%s+%s",name1.Data(),name2.Data()),GammaKForDrawing,mMin[i],mMax[i],3));
        gammaChannels[i][j]->SetParameters(1.0,i,j);
        gammaChannels[i][j]->SetLineColor(colors[j]);
        gammaChannels[i][j]->SetLineStyle(2);
        gammaChannels[i][j]->SetMinimum(0);

        breitWignerChannels[i].push_back(new TF1(Form("%s+%s",name1.Data(),name2.Data()),BreitWignerKForDrawing,mMin[i],mMax[i],3));
        breitWignerChannels[i][j]->SetParameters(1.0,i,j);
        breitWignerChannels[i][j]->SetLineColor(colors[j]);
        breitWignerChannels[i][j]->SetLineStyle(2);
        breitWignerChannels[i][j]->SetMinimum(0);


        can[i]->cd();
       if(j==0)
          gammaChannels[i][j]->Draw();
        else
          gammaChannels[i][j]->Draw("same");
        canBW[i]->cd();
        if(j==0)
          breitWignerChannels[i][j]->Draw();
        else
          breitWignerChannels[i][j]->Draw("same");
        //hist[i] = (TH1D*)gammaChannels[i][j]->GetHistogram(); 

          
    }
    

    fileOut->cd();
    can[i]->Write();
    fileOutBW->cd();
    canBW[i]->Write();
    cout<<endl;


  }
  fileOut->Close();
  fileOut->Save();
    fileOutBW->Close();
  fileOutBW->Save();


}

void generateMassSpectraDiffGamma() {
    particlesDB = new ParticleDB();
    ReadSHAREParticles(particlesDB);
    ReadSHAREDecays(particlesDB);

  generateMassSpectra();
  //drawGamma();

}