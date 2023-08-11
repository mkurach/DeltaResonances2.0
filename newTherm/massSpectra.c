//#include "ParticleType.h"
//#include "ParticleDB.h"
//#include "Parser.h"
#include "moje.h"
#include <string>
#include <fstream>
#include <sstream>

#define _N_PARTICLES_ 8



void ReadSHAREParticles(ParticleDB* aDB) {
    std::ifstream mFile;
    mFile.open("particles.data");
    if((!mFile) && (!mFile.is_open())) 
        cout<<"Nie dziala"<<endl;
    //else
        //cout<<"dziala"<<endl;


  istringstream* iss;
  ParticleType*	 tPartBuf;
  char   buff[200];
  char   name[20];
  double mass, gamma, spin, I, I3, Nq, Ns, Naq, Nas, Nc, Nac, MC;
  int    number = 0;   

  while (!mFile.eof()) {
    mFile.getline(buff,200);
    if (!(*buff) || (*buff == '#'))
      continue;
    iss = new istringstream(buff);   
    (*iss) >> name >> mass >> gamma >> spin >> I >> I3 >> Nq >> Ns >> Naq >> Nas >> Nc >> Nac >> MC;
    number++;
    //cout<<'\t'<<number<<" "<<name<<" "<<mass<<" "<<gamma<<" "<<spin<<" "<<I<<" "<<I3<<" "<<Nq<<" "<<Naq<<" "<<Ns<<" - "<<static_cast<int>(Nas)<<" - "<<Nc<<" "<<Nac<<" "<<MC<<endl;
    //cout<<'\t'<<Form("%i",number)<<" "<<name<<" "<<mass<<" "<<gamma<<" "<<spin<<" "<<I<<" "<<I3<<" "<<Nq<<" "<<Naq<<" "<<Ns<<" "<<Nas<<" "<<Nc<<" "<<Nac<<" "<<MC);
    tPartBuf = new ParticleType();
    tPartBuf->SetNumber(number);
    tPartBuf->SetName(name);
    tPartBuf->SetMass(mass);
    tPartBuf->SetGamma(gamma);
    tPartBuf->SetSpin(spin);
    tPartBuf->SetBarionN(static_cast<int> ((Nq + Ns + Nc)/3. - (Naq + Nas + Nac)/3.) );
    tPartBuf->SetI(I);
    tPartBuf->SetI3(I3);
    tPartBuf->SetStrangeN(static_cast<int> (Nas - Ns));
    tPartBuf->SetCharmN(static_cast<int> (Nc - Nac));
    tPartBuf->SetNumberQ(static_cast<int> (Nq));
    tPartBuf->SetNumberAQ(static_cast<int> (Naq));
    tPartBuf->SetNumberS(static_cast<int> (Ns));
    tPartBuf->SetNumberAS(static_cast<int> (Nas));
    tPartBuf->SetNumberC(static_cast<int> (Nc));
    tPartBuf->SetNumberAC(static_cast<int> (Nac));
    tPartBuf->SetPDGCode(static_cast<int> (MC));
    aDB->AddParticleType(tPartBuf);
    delete iss;
  }
}


void ReadSHAREDecays(ParticleDB* aDB) {

        std::ifstream mFile;
    mFile.open("decays.data");
    if((!mFile) && (!mFile.is_open())) 
        cout<<"Nie dziala"<<endl;
    //else
        //cout<<"dziala"<<endl;

  istringstream* iss;
  char   buff[200];
  char   tFather[20], tDaughter1[20], tDaughter2[20], tDaughter3[20];
  double tBRatio, tRatio;
  int    CGcoeff=0; // complete branching ratio by Clebsch-Gordan coefficient: 0-no 1-yes

  while (!mFile.eof()) {
    mFile.getline(buff,200);
    if (!(*buff) || (*buff == '#'))
      continue;
    iss = new istringstream(buff);
    *iss >> tFather >> tDaughter1 >> tDaughter2 >> tDaughter3;
    if (!aDB->ExistsParticleType(tFather)) {
      cout<<"<Parser::ReadSHAREDecay>\tDid not find the father particle: " << tFather;
      cout<<"\tNot adding channel";
      delete iss;
      continue;
    }
    if (!aDB->ExistsParticleType(tDaughter1)) {
      cout<<"<Parser::ReadSHAREDecay>\tDid not find the daughter 1 particle: " << tDaughter1;
      cout<<"\tNot adding channel";
      delete iss;
      continue;
    }
    if (!aDB->ExistsParticleType(tDaughter2)) {
      cout<<"<Parser::ReadSHAREDecay>\tDid not find the daughter 2 particle: " << tDaughter2;
      cout<<"\tNot adding channel";
      delete iss;
      continue;
    }
    if ((*tDaughter3 > 65) && (*tDaughter3 < 122) && (!aDB->ExistsParticleType(tDaughter3))) {
      cout<<"<Parser::ReadSHAREDecay>\tDid not find the daughter 3 particle: " << tDaughter3;
      cout<<"\tNot adding channel";
      delete iss;
      continue;
    }
   // cout<<"\tDecay channel for "<<tFather;
    if ( (*tDaughter3 > 65) && (*tDaughter3 < 122) ) {
      // check if first char is a letter - if yes then 3-body decay
      *iss >> tBRatio >> CGcoeff;
      //cout<<"\t\tBR ("<< tBRatio <<")";
      if (aDB->GetParticleType(tDaughter1)->GetMass() + aDB->GetParticleType(tDaughter2)->GetMass() + aDB->GetParticleType(tDaughter3)->GetMass() < aDB->GetParticleType(tFather)->GetMass()) {
	DecayChannel* newChannel = new DecayChannel(tBRatio, aDB->GetParticleTypeIndex(tDaughter1), aDB->GetParticleTypeIndex(tDaughter2), aDB->GetParticleTypeIndex(tDaughter3));
	aDB->GetParticleType(tFather)->AddDecayChannel(*newChannel);
	aDB->GetParticleType(tFather)->SetDecayChannelCount3( aDB->GetParticleType(tFather)->GetDecayChannelCount3() + 1 );
	//cout<<"\t\tAdding 3-body decay channel:     "<< tDaughter1 <<" + "<< tDaughter2 <<" + "<< tDaughter3;
      } //else
        //cout<<"\t\tNOT adding 3-body decay channel: "<< tDaughter1 <<" + "<< tDaughter2 <<" + "<< tDaughter3 <<", extra mass = "<< aDB->GetParticleType(tDaughter1)->GetMass() + aDB->GetParticleType(tDaughter2)->GetMass() + aDB->GetParticleType(tDaughter3)->GetMass() - aDB->GetParticleType(tFather)->GetMass();
    } else {
    // 2-body decay
      tBRatio=atof(tDaughter3);
      *iss >> CGcoeff;
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      /*if (CGcoeff) {// complete branching ratio by Clebsch-Gordan coefficient
        double j1, m1, j2, m2, J, M, CB;
        J  = aDB->GetParticleType(tFather)->GetI();
        M  = aDB->GetParticleType(tFather)->GetI3();
        j1 = aDB->GetParticleType(tDaughter1)->GetI();
        m1 = aDB->GetParticleType(tDaughter1)->GetI3();
        j2 = aDB->GetParticleType(tDaughter2)->GetI();
        m2 = aDB->GetParticleType(tDaughter2)->GetI3();
        cout<<"\t\tCG coeff for: J("<<J<<") M("<<M<<") j1("<<j1<<") m1("<<m1<<") j2("<<j2<<") m2("<<m2<<")";
        CB = SHAREClebschGordan(J, M, j1, m1, j2, m2);
        tRatio = CB*CB * tBRatio;
        // Multiply the Clebsh by two?
        // The same spin, mass, strangeness (charm?)
        // and different I3?
        if(  (TMath::Abs(aDB->GetParticleType(tDaughter1)->GetSpin() - aDB->GetParticleType(tDaughter2)->GetSpin()) < 0.01) 
	  && (TMath::Abs(aDB->GetParticleType(tDaughter1)->GetMass() - aDB->GetParticleType(tDaughter2)->GetMass()) < 0.01)
	  && (TMath::Abs(aDB->GetParticleType(tDaughter1)->GetI3()   - aDB->GetParticleType(tDaughter2)->GetI3())   > 0.01)
	  && (aDB->GetParticleType(tDaughter1)->GetStrangeN()        - aDB->GetParticleType(tDaughter2)->GetStrangeN() == 0)
	  && (aDB->GetParticleType(tDaughter1)->GetCharmN()          - aDB->GetParticleType(tDaughter2)->GetCharmN() == 0)
	  ) {
          cout<<"\t\tMultiplying CG coeff by two";
          tRatio *= 2.0;
        }
        cout<<"\t\tCG("<< CB <<") BR("<< tBRatio <<") -> BR("<< tRatio<<")";
      } else {*/
        tRatio=tBRatio;
        //cout<<"\t\tBR("<< tBRatio <<")";
      //}
      if (aDB->GetParticleType(tDaughter1)->GetMass() + aDB->GetParticleType(tDaughter2)->GetMass() < aDB->GetParticleType(tFather)->GetMass()) {
	DecayChannel *newChannel = new DecayChannel(tRatio, aDB->GetParticleTypeIndex(tDaughter1), aDB->GetParticleTypeIndex(tDaughter2), -1);
        aDB->GetParticleType(tFather)->AddDecayChannel(*newChannel);
        aDB->GetParticleType(tFather)->SetDecayChannelCount2( aDB->GetParticleType(tFather)->GetDecayChannelCount2() + 1 );
        //cout<<"\t\tAdding 2-body decay channel:     "<< tDaughter1 <<" + "<< tDaughter2;
      } //else
        //cout<<"\t\tNOT adding 2-body decay channel: "<< tDaughter1 <<" + "<< tDaughter2 <<", extra mass = "<< aDB->GetParticleType(tDaughter1)->GetMass() + aDB->GetParticleType(tDaughter2)->GetMass() - aDB->GetParticleType(tFather)->GetMass();
    }
    delete iss;
  }
}


Double_t breitWigner(Double_t *x, Double_t *par) { //0 - normalization, 1 - Gamma, 2 - M_R
    return par[0]*x[0]*x[0]*par[1]/((par[2]*par[2]-x[0]*x[0])*(par[2]*par[2]-x[0]*x[0]) + x[0]*x[0]*par[1]*par[1]);
}


void createBWSpectra() { //set range

    Int_t pid[_N_PARTICLES_] = {12112, 1214, 22112, 32214, 2122, 32212, 2216, 12216};

    ParticleDB *particles = new ParticleDB();
    ReadSHAREParticles(particles);


    ParticleType *particle;
    Double_t mMin = 0.5;
    Double_t mMax = 3.0;
   
    TF1* fun[_N_PARTICLES_]; 
    TCanvas *can[_N_PARTICLES_];
    TH1D *hist[_N_PARTICLES_];
    TH1D *histTest[_N_PARTICLES_];

    for(int i = 0; i < _N_PARTICLES_; i++) {
        particle = particles->FindByPID(pid[i]);
        cout<<particle->GetName()<<endl;
        
        fun[i] = new TF1(Form("fun%s",particle->GetName()),breitWigner,mMin,mMax,3); //mass from 0-3, 3 parameters
        fun[i]->SetParameter(0,1); 
        fun[i]->SetParameter(1,particle->GetGamma());
        fun[i]->SetParameter(2,particle->GetMass());
        fun[i]->SetParameter(0,1./fun[i]->Integral(mMin,mMax));
        

        hist[i] = (TH1D*)fun[i]->GetHistogram();
        hist[i]->SetTitle(particle->GetName());
        histTest[i] = new TH1D(Form("histTest%s",particle->GetName()),Form("histTest%s",particle->GetName()),100,mMin,mMax);

        //cout<<hist[i]->GetNbinsX()<<endl;
        for(int j = 0; j < 10000; j ++){
            histTest[i]->Fill(hist[i]->GetRandom());
        }
        


        can[i]= new TCanvas(Form("can%s",particle->GetName()),Form("can%s",particle->GetName()),1000,1000);
        can[i]->cd();
        //fun[i]->Draw();
        //hist[i]->Draw();
        histTest[i]->SetLineColor(kBlue);
        histTest[i]->Draw();
        //can[i]->SaveAs(Form("%sTest.png",particle->GetName()));


    }

    //PRINT HISTOGRAMS
    ofstream fileTxt;
    fileTxt.open("histograms.txt");


    fileTxt<<"TH1F *getBreitWigner(int pdg) {"<<endl;
    fileTxt<<"\tstatic std::map<int, Double_t[]> histograms;"<<endl;
    for(int i=0; i < _N_PARTICLES_; i++) {
        
        
        fileTxt<<"\thistograms["<<pid[i]<<"] = {";
        for(int j = 1; j < hist[i]->GetNbinsX()+1; j++) {
            if (j != hist[i]->GetNbinsX())
                fileTxt<<hist[i]->GetBinContent(j)<<",";
            else    
                fileTxt<<hist[i]->GetBinContent(j)<<"};"<<endl;
        }

    }
    fileTxt<<"}"<<endl;
      //fileTxt<<"\tstatic TH1F *h = new TH1F(\""<<hist[i]->GetTitle()<<"\",\""<<hist[i]->GetTitle()<<"\",100,0.5,3.0);"<<endl;
       // fileTxt<<"\tfor(int i = 1; i < h"<<hist[i]->GetTitle()<<"->GetNbinsX(); i++)"<<endl;
        //fileTxt<<"\t\th"<<hist[i]->GetTitle()<<"->SetBinContent(i,y"<<hist[i]->GetTitle()<<"[i-1]);"<<endl;
        //fileTxt<<endl<<"\treturn h"<<hist[i]->GetTitle()<<";"<<endl<<"}"<<endl<<endl;
    fileTxt.close();








}


void createBWSpectraDiffRange() { 

    Int_t pid[_N_PARTICLES_] = {12112, 1214, 22112, 32214, 2122, 32212, 2216, 12216};

    ParticleDB *particles = new ParticleDB();
    ReadSHAREParticles(particles);
    ReadSHAREDecays(particles);

    DecayTable* decaysTables[_N_PARTICLES_];
    std::vector<const DecayChannel*> decaysChannels[_N_PARTICLES_];


    ParticleType *particle[_N_PARTICLES_];
    Double_t mMin[_N_PARTICLES_];
    Double_t mMax[_N_PARTICLES_];
    Double_t gamma, mass;
   
    TF1* fun[_N_PARTICLES_]; 
    TCanvas *can[_N_PARTICLES_];
    TH1D *hist[_N_PARTICLES_];
    TH1D *histTest[_N_PARTICLES_];

    ParticleType* part1;
    ParticleType* part2;
    ParticleType* part3;
    Double_t massSum, massTresh[_N_PARTICLES_];


    for(int i = 0; i < _N_PARTICLES_; i++) {
        particle[i] = particles->FindByPID(pid[i]);
        decaysTables[i] = particle[i]->GetTable();
       // cout<<particle[i]->GetName()<<endl;
        //cout<<decaysTables[i]->GetChannelCount()<<endl;
        for(int j = 0; j < decaysTables[i]->GetChannelCount(); j++) {
          decaysChannels[i].push_back(decaysTables[i]->GetDecayChannel(j));
          part1 = particles->GetParticleType(decaysChannels[i][j]->GetParticle1());
          part2 = particles->GetParticleType(decaysChannels[i][j]->GetParticle2());
          part3 = particles->GetParticleType(decaysChannels[i][j]->GetParticle3());
          //cout<<part1->GetName()<<" + "<<part2->GetName();
          massSum = part1->GetMass() + part2->GetMass();
          if (j == 0)
            massTresh[i] = massSum;
          else {
            if(massSum < massTresh[i])
              massTresh[i] = massSum;
          }
          //cout<<"\t"<<massSum<<endl;;
          //cout<<decaysChannels[i][j]->GetParticle1()-><<" + "<<decaysChannels[i][j]->GetParticle2()<<decaysChannels[i][j]->GetParticle3()<<endl;

        }
        //cout<<"\t"<<massTresh[i]<<endl<<endl;

        gamma = particle[i]->GetGamma();
        //cout<<gamma<<endl;
        mass = particle[i]->GetMass();
        //cout<<mass<<endl;
        mMin[i] = massTresh[i];
        mMax[i] = mass + 3.0*gamma;
        
        fun[i] = new TF1(Form("fun%s",particle[i]->GetName()),breitWigner,mMin[i],mMax[i],3); //mass from 0-3, 3 parameters
        fun[i]->SetParameter(0,1); 
        fun[i]->SetParameter(1,particle[i]->GetGamma());
        fun[i]->SetParameter(2,particle[i]->GetMass());
        fun[i]->SetParameter(0,1./fun[i]->Integral(mMin[i],mMax[i]));
        

        hist[i] = (TH1D*)fun[i]->GetHistogram();
        hist[i]->SetTitle(particle[i]->GetName());
        //cout<<hist[i]->Integral("width")<<endl;
        //cout<<hist[i]->GetNbinsX()<<endl;
        //cout<<mMin[i]<<endl;


        can[i]= new TCanvas(Form("can%s",particle[i]->GetName()),Form("can%s",particle[i]->GetName()),1000,1000);
        can[i]->cd();
        //fun[i]->Draw();
        hist[i]->Draw();
        can[i]->SaveAs(Form("outputMassSpectra/%sMassTresh.png",particle[i]->GetName()));


    }

    //PRINT HISTOGRAMS
    ofstream fileTxt;
    fileTxt.open("histograms.txt");


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
    fileTxt<<"\tstatic TH1F *hist = new TH1F(\"hist\",\"hist\",100,mMin[pdg],mMax[pdg]);"<<endl;
    fileTxt<<"\tfor(int i = 1; i < hist->GetNbinsX(); i++)"<<endl;
    fileTxt<<"\t\thist->SetBinContent(i,histograms[pdg][i-1]);"<<endl;
    fileTxt<<endl<<"\treturn hist;"<<endl<<"}"<<endl<<endl;
    fileTxt.close();








}

void test() {
    /*TF1 *f = ((TF1*)(gROOT->GetFunction("pol0")));
  // f->Print();
  f->SetRange(-100., 100.); // xmin, xmax
  f->SetParameter(0, 10.); // [0] = y
  TCanvas *can = new TCanvas("can","can",1000,1000);
  can->cd();
  f->Draw();
  can->SaveAs("test.png");*/

  TRandom *random = new TRandom();
  TH1D * hist = new TH1D("hist","hist",100,0,5);
  TCanvas *can = new TCanvas("can","can",1000,1000);
   for(int i = 0; i < 1000000; i ++)
    hist->Fill(random->Uniform(0,5));
  can->cd();
  hist->Draw();
  can->SaveAs("uniform.png");


}

Double_t momentum(Double_t m, Double_t m1, Double_t m2) {
    if(m != 0)
        return TMath::Sqrt((m*m - (m1+m2)*(m1+m2))*(m*m - (m1-m2)*(m1-m2))/(4.0*m*m));
    else 
        return 0;
}

// trzeba to przerobic że liczy cale gamma total
Double_t Gammak(Double_t m,ParticleType* mother,ParticleType* kid1, ParticleType* kid2,Double_t branchingRatio) {
  //Double_t czynnik = mother->GetMass() / m;
  Double_t czynnik = 1.0;
  Double_t Ltr = TMath::Abs(mother->GetSpin()-(kid1->GetSpin() + kid2->GetSpin())); 
  Double_t pierwszyNawias = TMath::Power( momentum(m,kid1->GetMass(),kid2->GetMass()) / momentum(mother->GetMass(),kid1->GetMass(),kid2->GetMass()) , 2.0*Ltr+1);
  Double_t beta = TMath::Power(mother->GetMass() - kid1->GetMass() - kid2->GetMass(),2) + TMath::Power(mother->GetGamma(),2)/4.0;
  Double_t drugiNawias = (beta*beta + momentum(mother->GetMass(),kid1->GetMass(),kid2->GetMass()))/(beta*beta + momentum(m,kid1->GetMass(),kid2->GetMass()));
  Double_t gammaK = branchingRatio * mother->GetGamma();

  return czynnik * pierwszyNawias * drugiNawias *gammaK;

}

Double_t BWwithGamma(Double_t *x, Double_t *par) {
    Double_t gamma = 
}


void createBWSpectraDiffGamma() {

    Int_t pid[_N_PARTICLES_] = {12112, 1214, 22112, 32214, 2122, 32212, 2216, 12216};

    ParticleDB *particles = new ParticleDB();
    ReadSHAREParticles(particles);
    ReadSHAREDecays(particles);

    DecayTable* decaysTables[_N_PARTICLES_];
    std::vector<const DecayChannel*> decaysChannels[_N_PARTICLES_];


    ParticleType *particle[_N_PARTICLES_];
    Double_t mMin[_N_PARTICLES_];
    Double_t mMax[_N_PARTICLES_];
    Double_t gamma, mass;
   
    TF1* fun[_N_PARTICLES_]; 
    TCanvas *can[_N_PARTICLES_];
    TH1D *hist[_N_PARTICLES_];
    TH1D *histTest[_N_PARTICLES_];

    ParticleType* part1;
    ParticleType* part2;
    //ParticleType* part3;
    Double_t massSum, massTresh[_N_PARTICLES_];

    for(int i = 0; i < _N_PARTICLES_; i++) {


        TF1* funtmp;
    
    particle[i] = particles->FindByPID(pid[i]);
    decaysTables[i] = particle[i]->GetTable();
    cout<<particle[i]->GetName()<<"\t"<<particle[i]->GetSpin()<<endl;
    //cout<<decaysTables[i]->GetChannelCount()<<endl;
    for(int j = 0; j < decaysTables[i]->GetChannelCount(); j++) {
      decaysChannels[i].push_back(decaysTables[i]->GetDecayChannel(j));
      part1 = particles->GetParticleType(decaysChannels[i][j]->GetParticle1());
      part2 = particles->GetParticleType(decaysChannels[i][j]->GetParticle2());
      //part3 = particles->GetParticleType(decaysChannels[i][j]->GetParticle3());

      //cout<<part1->GetName()<<" + "<<part2->GetName();
      //cout<<"\t"<<part1->GetSpin()<<" + "<<part2->GetSpin()<<endl;
      massSum = part1->GetMass() + part2->GetMass();
      if (j == 0)
        massTresh[i] = massSum;
      else {
        if(massSum < massTresh[i])
          massTresh[i] = massSum;
      }
      //cout<<"\t"<<massSum<<endl;;
      //cout<<decaysChannels[i][j]->GetParticle1()-><<" + "<<decaysChannels[i][j]->GetParticle2()<<decaysChannels[i][j]->GetParticle3()<<endl;
      funtmp = new TF1(Form("gamma%s%i",particle[i]->GetName(),j),Gammak,0.5,2.0,0);

        
    }
    //cout<<"\t"<<massTresh[i]<<endl<<endl;

    gamma = particle[i]->GetGamma();
    //cout<<gamma<<endl;
    mass = particle[i]->GetMass();
    //cout<<mass<<endl;
    mMin[i] = massTresh[i];
    mMax[i] = mass + 3.0*gamma;
    
    fun[i] = new TF1(Form("fun%s",particle[i]->GetName()),breitWigner,mMin[i],mMax[i],3); //mass from 0-3, 3 parameters
    fun[i]->SetParameter(0,1); 
    fun[i]->SetParameter(1,particle[i]->GetGamma());
    fun[i]->SetParameter(2,particle[i]->GetMass());
    fun[i]->SetParameter(0,1./fun[i]->Integral(mMin[i],mMax[i]));
    

    hist[i] = (TH1D*)fun[i]->GetHistogram();
    hist[i]->SetTitle(particle[i]->GetName());
    //cout<<hist[i]->Integral("width")<<endl;
    //cout<<hist[i]->GetNbinsX()<<endl;
    //cout<<mMin[i]<<endl;


    can[i]= new TCanvas(Form("can%s",particle[i]->GetName()),Form("can%s",particle[i]->GetName()),1000,1000);
    can[i]->cd();
    //fun[i]->Draw();
    hist[i]->Draw();
    //can[i]->SaveAs(Form("outputMassSpectra/%sMassTresh.png",particle[i]->GetName()));


}




}


void massSpectra() {
   
    //createBWSpectra();
    //createBWSpectraDiffRange();
    //test();
    createBWSpectraDiffGamma();
  
  
  
   // ParticleDB *particles = new ParticleDB();
    //ReadSHAREParticles(particles);
    //ReadSHAREDecays(particles);
}