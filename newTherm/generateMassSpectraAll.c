#include "moje.h"
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>

#define _N_PARTICLES_ 15


ParticleDB *particlesDB;


std::map<int, ParticleType*> *particlesMap;
//std::map<int, DecayTable*> *decaysTablesMap;

ParticleType *particles[_N_PARTICLES_];
DecayTable* decaysTables[_N_PARTICLES_];

 Int_t pid[_N_PARTICLES_] = {12112, 1214, 22112, 32214, 2122, 32212, 2216, 12216,5218, 2128,32124,22124, 9299,2218,12214};//12214
 //Int_t pid[_N_PARTICLES_] = {12112};
 int colors[]={1, 600, 629, 414, 802, 880, 819, 922,433,618};
 Double_t unstableTresh = 1e-6;

Double_t BWFinalSpectra(Double_t *x, Double_t *par);

Double_t momentum(Double_t m, Double_t m1, Double_t m2) {
    if(m != 0 && (m > (m1 + m2)))
        return TMath::Sqrt(TMath::Abs((m*m - (m1+m2)*(m1+m2))*(m*m - (m1-m2)*(m1-m2))/(4.0*m*m)));
    else 
        return 0;
}
Double_t momentum2D(Double_t *x, Double_t *par) {// 0 - normalizacja, 1 - masa

    if(par[1] != 0 && (par[1] > (x[0] + x[1])))
       return  par[0]*TMath::Sqrt(TMath::Abs((par[1]*par[1] - (x[0]+x[1])*(x[0]+x[1]))*(par[1]*par[1] - (x[0]-x[1])*(x[0]-x[1]))/(4.0*par[1]*par[1])));
    else 
        return 0;


}

Double_t BWConstGamma(Double_t *x, Double_t *par) { //0 - normalization, 1 - Gamma, 2 - M_R
    if((x[0] > (par[2] - 2.0*par[1])) && (x[0] < (par[2] + 12.0*par[1])))
        return par[0]*x[0]*x[0]*par[1]/((par[2]*par[2]-x[0]*x[0])*(par[2]*par[2]-x[0]*x[0]) + x[0]*x[0]*par[1]*par[1]);
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

Double_t GammaK2StableForDrawing(Double_t *x, Double_t *par) {//0-normalizacja, 1 -pid matki, 2 - pid kid1, 3 - pid kid2, 4 - branching ratio
    return GammaK2Stable(x[0],particlesDB->FindByPID((int)par[1]),particlesDB->FindByPID((int)par[2]),particlesDB->FindByPID((int)par[3]),par[4]);


}


Double_t GammaK2Unstable(Double_t m,ParticleType* mother,ParticleType* kid1, ParticleType* kid2,Double_t branchingRatio) {

    
    Double_t mMin1 = kid1->GetMass() - 2.0*kid1->GetGamma();
    Double_t mMax1 = kid1->GetMass() + 12.0*kid1->GetGamma();


    Double_t mMin2 = kid2->GetMass() - 2.0*kid2->GetGamma();
    Double_t mMax2 = kid2->GetMass() + 12.0*kid2->GetGamma();


    TF2 *fun = new TF2(Form("momentum%s+%s",kid1->GetName(),kid2->GetName()),momentum2D,mMin1,mMax1,mMin2,mMax2,2);
    fun->SetParameters(1.0,m);
    fun->SetNpx(100);
    fun->SetNpy(100);
    TH2D* hist = (TH2D*)fun->GetHistogram();


    TF1* bw1 = new TF1(kid1->GetName(),BWConstGamma,mMin1,mMax1,3);
    bw1->SetParameters(1.0,kid1->GetGamma(),kid1->GetMass());
    bw1->SetParameter(0,1.0/bw1->Integral(mMin1,mMax1));

    TF1* bw2 = new TF1(kid2->GetName(),BWConstGamma,mMin2,mMax2,3);
    bw2->SetParameters(1.0,kid2->GetGamma(),kid2->GetMass());
    bw2->SetParameter(0,1.0/bw2->Integral(mMin2,mMax2));

    for(int i = 1; i <= hist->GetNbinsX(); i++) {
        for(int j = 1; j <= hist->GetNbinsY(); j++) {
           hist->SetBinContent(i,j,hist->GetBinContent(i,j)*bw1->Eval(hist->GetXaxis()->GetBinCenter(i))*bw2->Eval(hist->GetYaxis()->GetBinCenter(j)));
        }
    }

    Double_t integral = hist->Integral(1,hist->GetNbinsX(),1,hist->GetNbinsY(),"width");
    
    if(m > (kid1->GetMass() + kid2->GetMass()))
        return integral*GammaK2Stable(m,mother,kid1,kid2,branchingRatio);
    else
        return 0;


}

Double_t GammaK2UnstableForDrawing(Double_t *x, Double_t *par) { //0-normalizacja, 1 -pid matki, 2 - pid kid1, 3 - pid kid2, 4 - branching ratio
    return GammaK2Unstable(x[0],particlesDB->FindByPID((int)par[1]),particlesDB->FindByPID((int)par[2]),particlesDB->FindByPID((int)par[3]),par[4]);

}

Double_t GammaK1Unstable(Double_t m,ParticleType* mother,ParticleType* kid1, ParticleType* kid2,Double_t branchingRatio) { // kid1 -> unstable
    Double_t min = kid1->GetMass() - 2.0*kid1->GetGamma();
    Double_t max = kid1->GetMass() + 12.0*kid1->GetGamma();
    TF1 *fun = new TF1(kid1->GetName(),BWConstGamma,min,max,3);
    //TF1 *fun = new TF1(kid1->GetName(),BWFinalSpectra,min,max,2);
    fun->SetParameters(1.0,kid1->GetGamma(),kid1->GetMass());
    fun->SetNpx(500);
    TH1D* hist = (TH1D*)fun->GetHistogram();
    hist->Scale(1.0/hist->Integral("width"));
    for(int i = 1; i < hist->GetNbinsX(); i++)
        hist->SetBinContent(i,hist->GetBinContent(i)*momentum(m,hist->GetBinCenter(i),kid2->GetMass()));
    hist->Scale(GammaK2Stable(m,mother,kid1,kid2,branchingRatio));
    return hist->Integral("width");
}


Double_t GammaK1UnstableForDrawing(Double_t *x, Double_t *par) { //0-normalizacja, 1 -pid matki, 2 - pid kid1, 3 - pid kid2, 4 - branching ratio
    return GammaK1Unstable(x[0],particlesDB->FindByPID((int)par[1]),particlesDB->FindByPID((int)par[2]),particlesDB->FindByPID((int)par[3]),par[4]);

}
Double_t ThreebodyDecayForDrawing(Double_t *x, Double_t *par) {//0 - y value
    return 0.0*x[0] + par[0];
}



Double_t BWFinalSpectra(Double_t *x, Double_t *par) { // 0 - normalizacja, 1 - pdg czastki
    

    ParticleType* mother = particlesMap->find((Int_t)par[1])->second;
    DecayTable* decayTable = mother->GetTable();
    ParticleType* kid1;
    ParticleType* kid2;
    //ParticleType* kid3;
    Double_t gammaTot = 0;
    for(int i = 0; i <= decayTable->GetChannelCount(); i++) { //pętla po kanałach rozpadów
        kid1 = particlesDB->GetParticleType(decayTable->GetDecayChannel(i)->GetParticle1());
        kid2 = particlesDB->GetParticleType(decayTable->GetDecayChannel(i)->GetParticle2());

        if(!decayTable->GetDecayChannel(i)->Is3Particle()) { //2 ciałowy
            if ( kid1->GetTable()->GetChannelCount() < 0 && kid2->GetTable()->GetChannelCount() < 0) { //2 stabilne
                gammaTot += GammaK2Stable(x[0],mother,kid1,kid2,decayTable->GetDecayChannel(i)->GetBranchingRatio());
            }
            else if (kid1->GetTable()->GetChannelCount() >= 0 && kid2->GetTable()->GetChannelCount() >= 0 ) { // 2 niestabilne

                if(kid1->GetGamma() > unstableTresh && kid2->GetGamma() > unstableTresh) {
                    gammaTot += GammaK2Unstable(x[0],mother,kid1,kid2,decayTable->GetDecayChannel(i)->GetBranchingRatio());
                }
                else if(kid1->GetGamma() < unstableTresh && kid2->GetGamma() < unstableTresh) { //jednak stabilne
                    gammaTot += GammaK2Stable(x[0],mother,kid1,kid2,decayTable->GetDecayChannel(i)->GetBranchingRatio());
                }
                else if (kid1->GetGamma() > unstableTresh) { //1 - niestabilna
                    gammaTot += GammaK1Unstable(x[0],mother,kid1,kid2,decayTable->GetDecayChannel(i)->GetBranchingRatio());
                }
                else if(kid2->GetGamma() > unstableTresh) {//2 - niestabilna
                    gammaTot += GammaK1Unstable(x[0],mother,kid2,kid1,decayTable->GetDecayChannel(i)->GetBranchingRatio());
                }



                //gammaTot += GammaK2Stable(x[0],mother,kid1,kid2,decayTable->GetDecayChannel(i)->GetBranchingRatio());
            }
            else { // 1 niestabilna
                if(kid1->GetTable()->GetChannelCount() >= 0) //zawsze pierwsza podana jest niestabilna
                    gammaTot += GammaK1Unstable(x[0],mother,kid1,kid2,decayTable->GetDecayChannel(i)->GetBranchingRatio());
                else if(kid2->GetTable()->GetChannelCount() >= 0)
                    gammaTot += GammaK1Unstable(x[0],mother,kid2,kid1,decayTable->GetDecayChannel(i)->GetBranchingRatio());
            
            }
        }
        else { //3 ciałowy
            //kid3 = particlesDB->GetParticleType(decaysTables[index]->GetDecayChannel(i)->GetParticle3());
            gammaTot += decayTable->GetDecayChannel(i)->GetBranchingRatio() * mother->GetGamma();

        }
        
    }
    Double_t poleMass = mother->GetMass();
    if (((poleMass*poleMass-x[0]*x[0])*(poleMass*poleMass-x[0]*x[0]) + x[0]*x[0]*gammaTot*gammaTot) == 0 )
        return 0;
    else 
        return par[0]*x[0]*x[0]*gammaTot/((poleMass*poleMass-x[0]*x[0])*(poleMass*poleMass-x[0]*x[0]) + x[0]*x[0]*gammaTot*gammaTot);

}

void massSpectraFinal() {

    ParticleType* particle;

    std::vector<Double_t> mMax;
    std::vector<Double_t> mTresh;
    std::vector<TH1D*> hist;
    std::vector<TF1*> fun;


    //CHOOSING PARTICLES
    particlesMap = new std::map<int, ParticleType*>;

    for(int i = 0; i < particlesDB->GetParticleTypeCount(); i++) {
        particle = particlesDB->GetParticleType(i);
        
        if(particle->GetTable()->GetChannelCount()>=0 && particle->GetGamma() > unstableTresh && !(TMath::Abs(particle->GetPDGCode()) == 2224 || TMath::Abs(particle->GetPDGCode()) == 2214 || TMath::Abs(particle->GetPDGCode()) == 2114 || TMath::Abs(particle->GetPDGCode()) == 1114)) {
            if(particlesMap->find(TMath::Abs(particle->GetPDGCode())) == particlesMap->end())
                (*particlesMap)[TMath::Abs(particle->GetPDGCode())] = particlesDB->GetParticleType(i);
        }
    }

    //cout<<particlesMap->size()<<endl;


    //CALCULATING M TRESH AND BREIT WIGNER
    Double_t massSum;
    Double_t massTresh;
    ParticleType* tmp1;
    ParticleType* tmp2;
    ParticleType* tmp3;
    TFile *fileOut = new TFile("outputAllBW/BW.root","RECREATE");

    Int_t licznik = 1;

    for (auto const& [key, part] : (*particlesMap)) {
        //ParticleType* part = (*particlesMap)[1214];
        cout<<licznik<<"/"<<particlesMap->size()<<":\t"<<part->GetName()<<endl;
        for(int i = 0; i <= part->GetTable()->GetChannelCount(); i++) {

            tmp1 = (ParticleType*)particlesDB->GetParticleType(part->GetTable()->GetDecayChannel(i)->GetParticle1());
            tmp2 = (ParticleType*)particlesDB->GetParticleType(part->GetTable()->GetDecayChannel(i)->GetParticle2());
          //cout<<part1->GetName()<<" + "<<part2->GetName();
            if(!part->GetTable()->GetDecayChannel(i)->Is3Particle()) 
                massSum = tmp1->GetMass() + tmp2->GetMass();
            else {
                tmp3 = (ParticleType*)particlesDB->GetParticleType(part->GetTable()->GetDecayChannel(i)->GetParticle3());
                massSum = tmp1->GetMass() + tmp2->GetMass() + tmp3->GetMass();
            }
            
            if (i == 0)
                massTresh = massSum;
            else {
                if(massSum < massTresh)
                    massTresh = massSum;
            }
        }
        mTresh.push_back(massTresh);
        mMax.push_back(part->GetMass() + 12.0 * part->GetGamma());
        //cout<<mTresh.back()<<"\t"<<mMax.back()<<endl;

        fun.push_back( new TF1(Form("%i",key),BWFinalSpectra,mTresh.back(),mMax.back(),2));
        fun.back()->SetParameter(0,1);
        fun.back()->SetParameter(1,key);
        fun.back()->SetNpx(500);

        
        hist.push_back((TH1D*)fun.back()->GetHistogram()->Clone(Form("%i",key)));
        hist.back()->Scale(1.0/hist.back()->Integral("width"));
        hist.back()->SetOption("L");
        fileOut->cd();
        //fun.back()->Write();
        hist.back()->Write();
        licznik++;




    }

    fileOut->Close();
    fileOut->Save();


}

void test2Unstable() {
    cout<<particlesDB->GetParticleTypeCount()<<endl;
    DecayTable *table;
    ParticleType* kid1;
    ParticleType* kid2;
    int nTot = 0;


    for(int i = 0; i < particlesDB->GetParticleTypeCount(); i++) {
        table = particlesDB->GetParticleType(i)->GetTable();
        for(int j = 0; j <= table->GetChannelCount(); j++) {
            if(!table->GetDecayChannel(j)->Is3Particle()) {
                kid1 = particlesDB->GetParticleType(table->GetDecayChannel(j)->GetParticle1());
                kid2 = particlesDB->GetParticleType(table->GetDecayChannel(j)->GetParticle2());
                if (kid1->GetTable()->GetChannelCount() >= 0 && kid2->GetTable()->GetChannelCount() >= 0 ) {
                    if(kid1->GetGamma() > 1e-5 && kid2->GetGamma() > 1e-5) {
                        cout<<"2 niestabilne"<<endl;
                        cout<<particlesDB->GetParticleType(i)->GetName()<<" = "<<kid1->GetName()<<" + "<<kid2->GetName()<<endl;
                        cout<<"********"<<endl;
                        nTot++;
                    }
                }
            }
            else {
                cout<<particlesDB->GetParticleType(i)->GetName()<<" = "<<kid1->GetName()<<" + "<<kid2->GetName() <<" + "<<particlesDB->GetParticleType(table->GetDecayChannel(j)->GetParticle3())->GetName()<<endl;
            }
        }



    }

    cout<<nTot<<endl;
}

void drawingGamma() {
    ParticleType* particle;

    std::vector<Double_t> mMax;
    std::vector<Double_t> mTresh;
    //std::vector<TH1D*> hist;
    // std::vector<TF1*> fun;


    //CHOOSING PARTICLES
    particlesMap = new std::map<int, ParticleType*>;

    for(int i = 0; i < particlesDB->GetParticleTypeCount(); i++) {
        particle = particlesDB->GetParticleType(i);
        
        if(particle->GetTable()->GetChannelCount()>=0 && particle->GetGamma() > unstableTresh && !(TMath::Abs(particle->GetPDGCode()) == 2224 || TMath::Abs(particle->GetPDGCode()) == 2214 || TMath::Abs(particle->GetPDGCode()) == 2114 || TMath::Abs(particle->GetPDGCode()) == 1114)) {
            if(particlesMap->find(TMath::Abs(particle->GetPDGCode())) == particlesMap->end())
                (*particlesMap)[TMath::Abs(particle->GetPDGCode())] = particlesDB->GetParticleType(i);
        }
    }
    //cout<<particlesMap->size()<<endl;

    std::vector<TF1*> gammaK;
    ParticleType* tmp1;
    ParticleType* tmp2;
    ParticleType* tmp3;
    std::vector<TCanvas*> can;
   // TCanvas *can[_N_PARTICLES_];

    TH1D *gammaTot;
    Double_t massSum;
    Double_t massTresh;
    Int_t licznik = 1;

    TFile *fileOut = new TFile("outputAllBW/diffGammaFinal.root","RECREATE");

    for (auto const& [key, part] : (*particlesMap)) {
         cout<<licznik<<"/"<<particlesMap->size()<<":\t"<<part->GetName()<<endl;
        for(int i = 0; i <= part->GetTable()->GetChannelCount(); i++) {

            tmp1 = (ParticleType*)particlesDB->GetParticleType(part->GetTable()->GetDecayChannel(i)->GetParticle1());
            tmp2 = (ParticleType*)particlesDB->GetParticleType(part->GetTable()->GetDecayChannel(i)->GetParticle2());
          //cout<<part1->GetName()<<" + "<<part2->GetName();
            if(!part->GetTable()->GetDecayChannel(i)->Is3Particle()) 
                massSum = tmp1->GetMass() + tmp2->GetMass();
            else {
                tmp3 = (ParticleType*)particlesDB->GetParticleType(part->GetTable()->GetDecayChannel(i)->GetParticle3());
                massSum = tmp1->GetMass() + tmp2->GetMass() + tmp3->GetMass();
            }
            
            if (i == 0)
                massTresh = massSum;
            else {
                if(massSum < massTresh)
                    massTresh = massSum;
            }
        }
        mTresh.push_back(massTresh);
        mMax.push_back(part->GetMass() + 12.0 * part->GetGamma());
        can.push_back(new TCanvas(Form("%ican%s",licznik,part->GetName()),Form("%ican%s",licznik,part->GetName()),1000,1000));

        for(int i = 0; i <= part->GetTable()->GetChannelCount(); i++) {
            tmp1 = (ParticleType*)particlesDB->GetParticleType(part->GetTable()->GetDecayChannel(i)->GetParticle1());
            tmp2 = (ParticleType*)particlesDB->GetParticleType(part->GetTable()->GetDecayChannel(i)->GetParticle2());
            if(!part->GetTable()->GetDecayChannel(i)->Is3Particle()) { //2 ciałowy
                if ( tmp1->GetTable()->GetChannelCount() < 0 && tmp2->GetTable()->GetChannelCount() < 0) { //2 stabilne
                
                    gammaK.push_back(new TF1(Form("%s+%s",tmp1->GetName(),tmp2->GetName()),GammaK2StableForDrawing,mTresh.back(),mMax.back(),5));
                    gammaK.back()->SetParameters(1.0,key,tmp1->GetPDGCode(),tmp2->GetPDGCode(),part->GetTable()->GetDecayChannel(i)->GetBranchingRatio());

                }
                else if (tmp1->GetTable()->GetChannelCount() >= 0 && tmp2->GetTable()->GetChannelCount() >= 0) { // 2 niestabilne

                    if(tmp1->GetGamma() > unstableTresh && tmp2->GetGamma() > unstableTresh) {
                        //cout<<"2 na SERIO niestabilne: "<<tmp1->GetName()<<" +"<<tmp2->GetName()<<endl;
                        gammaK.push_back(new TF1(Form("%s+%s",tmp1->GetName(),tmp2->GetName()),GammaK2UnstableForDrawing,mTresh.back(),mMax.back(),5));
                        gammaK.back()->SetParameters(1.0,key,tmp1->GetPDGCode(),tmp2->GetPDGCode(),part->GetTable()->GetDecayChannel(i)->GetBranchingRatio());

                    }
                    else if(tmp1->GetGamma() < unstableTresh && tmp2->GetGamma() < unstableTresh) { //jednak stabilne
                        //cout<<"2 nie na serio niestabilne: "<<tmp1->GetName()<<" +"<<tmp2->GetName()<<endl;
                        gammaK.push_back(new TF1(Form("%s+%s",tmp1->GetName(),tmp2->GetName()),GammaK2StableForDrawing,mTresh.back(),mMax.back(),5));
                        gammaK.back()->SetParameters(1.0,key,tmp1->GetPDGCode(),tmp2->GetPDGCode(),part->GetTable()->GetDecayChannel(i)->GetBranchingRatio());
                    }
                    else if (tmp1->GetGamma() > unstableTresh) { //1 - niestabilna
                        //cout<<"2 niestabilne, ale 1. serio niestabilna: "<<tmp1->GetName()<<" +"<<tmp2->GetName()<<endl;
                        gammaK.push_back(new TF1(Form("%s+%s",tmp1->GetName(),tmp2->GetName()),GammaK1UnstableForDrawing,mTresh.back(),mMax.back(),5));
                        gammaK.back()->SetParameters(1.0,key,tmp1->GetPDGCode(),tmp2->GetPDGCode(),part->GetTable()->GetDecayChannel(i)->GetBranchingRatio());

                    }
                    else if(tmp2->GetGamma() > unstableTresh) {//2 - niestabilna
                        //cout<<"2 niestabilne, ale 2. serio niestabilna: "<<tmp1->GetName()<<" +"<<tmp2->GetName()<<endl;
                        gammaK.push_back(new TF1(Form("%s+%s",tmp2->GetName(),tmp1->GetName()),GammaK1UnstableForDrawing,mTresh.back(),mMax.back(),5));
                        gammaK.back()->SetParameters(1.0,key,tmp2->GetPDGCode(),tmp1->GetPDGCode(),part->GetTable()->GetDecayChannel(i)->GetBranchingRatio());
                    }
                }
                else { // 1 niestabilna
                    if(tmp1->GetTable()->GetChannelCount() >= 0){ //zawsze pierwsza podana jest niestabilna
                        gammaK.push_back(new TF1(Form("%s+%s",tmp1->GetName(),tmp2->GetName()),GammaK1UnstableForDrawing,mTresh.back(),mMax.back(),5));
                        gammaK.back()->SetParameters(1.0,key,tmp1->GetPDGCode(),tmp2->GetPDGCode(),part->GetTable()->GetDecayChannel(i)->GetBranchingRatio());
                    }
                    else if(tmp2->GetTable()->GetChannelCount() >= 0) {
                        gammaK.push_back(new TF1(Form("%s+%s",tmp2->GetName(),tmp1->GetName()),GammaK1UnstableForDrawing,mTresh.back(),mMax.back(),5));
                        gammaK.back()->SetParameters(1.0,key,tmp2->GetPDGCode(),tmp1->GetPDGCode(),part->GetTable()->GetDecayChannel(i)->GetBranchingRatio());
                    }
                }
            
            }  
            else {
                tmp3 = (ParticleType*)particlesDB->GetParticleType(part->GetTable()->GetDecayChannel(i)->GetParticle3());
               // cout<<"3 ciałowy: "<<kid1->GetName()<<" + "<<kid2->GetName()<<" + "<<kid3->GetName()<<endl;
                gammaK.push_back(new TF1(Form("%s+%s+%s",tmp1->GetName(),tmp2->GetName(),tmp3->GetName()),ThreebodyDecayForDrawing,mTresh.back(),mMax.back(),1));
                gammaK.back()->SetParameter(0,part->GetTable()->GetDecayChannel(i)->GetBranchingRatio() * part->GetGamma());

            }

            gammaK.back()->SetLineColor(colors[gammaK.size()-1]);
            gammaK.back()->SetLineStyle(2);
            gammaK.back()->SetMinimum(0);

        }

        gammaTot = new TH1D(part->GetName(),part->GetName(),100,mTresh.back(),mMax.back());
        for (auto gammaChannels : gammaK){
            // gammaChannels->SetNpx(200);
                gammaTot->Add(gammaChannels->GetHistogram());
        }
        gammaTot->SetMinimum(0);
        gammaTot->SetLineColor(kRed);
        gammaTot->SetLineWidth(2);

        can.back()->cd();
        gammaTot->Draw("L");
        for (auto gammaChannels : gammaK)
            gammaChannels->Draw("same");
        fileOut->cd();
        can.back()->Write();

        gammaK.clear();


        licznik++;
    }

    fileOut->Close();
    fileOut->Save();

    cout<<particlesMap->size()<<endl;
    

}

void testGammaTresh(){  
    TCanvas* can = new TCanvas("can","can",1000,1000);
    ParticleType *part = particlesDB->FindByPID(3212);
    Double_t mass = part->GetMass();
    Double_t gamma =  part->GetGamma();
    TF1 *fun = new TF1(part->GetName(),BWConstGamma,mass - 2.0 * gamma,mass + 12.0 * gamma,3);
    fun->SetParameters(1.0,gamma,mass);
    can->cd();
    fun->Draw();
    can->SaveAs("./outputAllBW/test.png");


}

void generateMassSpectraAll() {
        
    particlesDB = new ParticleDB();
    ReadSHAREParticles(particlesDB);
    ReadSHAREDecays(particlesDB);

    gStyle->SetOptStat(0);
   massSpectraFinal();


    drawingGamma();
   //GammaK2Unstable(2.5,particlesDB->FindByPID(40225),particlesDB->FindByPID(333), particlesDB->FindByPID(333),1.0);
   //test2Unstable();
    //testGammaTresh();




}