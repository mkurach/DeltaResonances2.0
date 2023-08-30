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
    if((x[0] > par[2] - 2.0*par[1]) && (x[0] < par[2] + 12.0*par[1]))
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

    
   /*
    
    Double_t mMin1 = kid1->GetMass() - 2.0*kid1->GetGamma();
    Double_t mMax1 = kid1->GetMass() + 12.0*kid1->GetGamma();


    Double_t mMin2 = kid2->GetMass() - 2.0*kid2->GetGamma();
    Double_t mMax2 = kid2->GetMass() + 12.0*kid2->GetGamma();

    TFile* file = new TFile("outputFinalBW/test.root","RECREATE");


    //TF2 *fun = new TF2(Form("momentum%s+%s",kid1->GetName(),kid2->GetName()),momentum2D,mMin1,mMax1,mMin2,mMax2,2);
    TF2 *fun = new TF2(Form("momentum%s+%s",kid1->GetName(),kid2->GetName()),momentum2D,0.497,5,0,5,2);
    fun->SetParameters(1.0,m);
    //fun->SetNpx(5);
    //fun->SetNpy(5);
    //TH2D* hist = (TH2D*)fun->GetHistogram();
    //hist->Scale(1.0/hist->Integral("width"));
    //cout<<hist->GetNbinsX()<<endl;
    //cout<<hist->GetNbinsY()<<endl;
    Double_t x,y;
    for(int i = 0; i <= hist->GetNbinsX(); i++) {
        for(int j = 0; j <= hist->GetNbinsY(); j++) {
            
           // cout<<"GetBin(i,j): "<<hist->GetBin(i,j)<<endl;
           cout<<"i: "<<i<<"\tj: "<<j<<endl;
           //cout<<"GetBin(i,j): "<<hist->GetBin(i,j)<<endl;
           //cout<<"FindBin x axis: "<<hist->GetXaxis()->FindBin(i)<<endl;
           //cout<<"FindBin y axis: "<<hist->GetYaxis()->FindBin(j)<<endl;
           //x = hist->GetXaxis()->GetBinCenter(hist->GetXaxis()->FindBin(i));
           //y = hist->GetYaxis()->GetBinCenter(hist->GetYaxis()->FindBin(j));
            x = hist->GetXaxis()->GetBinCenter(i);
            y = hist->GetYaxis()->GetBinCenter(j);

           cout<<"x: "<<x<<"\ty: "<<y<<endl;
           //hist->SetBinContent(i,j,hist->GetBinContent(i,j)*momentum(m,hist->GetBinCenter(hist->GetBin(i,1)),hist->GetBinCenter(hist->GetBin(1,j))));

        }
    }

    TCanvas *can = new TCanvas("can","can",1000,1000);
    can->cd();
    //hist->Draw("colz");
    fun->Draw();
    can->SaveAs("outputFinalBW/test.png");
    file->cd();
    fun->Write();
    file->Save();
    file->Close();
    
    cout<<mMin2<<endl;
    cout<<mMax2<<endl;
    cout<<mMin1<<endl;
    cout<<mMax1<<endl;
    //if(m > (kid1->GetMass() + kid2->GetMass()))
       // return hist->Integral("width")*GammaK2Stable(m,mother,kid1,kid2,branchingRatio);
   // else   

   */ 
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
                gammaTot += GammaK2Stable(x[0],mother,kid1,kid2,decayTable->GetDecayChannel(i)->GetBranchingRatio());
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
        
        if(particle->GetTable()->GetChannelCount()>=0 && particle->GetGamma() > 1e-5 && !(TMath::Abs(particle->GetPDGCode()) == 2224 || TMath::Abs(particle->GetPDGCode()) == 2214 || TMath::Abs(particle->GetPDGCode()) == 2114 || TMath::Abs(particle->GetPDGCode()) == 1114)) {
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

    for (auto const& [key, part] : (*particlesMap)) {
        cout<<key<<"\t"<<part->GetName()<<endl;
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

        fun.push_back( new TF1(part->GetName(),BWFinalSpectra,mTresh.back(),mMax.back(),2));
        fun.back()->SetParameter(0,1);
        fun.back()->SetParameter(1,key);
        fun.back()->SetNpx(500);

        
        hist.push_back((TH1D*)fun.back()->GetHistogram());
        hist.back()->Scale(1.0/hist.back()->Integral("width"));
        fileOut->cd();
        fun.back()->Write();




    }

    fileOut->Close();
    fileOut->Save();
    //for(auto m : mTresh)
        //cout<<m<<endl;
    //cout<<particlesMap->find(1214)->second->GetName()<<endl;

    //CALCULATING BREIT WIGNER   
    //TF1* fun[_N_PARTICLES_];
    //TH1D* hist[_N_PARTICLES_];





    /*for(int i = 0; i < particlesDB->GetParticleTypeCount(); i++) {
        particle = particlesDB->GetParticleType(i);
        decaysTable = particlesDB->GetParticleType(i)->GetTable();
        
        if(decaysTable->GetChannelCount()>=0 && particle->GetGamma() > 1e-5 && !(TMath::Abs(particle->GetPDGCode()) == 2224 || TMath::Abs(particle->GetPDGCode()) == 2214 || TMath::Abs(particle->GetPDGCode()) == 2114 || TMath::Abs(particle->GetPDGCode()) == 1114)) {
            cout<<particle->GetName()<<endl;
            //cout<<"Liczba kanałow "<<decaysTable->GetChannelCount()+1<<endl;
            nTot++;
            for(int j = 0; j <= decaysTable->GetChannelCount(); j++) {
                if(!decaysTable->GetDecayChannel(j)->Is3Particle()) {
                    tmp1 = (ParticleType*)particlesDB->GetParticleType(decaysTable->GetDecayChannel(j)->GetParticle1());
                    tmp2 = (ParticleType*)particlesDB->GetParticleType(decaysTable->GetDecayChannel(j)->GetParticle2());
                   // cout<<Form("%s+%s",tmp1->GetName(),tmp2->GetName());
                // cout<<tmp1->GetTable()->GetChannelCount()<<endl;
                    //cout<<tmp2->GetTable()->GetChannelCount()<<endl;
                   // cout<<"\t"<<decaysTable->GetDecayChannel(j)->GetBranchingRatio()<<endl;
                }
                else {
                    //out<<"3ciałowy!"<<endl;
                    tmp3 = (ParticleType*)particlesDB->GetParticleType(decaysTable->GetDecayChannel(j)->GetParticle3());
                    //cout<<Form("%s+%s+%s",tmp1->GetName(),tmp2->GetName(),tmp3->GetName());
                    //cout<<"branching ratio: "<<decaysTable->GetDecayChannel(j)->GetBranchingRatio()<<endl;
                }
            }
        }*/

        /*for(int j = 0; j <= decaysTables[i]->GetChannelCount(); j++) {

            tmp1 = (ParticleType*)particlesDB->GetParticleType(decaysTables[i]->GetDecayChannel(j)->GetParticle1());
            tmp2 = (ParticleType*)particlesDB->GetParticleType(decaysTables[i]->GetDecayChannel(j)->GetParticle2());
          //cout<<part1->GetName()<<" + "<<part2->GetName();
            if(!decaysTables[i]->GetDecayChannel(j)->Is3Particle()) 
                massSum = tmp1->GetMass() + tmp2->GetMass();
            else {
                tmp3 = (ParticleType*)particlesDB->GetParticleType(decaysTables[i]->GetDecayChannel(j)->GetParticle3());
                massSum = tmp1->GetMass() + tmp2->GetMass() + tmp3->GetMass();
            }
                
            
                if (j == 0)
                    massTresh[i] = massSum;
                else {
                    if(massSum < massTresh[i])
                        massTresh[i] = massSum;
                }
          }
          //cout<<"\t"<<massTresh[i]<<endl;;
          //cout<<decaysChannels[i][j]->GetParticle1()-><<" + "<<decaysChannels[i][j]->GetParticle2()<<decaysChannels[i][j]->GetParticle3()<<endl;

        
        mMin[i] = particles[i]->GetMass() - 2.0*particles[i]->GetGamma();
        //mMin[i] = massTresh[i];
        mMax[i] = particles[i]->GetMass() + 12.0*particles[i]->GetGamma();
        fun[i] = new TF1(particles[i]->GetName(),BWFinalSpectra,mMin[i],mMax[i],2);
        fun[i]->SetParameter(0,1);
        fun[i]->SetParameter(1,i);
        fun[i]->SetNpx(500);

        can[i]= new TCanvas(Form("can%s",particles[i]->GetName()),Form("can%s",particles[i]->GetName()),1000,1000);
        can[i]->cd();


        hist[i] = (TH1D*)fun[i]->GetHistogram();
        hist[i]->Scale(1.0/hist[i]->Integral("width"));
                hist[i]->Draw("hist");
        can[i]->SaveAs(Form("outputFinalBW/%sFinalBW.png",particles[i]->GetName()));




        cout<<endl<<"******************"<<endl<<endl;
    }*/

    //PRINT HISTOGRAMS
   /*ofstream fileTxt;
    fileTxt.open("histogramsFinal-2_12.txt");


    fileTxt<<"TH1D *getBreitWigner(int pdg) {"<<endl;
    fileTxt<<"\tstatic std::map<int, std::vector<Double_t>> histograms;"<<endl;
    fileTxt<<"\tstatic std::map<int, Double_t> mMin;"<<endl;
    fileTxt<<"\tstatic std::map<int, Double_t> mMax;"<<endl;
    for(int i=0; i < _N_PARTICLES_; i++) {
        
        
        fileTxt<<"\thistograms["<<pid[i]<<"] = {";
        for(int j = 1; j <= hist[i]->GetNbinsX(); j++) {
            if (j != hist[i]->GetNbinsX())
                fileTxt<<hist[i]->GetBinContent(j)<<",";
            else    
                fileTxt<<hist[i]->GetBinContent(j)<<"};"<<endl;
        }
        fileTxt<<"\tmMin["<<pid[i]<<"] = "<<mMin[i]<<";"<<endl;
        fileTxt<<"\tmMax["<<pid[i]<<"] = "<<mMax[i]<<";"<<endl;

    }
    fileTxt<<"\tTH1D *hist = new TH1D(Form(\"hist%i\",pdg),Form(\"hist%i\",pdg),500,mMin[pdg],mMax[pdg]);"<<endl;
    fileTxt<<"\tfor(int i = 1; i <= hist->GetNbinsX(); i++)"<<endl;
    fileTxt<<"\t\thist->SetBinContent(i,histograms[pdg][i]);"<<endl;
    fileTxt<<endl<<"\treturn hist;"<<endl<<"}"<<endl<<endl;
    fileTxt.close();*/

        //PRINT FOR MATHEMATICA

    /*hist[0]->Rebin(5);
    hist[0]->Scale(1.0/5);
    ofstream fileMath;
    fileMath.open("Ns1440zerFinal.dat");
    for(int i = 1; i <= hist[0]->GetNbinsX(); i++ ) {
      fileMath<<hist[0]->GetBinCenter(i)<<"\t"<<hist[0]->GetBinContent(i)<<endl;
    }
  fileMath.close();*/

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
        
        if(particle->GetTable()->GetChannelCount()>=0 && particle->GetGamma() > 1e-5 && !(TMath::Abs(particle->GetPDGCode()) == 2224 || TMath::Abs(particle->GetPDGCode()) == 2214 || TMath::Abs(particle->GetPDGCode()) == 2114 || TMath::Abs(particle->GetPDGCode()) == 1114)) {
            if(particlesMap->find(TMath::Abs(particle->GetPDGCode())) == particlesMap->end())
                (*particlesMap)[TMath::Abs(particle->GetPDGCode())] = particlesDB->GetParticleType(i);
        }
    }

    std::vector<TF1*> gammaK[_N_PARTICLES_];
    ParticleType* kid1;
    ParticleType* kid2;
    ParticleType* kid3;
    TCanvas *can[_N_PARTICLES_];

    TH1D *gammaTot[_N_PARTICLES_];

    TFile *fileOut = new TFile("outputFinalBW/diffGammaFinal.root","RECREATE");

    for(int i = 0; i < _N_PARTICLES_; i++) {
        particles[i] = particlesDB->FindByPID(pid[i]);
        decaysTables[i] = particles[i]->GetTable();
        cout<<particles[i]->GetName()<<endl;
        mMin[i] = particles[i]->GetMass() - 2.0*particles[i]->GetGamma();
        mMax[i] = particles[i]->GetMass() + 12.0*particles[i]->GetGamma();
        can[i]= new TCanvas(Form("can%s",particles[i]->GetName()),Form("can%s",particles[i]->GetName()),1000,1000);
        for(int j = 0; j <= decaysTables[i]->GetChannelCount(); j++) {
            kid1 = (ParticleType*)particlesDB->GetParticleType(decaysTables[i]->GetDecayChannel(j)->GetParticle1());
            kid2 = (ParticleType*)particlesDB->GetParticleType(decaysTables[i]->GetDecayChannel(j)->GetParticle2());
            if(!decaysTables[i]->GetDecayChannel(j)->Is3Particle()) {
                if ( kid1->GetTable()->GetChannelCount() < 0 && kid2->GetTable()->GetChannelCount() < 0) { //2 stabilne
                    cout<<"2 stabilne: "<<kid1->GetName()<<" + "<<kid2->GetName()<<endl;

                    gammaK[i].push_back(new TF1(Form("%s+%s",kid1->GetName(),kid2->GetName()),GammaK2StableForDrawing,mMin[i],mMax[i],5));
                    gammaK[i][j]->SetParameters(1.0,pid[i],kid1->GetPDGCode(),kid2->GetPDGCode(),decaysTables[i]->GetDecayChannel(j)->GetBranchingRatio());
                }
                else if (kid1->GetTable()->GetChannelCount() >= 0 && kid2->GetTable()->GetChannelCount() >= 0 ) { // 2 niestabilne
                    cout<<"2 NIESTABILNE: "<<kid1->GetName()<<" + "<<kid2->GetName()<<endl;
                    gammaK[i].push_back(new TF1(Form("%s+%s",kid1->GetName(),kid2->GetName()),GammaK2StableForDrawing,mMin[i],mMax[i],5));
                    gammaK[i][j]->SetParameters(1.0,pid[i],kid1->GetPDGCode(),kid2->GetPDGCode(),decaysTables[i]->GetDecayChannel(j)->GetBranchingRatio());
                }
                else { // 1 niestabilna
                    if(kid1->GetTable()->GetChannelCount() >= 0) {//zawsze pierwsza podana jest niestabilna
                        cout<<"1 niestabilna: "<<kid1->GetName()<<" + "<<kid2->GetName()<<endl;
                        gammaK[i].push_back(new TF1(Form("%s+%s",kid1->GetName(),kid2->GetName()),GammaK1UnstableForDrawing,mMin[i],mMax[i],5));
                        gammaK[i][j]->SetParameters(1.0,pid[i],kid1->GetPDGCode(),kid2->GetPDGCode(),decaysTables[i]->GetDecayChannel(j)->GetBranchingRatio());
                    }
                    else if(kid2->GetTable()->GetChannelCount() >= 0) {
                        cout<<"1 niestabilna: "<<kid2->GetName()<<" + "<<kid1->GetName()<<endl;
                        gammaK[i].push_back(new TF1(Form("%s+%s",kid1->GetName(),kid2->GetName()),GammaK1UnstableForDrawing,mMin[i],mMax[i],5));
                        gammaK[i][j]->SetParameters(1.0,pid[i],kid2->GetPDGCode(),kid1->GetPDGCode(),decaysTables[i]->GetDecayChannel(j)->GetBranchingRatio());
                    }
                }
            }
            else {
                kid3 = (ParticleType*)particlesDB->GetParticleType(decaysTables[i]->GetDecayChannel(j)->GetParticle3());
                cout<<"3 ciałowy: "<<kid1->GetName()<<" + "<<kid2->GetName()<<" + "<<kid3->GetName()<<endl;
                gammaK[i].push_back(new TF1(Form("%s+%s+%s",kid1->GetName(),kid2->GetName(),kid3->GetName()),ThreebodyDecayForDrawing,mMin[i],mMax[i],1));
                gammaK[i][j]->SetParameter(0,decaysTables[i]->GetDecayChannel(j)->GetBranchingRatio() * particles[i]->GetGamma());
            }

            gammaK[i][j]->SetLineColor(colors[j]);
            gammaK[i][j]->SetLineStyle(2);
            gammaK[i][j]->SetMinimum(0);
            //gammaK[i][j]->SetMaximum(1);
            /*can[i]->cd();
            if(j==0)    
                gammaK[i][j]->Draw();
            else
                gammaK[i][j]->Draw("same");*/



        }

        gammaTot[i] = new TH1D(particles[i]->GetName(),particles[i]->GetName(),100,mMin[i],mMax[i]);
        for (auto gammaChannels : gammaK[i]){
           // gammaChannels->SetNpx(200);
            gammaTot[i]->Add(gammaChannels->GetHistogram());
        }
        gammaTot[i]->SetMinimum(0);
        gammaTot[i]->SetLineColor(kRed);
        gammaTot[i]->SetLineWidth(2);

        can[i]->cd();
        gammaTot[i]->Draw("L");
        for (auto gammaChannels : gammaK[i])
            gammaChannels->Draw("same");

        fileOut->cd();
        can[i]->Write();

        
    }
    fileOut->Close();
    fileOut->Save();


}

void generateMassSpectraAll() {
        
    particlesDB = new ParticleDB();
    ReadSHAREParticles(particlesDB);
    ReadSHAREDecays(particlesDB);

    gStyle->SetOptStat(0);
   massSpectraFinal();


    //drawingGamma();
   //GammaK2Unstable(1.8,particlesDB->FindByPID(8117),particlesDB->FindByPID(3122), particlesDB->FindByPID(311),0.08);
   //test2Unstable();
    




}