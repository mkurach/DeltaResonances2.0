#include <iostream>
#include <TFile.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TLine.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TMath.h>
#include <TStyle.h>
#include <TObjArray.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TString.h>
#include <TObjString.h>
#include <TSystem.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include "TSystemDirectory.h"
#include "TList.h"
#include <TF1.h>
#include <vector>
#include "drawStyle.C"

#define _N_PARTICLES_ 3
#define _N_CASES_ 3
#define _N_FIGURES_ 2

Double_t fitFunctionMtMidY(Double_t *x, Double_t *par) { // 0 - Teff, 1 - normalization
    if(x[0] >= 1 && x[0] <= 3)
        return par[1]*x[0]*x[0]*TMath::Exp(-x[0]/par[0]);
    else    
        return 0;
}
 
Double_t fitFunctionMtFullY(Double_t *x, Double_t *par) { // 0 - Teff, 1 - normalization
    if(x[0] >= 1 && x[0] <= 3)
        return par[1]*TMath::Power(x[0],1.5)*TMath::Exp(-x[0]/par[0]);
    else    
        return 0;
}

Double_t fitFunctionMt(Double_t *x, Double_t *par) {
    if(x[0] >= 0 && x[0] <= 3)
        return par[1]*TMath::Exp(-x[0]/par[0]);
    else    
        return 0;
}


void createHistograms(int i) {

    TH1D* histMt[_N_PARTICLES_][_N_FIGURES_];
    TString histNames[_N_PARTICLES_]= {"MtHistProton","MtHistKaons","MtHistPions"};
    TString particlesTitles[_N_PARTICLES_] = {"proton","kaons","pions"};
    TString filesNames[_N_CASES_] = {"outputPart/outCaseAPartMt.root","outputPart/outCaseBPartMt.root","outputPart/outCaseCPartMt.root"};
    TString casesNames[_N_CASES_] = {"CaseA","CaseB","CaseC"};
    TString figuresNames[_N_FIGURES_] = {"MidY","FullY"};

    TFile* file;

    //for(int i = 0; i < _N_CASES_; i++) {
    //for(int i = 0; i < 1; i++){
        file = new TFile(filesNames[i].Data());
        for(int j = 0; j < _N_PARTICLES_; j++) {
            for(int k= 0; k < _N_FIGURES_; k++) {
            
                histMt[j][k] = (TH1D*)file->Get(Form("%s%s",histNames[j].Data(),figuresNames[k].Data()));
                //histMt[j][k]->Sumw2();
                histMt[j][k] = (TH1D*)histMt[j][k]->Clone(Form("%s%s%s",histNames[j].Data(),casesNames[i].Data(),figuresNames[k].Data()));
                //can[i][j]->cd();
                //histMt[i][j]->Draw();
            }
        }
    //} 


    //FUNCTIONS

    Int_t nparams = 2;
    //TF1* fun[_N_PARTICLES_][_N_FIGURES_];
    Double_t teffParameter[_N_PARTICLES_] = {0.129759,0.101,0.076};
    Double_t normParameter[_N_PARTICLES_] = {2.2e6,1e3,1.4e5};
    //Double_t normParameter[_N_PARTICLES_] = {4.76e10,3.7e6,8e7};

    /*for(int j = 0; j < _N_PARTICLES_; j ++) {
        fun[j][0] = new TF1(Form("fun%d0",j),fitFunctionMtMidY,1,3,nparams);
        fun[j][1] = new TF1(Form("fun%d1",j),fitFunctionMtFullY,1,3,nparams);
        for (int k = 0; k < _N_FIGURES_; k++) {
            fun[j][k]->SetParameter(0,teffParameter[j]);
            fun[j][k]->SetParameter(1,normParameter[j]);
        }
    }*/

    /*TF1* fun[_N_PARTICLES_];
    fun[0] = new TF1("fun0",fitFunctionMtMidY,0,3,nparams);
    fun[1] = new TF1("fun1",fitFunctionMtFullY,0,3,nparams);

    for(int j = 0; j < _N_FIGURES_; j++) {
        fun[j]->SetParameter(0,0.00001);
        fun[j]->SetParameter(1,2e3);
        fun[j]->SetParLimits(0,0.00001,0.001);
        fun[j]->SetParLimits(1,100,10000);
    }*/

    TF1* fun[_N_PARTICLES_];
    for(int j = 0; j < _N_PARTICLES_; j++) {
        fun[j] = new TF1(Form("fun%d",j),fitFunctionMt,0,3,nparams);
        fun[j]->SetParameter(0,teffParameter[j]);
        fun[j]->SetParameter(1,normParameter[j]);
    }
    

   

    //FIITING

    Double_t fitMin[_N_CASES_][_N_PARTICLES_] = {{1.15,0.6,0.2},
                                                {1.15,0.6,0.2},
                                                {1.15,0.6,0.2}};

    Double_t fitMax[_N_CASES_][_N_PARTICLES_] = {{1.7,0.9,0.55},
                                                {1.7,0.9,0.55},
                                                {1.7,0.9,0.55}};

    TF1* funCop;
    Double_t minRange[_N_PARTICLES_] ={0,0,0};
    Double_t maxRange[_N_PARTICLES_] ={2,1.2,0.8};

    Double_t tEff[_N_PARTICLES_][_N_FIGURES_];
    Double_t tEffEr[_N_PARTICLES_][_N_FIGURES_];
    Double_t mXValues[_N_PARTICLES_] = {0.9395653,0.4936770,0.1395699};
    Double_t XErr[_N_PARTICLES_] = {0,0,0};

    //gROOT -> SetBatch(kTRUE);
    
    TFile* fileOut = new TFile(Form("outputPart/outFitMtPart%s.root",casesNames[i].Data()),"RECREATE");
	fileOut->cd();

    //TCanvas* can[_N_PARTICLES_][_N_FIGURES_];
    //for(int i = 0; i < _N_CASES_; i++) {
    //for(int i = 0; i < 1; i++) {
        for(int j = 0; j < _N_PARTICLES_; j++) { 
            for(int k = 0; k < _N_FIGURES_; k++) {
                //can[j][k] = new TCanvas(Form("%s%s%s",casesNames[i].Data(),particlesTitles[j].Data(),figuresNames[k].Data()),Form("%s%s%s",casesNames[i].Data(),particlesTitles[j].Data(),figuresNames[k].Data()),1000,1000);
                histMt[j][k]->Fit(fun[j],"","",fitMin[i][j],fitMax[i][j]);
//Form("fun%i%i",j,k)
                funCop = (TF1*)fun[j]->Clone("funCop");
                funCop->SetRange(mXValues[j],maxRange[j]);
                funCop->SetLineStyle(2);
                histMt[j][k]->GetListOfFunctions()->Add(funCop);

                //can[j][k]->cd();
                //histMt[j]->Draw();
                //histMt[j]->GetFunction(Form("fun%i",j))->Draw("same");
                //funCop[j]->Draw("same");
                histMt[j][k]->Write();

                tEff[j][k] = fun[j]->GetParameter(0);
                tEffEr[j][k] = fun[j]->GetParError(0);
            }
        }
    //}

    //TEFF(M) PLOTS
    for(int j = 0; j < _N_PARTICLES_; j ++) {
        cout<<particlesTitles[j].Data()<<endl;
        cout<<"M: "<<mXValues[j]<<endl;
        for(int k = 0; k < _N_FIGURES_; k++) {
            cout<<figuresNames[k].Data()<<endl;
            cout<<"Teff: "<<tEff[j][k]<<endl;
            cout<<endl;
        }
    }
    /*TGraphErrors *grTeff = new TGraphErrors(_N_PARTICLES_,mXValues,tEff,XErr,tEffEr);
    TCanvas* can2 = new TCanvas("chuj","chuj",1000,1000);
    can2->cd();
    //grTeff->SetMarkerStyle(20);
    grTeff->SetMarkerSize(1);
    grTeff->Draw("ap");*/





    fileOut->Save();
	fileOut->Close();
    //gROOT -> SetBatch(kFALSE);

}

void compParticlesMt() {

    //createHistograms(0);
    createHistograms(1);
    //createHistograms(2);

}

