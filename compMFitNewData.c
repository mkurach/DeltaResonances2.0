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
#include "GPlotHandler.h"

#define _N_CASES_ 3
#define _N_PAIRS_ 2
#define _N_FIGURES_ 2
#define _N_CANVASES_ 4



void makePaveText(TVirtualPad* can, TString text, double x1, double y1, double x2, double y2, double size, Int_t color=kBlack) {
    TPaveText *description=new TPaveText(x1,y1,x2,y2,"NDC");
    description->SetLineWidth(0);
  description->AddText(text);
  description->SetTextSize(size);
  description->SetBorderSize(0.0);
  description->SetTextFont(42);
  description->SetTextColor(color);
  description->SetFillColor(0);
  description->SetFillStyle(0);
  description->SetTextAlign(13);
      can->cd();
    description->Draw("same");
}

TString pairsTitles[_N_PAIRS_] = {"PiPlusP","PiMinusP"};
TString pairsNames[_N_PAIRS_] = {"#pi^{+}p","#pi^{-}p"};
TString figuresNames[_N_FIGURES_] = {"Monitz","Manley"};
TString canvasesNames[_N_CANVASES_] = {"SetDeltaNoExp","ParDeltaNoExp","SetDeltaExp","ParDeltaExp"};
TString canvasesNamesNice[_N_CANVASES_] = {"#delta^{2} set, no exp","#delta^{2} as parameter, no exp","#delta^{2} set, with exp","#delta^{2} as parameter, with exp"};

int colors[_N_FIGURES_] = {kBlue, kGreen+2};
int lines[_N_FIGURES_] = {2,4};
Float_t XMin = 1.08;
Float_t XMax = 1.6;

Double_t fitMin = 1.125;
Double_t fitMax = 1.25;


Double_t mPi = 0.1395699; //GeV
Double_t mProt = 0.9382720;
Double_t delta2Mon = 0.09; //GeV^2
Double_t delta2Man = 0.09;//1.0
Double_t delta2K = 0.0516; //GeV^2

Double_t k(Double_t x) {
    if(x != 0)
        return (x*x-(mPi+mProt)*(mPi+mProt))*(x*x-(mPi-mProt)*(mPi-mProt))/(4.0*x*x);
    else 
        return 0;
}


Double_t monitzSetDeltaNoExp(Double_t *x, Double_t *par) {// 0 - norm, 1 - m_delta, 2 - gamma_delta 
    Double_t k2 = k(x[0]);

    Double_t gamma = par[2]*par[1]/x[0]*TMath::Power(TMath::Sqrt(k2)/TMath::Sqrt(delta2K),3)*TMath::Power((delta2K+delta2Mon)/(k2+delta2Mon),2);
    return par[0]*2.0/3.14*x[0]*x[0]*gamma/(TMath::Power(x[0]*x[0]-par[1]*par[1],2)+x[0]*x[0]*gamma*gamma);

}

Double_t monitzParDeltaNoExp(Double_t *x, Double_t *par) {// 0 - norm, 1 - m_delta, 2 - gamma_delta, 3 - delta 
    Double_t k2 = k(x[0]);

    Double_t gamma = par[2]*par[1]/x[0]*TMath::Power(TMath::Sqrt(k2)/TMath::Sqrt(delta2K),3)*TMath::Power((delta2K+par[3])/(k2+par[3]),2);
    return par[0]*2.0/3.14*x[0]*x[0]*gamma/(TMath::Power(x[0]*x[0]-par[1]*par[1],2)+x[0]*x[0]*gamma*gamma);

}

Double_t monitzSetDeltaExp(Double_t *x, Double_t *par) {// 0 - norm, 1 - m_delta, 2 - gamma_delta, 3 - temp
    Double_t k2 = k(x[0]);

    Double_t gamma = par[2]*par[1]/x[0]*TMath::Power(TMath::Sqrt(k2)/TMath::Sqrt(delta2K),3)*TMath::Power((delta2K+delta2Mon)/(k2+delta2Mon),2);
    return par[0]*2.0/3.14*TMath::Power(x[0]*par[3],1.5)*TMath::Exp(-x[0]/par[3])*x[0]*x[0]*gamma/(TMath::Power(x[0]*x[0]-par[1]*par[1],2)+x[0]*x[0]*gamma*gamma);

}

Double_t monitzParDeltaExp(Double_t *x, Double_t *par) {// 0 - norm, 1 - m_delta, 2 - gamma_delta, 3 - delta , 4 - temp
    Double_t k2 = k(x[0]);

    Double_t gamma = par[2]*par[1]/x[0]*TMath::Power(TMath::Sqrt(k2)/TMath::Sqrt(delta2K),3)*TMath::Power((delta2K+par[3])/(k2+par[3]),2);
    return par[0]*2.0/3.14*TMath::Power(x[0]*par[4],1.5)*TMath::Exp(-x[0]/par[4])*x[0]*x[0]*gamma/(TMath::Power(x[0]*x[0]-par[1]*par[1],2)+x[0]*x[0]*gamma*gamma);

}


Double_t manleySetDeltaNoExp(Double_t *x, Double_t *par) {// 0 - norm, 1 - m_delta, 2 - gamma_delta
    Double_t k2 = k(x[0]);

    Double_t gamma = par[2]*par[1]/x[0]*TMath::Power(TMath::Sqrt(k2)/TMath::Sqrt(delta2K),3)*(delta2K+delta2Man)/(k2+delta2Man);
    return par[0]*2.0/3.14*x[0]*x[0]*gamma/(TMath::Power(x[0]*x[0]-par[1]*par[1],2)+x[0]*x[0]*gamma*gamma);

}

Double_t manleyParDeltaNoExp(Double_t *x, Double_t *par) {// 0 - norm, 1 - m_delta, 2 - gamma_delta, 3 - delta
    Double_t k2 = k(x[0]);

    Double_t gamma = par[2]*par[1]/x[0]*TMath::Power(TMath::Sqrt(k2)/TMath::Sqrt(delta2K),3)*(delta2K+par[3])/(k2+par[3]);
    return par[0]*2.0/3.14*x[0]*x[0]*gamma/(TMath::Power(x[0]*x[0]-par[1]*par[1],2)+x[0]*x[0]*gamma*gamma);

}

Double_t manleySetDeltaExp(Double_t *x, Double_t *par) {// 0 - norm, 1 - m_delta, 2 - gamma_delta, 3 - temp
    Double_t k2 = k(x[0]);

    Double_t gamma = par[2]*par[1]/x[0]*TMath::Power(TMath::Sqrt(k2)/TMath::Sqrt(delta2K),3)*(delta2K+delta2Man)/(k2+delta2Man);
    return par[0]*2.0/3.14*TMath::Power(x[0]*par[3],1.5)*TMath::Exp(-x[0]/par[3])*x[0]*x[0]*gamma/(TMath::Power(x[0]*x[0]-par[1]*par[1],2)+x[0]*x[0]*gamma*gamma);

}

Double_t manleyParDeltaExp(Double_t *x, Double_t *par) {// 0 - norm, 1 - m_delta, 2 - gamma_delta, 3 - delta, 4 - temp
    Double_t k2 = k(x[0]);

    Double_t gamma = par[2]*par[1]/x[0]*TMath::Power(TMath::Sqrt(k2)/TMath::Sqrt(delta2K),3)*(delta2K+par[3])/(k2+par[3]);
    return par[0]*2.0/3.14*TMath::Power(x[0]*par[4],1.5)*TMath::Exp(-x[0]/par[4])*x[0]*x[0]*gamma/(TMath::Power(x[0]*x[0]-par[1]*par[1],2)+x[0]*x[0]*gamma*gamma);

}

Double_t fitFunctionM(Double_t *x, Double_t *par) {

    /*Double_t Ecm = 2.42; //GeV
    Double_t mN = 0.93;
    Double_t q = TMath::Sqrt(2*Ecm*mN+TMath::Power(Ecm,2));
    Double_t mu = 0.18;*/
    if(x[0] >= 1 && x[0] <= 1.6)
        //return TMath::Power(q,3)/(TMath::Power(q,3)+TMath::Power(mu,3))*par[0]/(1+4*TMath::Power((x[0]-par[1])/par[2],2));
        return par[0]/(1+4*TMath::Power((x[0]-par[1])/par[2],2));
    else    
        return 0;
}

/*Double_t gammaFunctionMon(Double_t *x, Double_t *par) {
    Double_t k2;
    if(x[0] != 0)
        k2 = (x[0]*x[0]-(mPi+mProt)*(mPi+mProt))*(x[0]*x[0]-(mPi-mProt)*(mPi-mProt))/(4.0*x[0]*x[0]);
    else 
        k2 = 0;

    return par[0]*par[2]*par[1]/x[0]*TMath::Power(TMath::Sqrt(k2)/TMath::Sqrt(delta2K),3)*TMath::Power((delta2K+delta2Mon)/(k2+delta2Mon),2);
   //return par[0]*par[2]*par[1]/x[0]*TMath::Power(TMath::Sqrt(k2)/TMath::Sqrt(delta2K),3)*TMath::Power((delta2K+par[3])/(k2+par[3]),2);
}

Double_t gammaFunctionMan(Double_t *x, Double_t *par) {
        Double_t k2;
    if(x[0] != 0)
        k2 = (x[0]*x[0]-(mPi+mProt)*(mPi+mProt))*(x[0]*x[0]-(mPi-mProt)*(mPi-mProt))/(4.0*x[0]*x[0]);
    else 
        k2 = 0;

    return par[2]*par[1]/x[0]*TMath::Power(TMath::Sqrt(k2)/TMath::Sqrt(delta2K),3)*(delta2K+delta2Man)/(k2+delta2Man);
    //return par[2]*par[1]/x[0]*TMath::Power(TMath::Sqrt(k2)/TMath::Sqrt(delta2K),3)*(delta2K+par[3])/(k2+par[3]);
}*/


void setDeltaNoExp(TH1D* dataHist[],TF1* dataFunction[],TF1* dataFunctionExtra[]) { // 0 - norm, 1 - m_delta, 2 - gamma_delta

    TF1* fun[_N_FIGURES_];
    fun[0] = new TF1(figuresNames[0].Data(),monitzSetDeltaNoExp,XMin,XMax,3);
    fun[1] = new TF1(figuresNames[1].Data(),manleySetDeltaNoExp,XMin,XMax,3);

    fun[0]->SetLineColor(colors[0]);
    fun[1]->SetLineColor(colors[1]);

    fun[0]->SetLineStyle(lines[0]);
    fun[1]->SetLineStyle(lines[1]);

    for(int i = 0; i < _N_FIGURES_; i++) {
        fun[i]->SetParameters(1.0,1.232,0.117);
        fun[i]->SetParLimits(1,1.1,1.4);
        //fun[i]->SetParLimits(2,0.1,0.2);
    }

    TFile* fileOut = new TFile("outputNewData/outSetDeltaNoExp.root","RECREATE");

    TCanvas* can[_N_PAIRS_];

    Double_t m0[_N_PAIRS_][_N_FIGURES_];
    Double_t gamma0[_N_PAIRS_][_N_FIGURES_];
    Double_t m0Err[_N_PAIRS_][_N_FIGURES_];
    Double_t gamma0Err[_N_PAIRS_][_N_FIGURES_];

    for(int i = 0; i < _N_PAIRS_; i++) {
        cout<<"\t\t"<<pairsTitles[i].Data()<<endl;
        cout<<"\t\t\t"<<canvasesNames[0].Data()<<endl;
        can[i] = new TCanvas(Form("%s%s",pairsTitles[i].Data(),canvasesNames[0].Data()),Form("%s%s",pairsTitles[i].Data(),canvasesNames[0].Data()),1000,1000);

        for(int j = 0; j < _N_FIGURES_; j++) {
            cout<<"\t\t\t\t"<<figuresNames[j].Data()<<endl;
            dataHist[i]->Fit(fun[j],"M+","M+",fitMin,fitMax);
            
            fun[j]->SetRange(XMin,XMax);
            dataHist[i]->GetListOfFunctions()->Add(fun[j]);

            m0[i][j] = fun[j]->GetParameter(1);
            gamma0[i][j] = fun[j]->GetParameter(2);
            m0Err[i][j] = fun[j]->GetParError(1);
            gamma0Err[i][j] = fun[j]->GetParError(2);

        }

        can[i]->cd();
        dataHist[i]->GetListOfFunctions()->Add(dataFunctionExtra[i]);
        dataHist[i]->Draw();


        makePaveText(can[i], canvasesNamesNice[0].Data(),0.45, 0.85, 0.7, 0.9, 0.04);
        makePaveText(can[i], "Monitz", 0.6, 0.75, 0.8, 0.8, 0.04, colors[0]);
        makePaveText(can[i], "Manley", 0.6, 0.7, 0.8, 0.75, 0.04, colors[1]);
        makePaveText(can[i], "Breit-Wigner", 0.6, 0.65, 0.8, 0.7 , 0.04, kRed);
        makePaveText(can[i], pairsNames[i].Data(),0.7, 0.5, 0.8, 0.55, 0.05);

        fileOut->cd();
        dataHist[i]->Write();
        can[i]->SaveAs(Form("outputNewData/%s%s.png",pairsTitles[i].Data(),canvasesNames[0].Data()));


    }

    fileOut->Save();
    fileOut->Close();


    ofstream fileTxt;
    fileTxt.open(Form("outputNewData/%sFitResults.txt",canvasesNames[0].Data()));


    for(int i = 0; i < _N_PAIRS_; i++) {
        fileTxt<<pairsTitles[i].Data()<<"\n";
        for(int j = 0; j < _N_FIGURES_; j++) {
            fileTxt<<"\t"<<figuresNames[j].Data()<<"\n";
            fileTxt<<Form("\t\tM=%.3f (%.3f)\t#Gamma = %.3f (%.3f)  \t#delta^{2} = %.2f",m0[i][j],m0Err[i][j],gamma0[i][j],gamma0Err[i][j],delta2Mon)<<"\n";

        }

    }
    fileTxt<<"\tB-W\n\t\t M = 1.215\t#Gamma = 0.117";

    fileTxt.close();






}
  
void parDeltaNoExp(TH1D* dataHist[],TF1* dataFunction[],TF1* dataFunctionExtra[]) { // 0 - norm, 1 - m_delta, 2 - gamma_delta, 3 - delta


    TF1* fun[_N_FIGURES_];

    fun[0] = new TF1(figuresNames[0].Data(),monitzParDeltaNoExp,XMin,XMax,4);
    fun[1] = new TF1(figuresNames[1].Data(),manleyParDeltaNoExp,XMin,XMax,4);

    fun[0]->SetLineColor(colors[0]);
    fun[1]->SetLineColor(colors[1]);

    fun[0]->SetLineStyle(lines[0]);
    fun[1]->SetLineStyle(lines[1]);

    for(int j = 0; j < _N_FIGURES_; j++) {
        fun[j]->SetParameters(1.0,1.232,0.117,delta2Mon);
        fun[j]->SetParLimits(1,1.1,1.4);
        fun[j]->SetParLimits(3,0.01,1.0);
        //fun[i]->SetParLimits(2,0.1,0.2);
    }


    TFile* fileOut = new TFile("outputNewData/outParDeltaNoExp.root","RECREATE");

    TCanvas* can[_N_PAIRS_];

    Double_t m0[_N_PAIRS_][_N_FIGURES_];
    Double_t gamma0[_N_PAIRS_][_N_FIGURES_];
    Double_t m0Err[_N_PAIRS_][_N_FIGURES_];
    Double_t gamma0Err[_N_PAIRS_][_N_FIGURES_];
    Double_t delta[_N_PAIRS_][_N_FIGURES_];
    Double_t deltaErr[_N_PAIRS_][_N_FIGURES_];

    for(int i = 0; i < _N_PAIRS_; i++) {
        cout<<"\t\t"<<pairsTitles[i].Data()<<endl;
        cout<<"\t\t\t"<<canvasesNames[1].Data()<<endl;
        can[i] = new TCanvas(Form("%s%s",pairsTitles[i].Data(),canvasesNames[1].Data()),Form("%s%s",pairsTitles[i].Data(),canvasesNames[1].Data()),1000,1000);


        for(int j = 0; j < _N_FIGURES_; j++) {
            cout<<"\t\t\t\t"<<figuresNames[j].Data()<<endl;
            dataHist[i]->Fit(fun[j],"M+","M+",fitMin,fitMax);

            fun[j]->SetRange(XMin,XMax);
            dataHist[i]->GetListOfFunctions()->Add(fun[j]);

            m0[i][j] = fun[j]->GetParameter(1);
            gamma0[i][j] = fun[j]->GetParameter(2);
            m0Err[i][j] = fun[j]->GetParError(1);
            gamma0Err[i][j] = fun[j]->GetParError(2);
            delta[i][j] = fun[j]->GetParameter(3);
            deltaErr[i][j] = fun[j]->GetParError(3);


        }

        can[i]->cd();
        dataHist[i]->GetListOfFunctions()->Add(dataFunctionExtra[i]);
        dataHist[i]->Draw();

        makePaveText(can[i], canvasesNamesNice[1].Data(),0.45, 0.85, 0.7, 0.9, 0.04);
        makePaveText(can[i], "Monitz", 0.6, 0.75, 0.8, 0.8, 0.04, colors[0]);
        makePaveText(can[i], "Manley", 0.6, 0.7, 0.8, 0.75, 0.04, colors[1]);
        makePaveText(can[i], "Breit-Wigner", 0.6, 0.65, 0.8, 0.7 , 0.04, kRed);
        makePaveText(can[i], pairsNames[i].Data(),0.7, 0.5, 0.8, 0.55, 0.05);

        fileOut->cd();
        dataHist[i]->Write();
        can[i]->SaveAs(Form("outputNewData/%s%s.png",pairsTitles[i].Data(),canvasesNames[1].Data()));


        


    }

    fileOut->Save();
    fileOut->Close();

    ofstream fileTxt;
    fileTxt.open(Form("outputNewData/%sFitResults.txt",canvasesNames[1].Data()));


    for(int i = 0; i < _N_PAIRS_; i++) {
        fileTxt<<pairsTitles[i].Data()<<"\n";
        for(int j = 0; j < _N_FIGURES_; j++) {
            fileTxt<<"\t"<<figuresNames[j].Data()<<"\n";
            fileTxt<<Form("\t\t  M = %.3f (%.3f)\t#Gamma = %.3f (%.3f)\t#delta^{2} = %.2f (%.2f)",m0[i][j],m0Err[i][j],gamma0[i][0],gamma0Err[i][j],delta[i][j],deltaErr[i][j])<<"\n";


        }

    }
    fileTxt<<"\tB-W\n\t\t M = 1.215\t#Gamma = 0.117";

    fileTxt.close();

}
   

void setDeltaExp(TH1D* dataHist[],TF1* dataFunction[],TF1* dataFunctionExtra[]){ // 0 - norm, 1 - m_delta, 2 - gamma_delta, 3 - temp

    TF1* fun[_N_FIGURES_];
    fun[0] = new TF1(figuresNames[0].Data(),monitzSetDeltaExp,XMin,XMax,4);
    fun[1] = new TF1(figuresNames[1].Data(),manleySetDeltaExp,XMin,XMax,4);

    fun[0]->SetLineColor(colors[0]);
    fun[1]->SetLineColor(colors[1]);

    fun[0]->SetLineStyle(lines[0]);
    fun[1]->SetLineStyle(lines[1]);

    for(int i = 0; i < _N_FIGURES_; i++) {
        fun[i]->SetParameters(10000,1.232,0.117,0.5);
        fun[i]->SetParLimits(1,1.1,1.5);
        fun[i]->SetParLimits(3,0.01,0.5);
        //fun[i]->SetParLimits(2,0.1,1.0);
        //fun[i]->SetParLimits(3,0.01,1.0);
    }

    TFile* fileOut = new TFile("outputNewData/outSetDeltaExp.root","RECREATE");

    TCanvas* can[_N_PAIRS_];

    Double_t m0[_N_PAIRS_][_N_FIGURES_];
    Double_t gamma0[_N_PAIRS_][_N_FIGURES_];
    Double_t m0Err[_N_PAIRS_][_N_FIGURES_];
    Double_t gamma0Err[_N_PAIRS_][_N_FIGURES_];
    Double_t temp[_N_PAIRS_][_N_FIGURES_];
    Double_t tempErr[_N_PAIRS_][_N_FIGURES_];

    for(int i = 0; i < _N_PAIRS_; i++) {
        cout<<"\t\t"<<pairsTitles[i].Data()<<endl;
        cout<<"\t\t\t"<<canvasesNames[2].Data()<<endl;
        can[i] = new TCanvas(Form("%s%s",pairsTitles[i].Data(),canvasesNames[2].Data()),Form("%s%s",pairsTitles[i].Data(),canvasesNames[2].Data()),1000,1000);

        for(int j = 0; j < _N_FIGURES_; j++) {
            cout<<"\t\t\t\t"<<figuresNames[j].Data()<<endl;
            dataHist[i]->Fit(fun[j],"M+","M+",fitMin,fitMax);

            fun[j]->SetRange(XMin,XMax);
            dataHist[i]->GetListOfFunctions()->Add(fun[j]);

            m0[i][j] = fun[j]->GetParameter(1);
            gamma0[i][j] = fun[j]->GetParameter(2);
            m0Err[i][j] = fun[j]->GetParError(1);
            gamma0Err[i][j] = fun[j]->GetParError(2);
            temp[i][j] = fun[j]->GetParameter(3);
            tempErr[i][j] = fun[j]->GetParError(3);

        }

        can[i]->cd();
        dataHist[i]->GetListOfFunctions()->Add(dataFunctionExtra[i]);
        dataHist[i]->Draw();


        makePaveText(can[i], canvasesNamesNice[2].Data(),0.45, 0.85, 0.7, 0.9, 0.04);
        makePaveText(can[i], "Monitz", 0.6, 0.75, 0.8, 0.8, 0.04, colors[0]);
        makePaveText(can[i], "Manley", 0.6, 0.7, 0.8, 0.75, 0.04, colors[1]);
        makePaveText(can[i], "Breit-Wigner", 0.6, 0.65, 0.8, 0.7 , 0.04, kRed);
        makePaveText(can[i], pairsNames[i].Data(),0.7, 0.5, 0.8, 0.55, 0.05);

        fileOut->cd();
        dataHist[i]->Write();
        can[i]->SaveAs(Form("outputNewData/%s%s.png",pairsTitles[i].Data(),canvasesNames[2].Data()));


    }

    fileOut->Save();
    fileOut->Close();


    ofstream fileTxt;
    fileTxt.open(Form("outputNewData/%sFitResults.txt",canvasesNames[2].Data()));


    for(int i = 0; i < _N_PAIRS_; i++) {
        fileTxt<<pairsTitles[i].Data()<<"\n";
        for(int j = 0; j < _N_FIGURES_; j++) {
            fileTxt<<"\t"<<figuresNames[j].Data()<<"\n";
            fileTxt<<Form("\t\tM=%.3f (%.3f)\t#Gamma = %.3f (%.3f)\t#delta^{2} = %.2f\tT = %.3f (%.3f) MeV",m0[i][j],m0Err[i][j],gamma0[i][j],gamma0Err[i][j],delta2Mon,temp[i][j]*1e3,tempErr[i][j]*1e3)<<"\n";


        }

    }
    fileTxt<<"\tB-W\n\t\t M = 1.215\t#Gamma = 0.117";

    fileTxt.close();



}

void parDeltaExp(TH1D* dataHist[],TF1* dataFunction[],TF1* dataFunctionExtra[]) { // 0 - norm, 1 - m_delta, 2 - gamma_delta, 3 - delta, 4 - temp
    TF1* fun[_N_FIGURES_];
    fun[0] = new TF1(figuresNames[0].Data(),monitzParDeltaExp,XMin,XMax,5);
    fun[1] = new TF1(figuresNames[1].Data(),manleyParDeltaExp,XMin,XMax,5);

    fun[0]->SetLineColor(colors[0]);
    fun[1]->SetLineColor(colors[1]);

    fun[0]->SetLineStyle(lines[0]);
    fun[1]->SetLineStyle(lines[1]);

    for(int i = 0; i < _N_FIGURES_; i++) {
        fun[i]->SetParameters(10.0,1.232,0.117,delta2Mon,0.5);
        fun[i]->SetParLimits(1,1.1,1.5);
        fun[i]->SetParLimits(3,0.01,1.0);
        fun[i]->SetParLimits(4,0.01,0.5);
    }

    TFile* fileOut = new TFile("outputNewData/outParDeltaExp.root","RECREATE");

    TCanvas* can[_N_PAIRS_];

    Double_t m0[_N_PAIRS_][_N_FIGURES_];
    Double_t gamma0[_N_PAIRS_][_N_FIGURES_];
    Double_t m0Err[_N_PAIRS_][_N_FIGURES_];
    Double_t gamma0Err[_N_PAIRS_][_N_FIGURES_];
    Double_t delta[_N_PAIRS_][_N_FIGURES_];
    Double_t deltaErr[_N_PAIRS_][_N_FIGURES_];
    Double_t temp[_N_PAIRS_][_N_FIGURES_];
    Double_t tempErr[_N_PAIRS_][_N_FIGURES_];



    for(int i = 0; i < _N_PAIRS_; i++) {
        cout<<"\t\t"<<pairsTitles[i].Data()<<endl;
        cout<<"\t\t\t"<<canvasesNames[3].Data()<<endl;
        can[i] = new TCanvas(Form("%s%s",pairsTitles[i].Data(),canvasesNames[3].Data()),Form("%s%s",pairsTitles[i].Data(),canvasesNames[3].Data()),1000,1000);

        for(int j = 0; j < _N_FIGURES_; j++) {
            cout<<"\t\t\t\t"<<figuresNames[j].Data()<<endl;
            dataHist[i]->Fit(fun[j],"M+","M+",fitMin,fitMax);

            fun[j]->SetRange(XMin,XMax);
            dataHist[i]->GetListOfFunctions()->Add(fun[j]);

            m0[i][j] = fun[j]->GetParameter(1);
            gamma0[i][j] = fun[j]->GetParameter(2);
            m0Err[i][j] = fun[j]->GetParError(1);
            gamma0Err[i][j] = fun[j]->GetParError(2);
            delta[i][j] = fun[j]->GetParameter(3);
            deltaErr[i][j] = fun[j]->GetParError(3);
            temp[i][j] = fun[j]->GetParameter(4);
            tempErr[i][j] = fun[j]->GetParError(4);


        }

        can[i]->cd();
        dataHist[i]->GetListOfFunctions()->Add(dataFunctionExtra[i]);
        dataHist[i]->Draw();


        makePaveText(can[i], canvasesNamesNice[3].Data(),0.45, 0.85, 0.7, 0.9, 0.04);
        makePaveText(can[i], "Monitz", 0.6, 0.75, 0.8, 0.8, 0.04, colors[0]);
        makePaveText(can[i], "Manley", 0.6, 0.7, 0.8, 0.75, 0.04, colors[1]);
        makePaveText(can[i], "Breit-Wigner", 0.6, 0.65, 0.8, 0.7 , 0.04, kRed);
        makePaveText(can[i], pairsNames[i].Data(),0.7, 0.5, 0.8, 0.55, 0.05);

        fileOut->cd();
        dataHist[i]->Write();
        can[i]->SaveAs(Form("outputNewData/%s%s.png",pairsTitles[i].Data(),canvasesNames[3].Data()));


    }

    fileOut->Save();
    fileOut->Close();

    ofstream fileTxt;
    fileTxt.open(Form("outputNewData/%sFitResults.txt",canvasesNames[3].Data()));


    for(int i = 0; i < _N_PAIRS_; i++) {
        fileTxt<<pairsTitles[i].Data()<<"\n";
        for(int j = 0; j < _N_FIGURES_; j++) {
            fileTxt<<"\t"<<figuresNames[j].Data()<<"\n";
            fileTxt<<Form("\t\tM=%.3f (%.3f)\t#Gamma = %.3f (%.3f)\t#delta^{2} = %.2f\tT = %.2f (%.2f) MeV",m0[i][j],m0Err[i][j],gamma0[i][j],gamma0Err[i][j],delta[i][j],temp[i][j]*1e3,tempErr[i][j]*1e3)<<"\n";

        }

    }
    fileTxt<<"\tB-W\n\t\t M = 1.215\t#Gamma = 0.117";

    fileTxt.close();

}


void compCanvases() {

    TH1D* hist[_N_PAIRS_][_N_CANVASES_];

    TFile* fileIn[_N_CANVASES_];
    for(int i = 0; i < _N_CANVASES_; i++) {
        fileIn[i] = new TFile(Form("outputNewData/out%s.root",canvasesNames[i].Data()));
        for(int j = 0; j < _N_PAIRS_; j++) {
            hist[j][i] = (TH1D*)fileIn[i]->Get(Form("%sDataHist",pairsTitles[j].Data()));
            hist[j][i]->SetTitle(canvasesNamesNice[i].Data());
        }
    }

    TCanvas* canOut[_N_PAIRS_];
    TFile* fileOut = new TFile("outputNewData/outComp.root","RECREATE");
    for(int i = 0; i < _N_PAIRS_; i++) {
        canOut[i] = new TCanvas(Form("%sComp",pairsTitles[i].Data()),Form("%sComp",pairsTitles[i].Data()),4000,4000);
        canOut[i]->Divide(2,2);
        for(int j = 0; j < _N_CANVASES_; j++) {
            canOut[i]->cd(j+1);
            hist[i][j]->Draw();
        }

        fileOut->cd();
        canOut[i]->Write();
        //canOut[i]->SaveAs(Form("outputNewData/%sComp.png",pairsTitles[i].Data()));



    }

    fileOut->Close();
    fileOut->Save();







}







void compMFitNewData() {
    
    //gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);


    //Experiment histograms

    TH1D *dataHist[_N_PAIRS_];

    TFile* filePlus = new TFile("/u/mkurach/figures_with_data/moje/exp/massPlus.root");
    TFile* fileMinus = new TFile("/u/mkurach/figures_with_data/moje/exp/massMinus.root");

    dataHist[0] = (TH1D*)filePlus->Get("h_circ_int_mass_mult0_case1__289")->Clone(Form("%sDataHist",pairsTitles[0].Data()));
    dataHist[1] = (TH1D*)fileMinus->Get("h_circ_int_mass_mult0_case2__302")->Clone(Form("%sDataHist",pairsTitles[1].Data()));  
    TString YTitles[_N_PAIRS_] = {"M_{#pi^{+} p} (GeV/c^{2})","M_{#pi^{-} p} (GeV/c^{2})"};

    for(int i = 0; i < _N_PAIRS_; i++) {
        dataHist[i]->SetTitle("");
        dataHist[i]->GetXaxis()->SetTitle(YTitles[i].Data());
        dataHist[i]->GetYaxis()->SetTitle("dN/dM (GeV/c^{2})^{-1}");
        dataHist[i]->GetYaxis()->SetRangeUser(0,25);
        dataHist[i]->GetYaxis()->SetTitleSize(0.04);
        dataHist[i]->GetXaxis()->SetTitleSize(0.04);
        dataHist[i]->GetYaxis()->SetLabelSize(0.04);
        dataHist[i]->GetXaxis()->SetLabelSize(0.04);

        dataHist[i]->GetXaxis()->SetRangeUser(1.05,XMax);
        dataHist[i]->SetMarkerSize(1.3);
    }

    dataHist[1]->SetMarkerStyle(21);
     dataHist[0]->SetMarkerStyle(8);

    //Breit Wigner from experiment

    TF1* dataFunction[_N_PAIRS_];
    filePlus = new TFile("/u/mkurach/figures_with_data/moje/exp/massFitPlus.root");  
    fileMinus = new TFile("/u/mkurach/figures_with_data/moje/exp/massFitMinus.root");
    dataFunction[0] = (TF1*)filePlus->Get("fbwr3")->Clone("fitPlus");
    dataFunction[1] = (TF1*)fileMinus->Get("fbwr3")->Clone("fitMinus");



    TF1* dataFunctionExtra[_N_PAIRS_];
    filePlus = new TFile("/u/mkurach/figures_with_data/moje/exp/massFitExtraPlus.root");  
    fileMinus = new TFile("/u/mkurach/figures_with_data/moje/exp/massFitExtraMinus.root");
    dataFunctionExtra[0] = (TF1*)filePlus->Get("fbwr2")->Clone("fitPlusExtra");
    dataFunctionExtra[1] = (TF1*)fileMinus->Get("fbwr2")->Clone("fitMinusExtra");
    dataFunctionExtra[0]->SetRange(XMin,XMax);
    dataFunctionExtra[1]->SetRange(XMin,XMax);
    dataFunctionExtra[0]->SetLineStyle(1);
    dataFunctionExtra[1]->SetLineStyle(1);



    //setDeltaNoExp(dataHist,dataFunction,dataFunctionExtra);
   //parDeltaNoExp(dataHist,dataFunction,dataFunctionExtra);
    //setDeltaExp(dataHist,dataFunction,dataFunctionExtra);
    //parDeltaExp(dataHist,dataFunction,dataFunctionExtra);

    compCanvases();





}
