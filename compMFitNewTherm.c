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
#define _N_FIGURES_ 3 
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
TString figuresNames[_N_FIGURES_] = {"Monitz","Manley","Breit-Wigner"};
TString casesNames[_N_CASES_] = {"CaseA","CaseB","CaseC"};
TString canvasesNames[_N_CANVASES_] = {"SetDeltaNoExp","ParDeltaNoExp","SetDeltaExp","ParDeltaExp"};
TString canvasesNamesNice[_N_CANVASES_] = {"#delta^{2} set, no exp","#delta^{2} as parameter, no exp","#delta^{2} set, with exp","#delta^{2} as parameter, with exp"};

int colors[_N_FIGURES_] = {kBlue, kGreen+2,kRed};
int lines[_N_FIGURES_] = {2,4,1};
Int_t XBins = 1000;
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
    //if(x[0] >= 1 && x[0] <= 1.6)
        //return TMath::Power(q,3)/(TMath::Power(q,3)+TMath::Power(mu,3))*par[0]/(1+4*TMath::Power((x[0]-par[1])/par[2],2));
        return par[0]/(1+4*TMath::Power((x[0]-par[1])/par[2],2));
    //else    
        //return 0;
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


void setDeltaNoExp(TH1D* dataHist[][_N_PAIRS_]) { // 0 - norm, 1 - m_delta, 2 - gamma_delta

    TF1* fun[_N_FIGURES_];
    fun[0] = new TF1(figuresNames[0].Data(),monitzSetDeltaNoExp,XMin,XMax,3);
    fun[1] = new TF1(figuresNames[1].Data(),manleySetDeltaNoExp,XMin,XMax,3);
    fun[2] = new TF1(figuresNames[2].Data(),fitFunctionM,XMin,XMax,3);

    for(int i = 0; i < _N_FIGURES_; i++) {
        fun[i]->SetLineColor(colors[i]);
        fun[i]->SetLineStyle(lines[i]);

    }

    for(int i = 0; i < _N_FIGURES_; i++) {
        fun[i]->SetParameters(0.1,1.232,0.117);
        fun[i]->SetParLimits(1,1.1,1.4);
        fun[i]->SetParLimits(2,0.0,1.0);
    }

    TFile* fileOut = new TFile("outputMNewTherm/outSetDeltaNoExp.root","RECREATE");

    TCanvas* can[_N_CASES_][_N_PAIRS_];

    Double_t m0[_N_CASES_][_N_PAIRS_][_N_FIGURES_];
    Double_t gamma0[_N_CASES_][_N_PAIRS_][_N_FIGURES_];
    Double_t m0Err[_N_CASES_][_N_PAIRS_][_N_FIGURES_];
    Double_t gamma0Err[_N_CASES_][_N_PAIRS_][_N_FIGURES_];


    for(int k = 0; k < _N_CASES_; k++) {
        cout<<"\t"<<casesNames[k].Data()<<endl;
        for(int i = 0; i < _N_PAIRS_; i++) {
            cout<<"\t\t"<<pairsTitles[i].Data()<<endl;
            cout<<"\t\t\t"<<canvasesNames[0].Data()<<endl;
            can[k][i] = new TCanvas(Form("%s%s%s",casesNames[k].Data(),pairsTitles[i].Data(),canvasesNames[0].Data()),Form("%s%s%s",casesNames[k].Data(),pairsTitles[i].Data(),canvasesNames[0].Data()),1000,1000);

            for(int j = 0; j < _N_FIGURES_; j++) {
                cout<<"\t\t\t\t"<<figuresNames[j].Data()<<endl;
                dataHist[k][i]->Fit(fun[j],"M+","M+",fitMin,fitMax);
                
                fun[j]->SetRange(XMin,XMax);
                dataHist[k][i]->GetListOfFunctions()->Add(fun[j]);

                m0[k][i][j] = fun[j]->GetParameter(1);
                gamma0[k][i][j] = fun[j]->GetParameter(2);
                m0Err[k][i][j] = fun[j]->GetParError(1);
                gamma0Err[k][i][j] = fun[j]->GetParError(2);

            }

            can[k][i]->cd();
            dataHist[k][i]->Draw();


            makePaveText(can[k][i], canvasesNamesNice[0].Data(),0.45, 0.85, 0.7, 0.9, 0.04);
            makePaveText(can[k][i], "Monitz", 0.6, 0.75, 0.8, 0.8, 0.04, colors[0]);
            makePaveText(can[k][i], "Manley", 0.6, 0.7, 0.8, 0.75, 0.04, colors[1]);
            makePaveText(can[k][i], "Breit-Wigner", 0.6, 0.65, 0.8, 0.7 , 0.04, kRed);
            makePaveText(can[k][i], pairsNames[i].Data(),0.7, 0.5, 0.8, 0.55, 0.05);

            fileOut->cd();
            dataHist[k][i]->Write();
            can[k][i]->SaveAs(Form("outputMNewTherm/%s%s%s.png",casesNames[k].Data(),pairsTitles[i].Data(),canvasesNames[0].Data()));


        }
    }

    fileOut->Save();
    fileOut->Close();


    ofstream fileTxt;
    fileTxt.open(Form("outputMNewTherm/%sFitResults.txt",canvasesNames[0].Data()));

    for(int k = 0; k < _N_CASES_; k++) {
        fileTxt<<casesNames[k].Data()<<"\n";
        for(int i = 0; i < _N_PAIRS_; i++) {
            fileTxt<<"\t"<<pairsTitles[i].Data()<<"\n";
            for(int j = 0; j < _N_FIGURES_; j++) {
                fileTxt<<"\t\t"<<figuresNames[j].Data()<<"\n";
                fileTxt<<Form("\t\t\tM=%.3f (%.3f)\t#Gamma = %.3f (%.3f) ",m0[k][i][j],m0Err[k][i][j],gamma0[k][i][j],gamma0Err[k][i][j])<<"\n";

            }

        }
        fileTxt<<"\n";
    }


    fileTxt.close();






}
  
void parDeltaNoExp(TH1D* dataHist[][_N_PAIRS_]) { // 0 - norm, 1 - m_delta, 2 - gamma_delta, 3 - delta


    TF1* fun[_N_FIGURES_];

    fun[0] = new TF1(figuresNames[0].Data(),monitzParDeltaNoExp,XMin,XMax,4);
    fun[1] = new TF1(figuresNames[1].Data(),manleyParDeltaNoExp,XMin,XMax,4);
    fun[2] = new TF1(figuresNames[2].Data(),fitFunctionM,XMin,XMax,3);

    for(int i = 0; i < _N_FIGURES_; i++) {
        fun[i]->SetLineColor(colors[i]);
        fun[i]->SetLineStyle(lines[i]);
    }

    for(int j = 0; j < _N_FIGURES_; j++) {
        if (j == 2) {
            fun[j]->SetParameters(1.0,1.232,0.117);
        }
        else {
            fun[j]->SetParameters(10.0,1.232,0.117,0.5);
            fun[j]->SetParLimits(3,0.01,1.0);
        }
        fun[j]->SetParLimits(1,1.1,1.4);
        fun[j]->SetParLimits(2,0.0,1.0);

    }


    TFile* fileOut = new TFile("outputMNewTherm/outParDeltaNoExp.root","RECREATE");

    TCanvas* can[_N_CASES_][_N_PAIRS_];

    Double_t m0[_N_CASES_][_N_PAIRS_][_N_FIGURES_];
    Double_t gamma0[_N_CASES_][_N_PAIRS_][_N_FIGURES_];
    Double_t m0Err[_N_CASES_][_N_PAIRS_][_N_FIGURES_];
    Double_t gamma0Err[_N_CASES_][_N_PAIRS_][_N_FIGURES_];
    Double_t delta[_N_CASES_][_N_PAIRS_][_N_FIGURES_];
    Double_t deltaErr[_N_CASES_][_N_PAIRS_][_N_FIGURES_];

    for(int k = 0; k < _N_CASES_; k++) {
        cout<<"\t"<<casesNames[k].Data()<<endl;
        for(int i = 0; i < _N_PAIRS_; i++) {
            cout<<"\t\t"<<pairsTitles[i].Data()<<endl;
            cout<<"\t\t\t"<<canvasesNames[1].Data()<<endl;
            can[k][i] = new TCanvas(Form("%s%s%s",casesNames[k].Data(),pairsTitles[i].Data(),canvasesNames[1].Data()),Form("%s%s%s",casesNames[k].Data(),pairsTitles[1].Data(),canvasesNames[0].Data()),1000,1000);

            for(int j = 0; j < _N_FIGURES_; j++) {
                cout<<"\t\t\t\t"<<figuresNames[j].Data()<<endl;
                dataHist[k][i]->Fit(fun[j],"M+","M+",fitMin,fitMax);
                
                fun[j]->SetRange(XMin,XMax);
                dataHist[k][i]->GetListOfFunctions()->Add(fun[j]);

                m0[k][i][j] = fun[j]->GetParameter(1);
                gamma0[k][i][j] = fun[j]->GetParameter(2);
                m0Err[k][i][j] = fun[j]->GetParError(1);
                gamma0Err[k][i][j] = fun[j]->GetParError(2);

                if (j != 2) {
                    delta[k][i][j] = fun[j]->GetParameter(3);
                    deltaErr[k][i][j] = fun[j]->GetParError(3);
                }

            }

            can[k][i]->cd();
            dataHist[k][i]->Draw();


            makePaveText(can[k][i], canvasesNamesNice[1].Data(),0.45, 0.85, 0.7, 0.9, 0.04);
            makePaveText(can[k][i], "Monitz", 0.6, 0.75, 0.8, 0.8, 0.04, colors[0]);
            makePaveText(can[k][i], "Manley", 0.6, 0.7, 0.8, 0.75, 0.04, colors[1]);
            makePaveText(can[k][i], "Breit-Wigner", 0.6, 0.65, 0.8, 0.7 , 0.04, kRed);
            makePaveText(can[k][i], pairsNames[i].Data(),0.7, 0.5, 0.8, 0.55, 0.05);

            fileOut->cd();
            dataHist[k][i]->Write();
            can[k][i]->SaveAs(Form("outputMNewTherm/%s%s%s.png",casesNames[k].Data(),pairsTitles[i].Data(),canvasesNames[1].Data()));


        }
    }


    fileOut->Save();
    fileOut->Close();

    ofstream fileTxt;
    fileTxt.open(Form("outputMNewTherm/%sFitResults.txt",canvasesNames[1].Data()));

    for(int k = 0; k < _N_CASES_; k++) {
        fileTxt<<casesNames[k].Data()<<"\n";
        for(int i = 0; i < _N_PAIRS_; i++) {
            fileTxt<<"\t"<<pairsTitles[i].Data()<<"\n";
            for(int j = 0; j < _N_FIGURES_; j++) {
                fileTxt<<"\t\t"<<figuresNames[j].Data()<<"\n";
                if(j == 2)
                    fileTxt<<Form("\t\t\tM=%.3f (%.3f)\t#Gamma = %.3f (%.3f)",m0[k][i][j],m0Err[k][i][j],gamma0[k][i][j],gamma0Err[k][i][j])<<"\n";
                else
                    fileTxt<<Form("\t\t\tM=%.3f (%.3f)\t#Gamma = %.3f (%.3f)  \t#delta^{2} = %.2f (%.2f)",m0[k][i][j],m0Err[k][i][j],gamma0[k][i][j],gamma0Err[k][i][j],delta[k][i][j],deltaErr[k][i][j])<<"\n";
            }

        }
        fileTxt<<"\n";
    }


    fileTxt.close();

}
   

void setDeltaExp(TH1D* dataHist[][_N_PAIRS_]){ // 0 - norm, 1 - m_delta, 2 - gamma_delta, 3 - temp

    TF1* fun[_N_FIGURES_];
    fun[0] = new TF1(figuresNames[0].Data(),monitzSetDeltaExp,XMin,XMax,4);
    fun[1] = new TF1(figuresNames[1].Data(),manleySetDeltaExp,XMin,XMax,4);
    fun[2] = new TF1(figuresNames[2].Data(),fitFunctionM,XMin,XMax,3);


    for(int i = 0; i < _N_FIGURES_; i++) {
        fun[i]->SetLineColor(colors[i]);
        fun[i]->SetLineStyle(lines[i]);
    }

    for(int j = 0; j < _N_FIGURES_; j++) {
        if (j == 2) {
            fun[j]->SetParameters(1.0,1.232,0.117);
        }
        else {
            fun[j]->SetParameters(10,1.232,0.117,0.5);
            fun[j]->SetParLimits(3,0.01,0.5);
        }
        fun[j]->SetParLimits(1,1.1,1.4);
        //fun[j]->SetParLimits(0,1,10000);
        fun[j]->SetParLimits(2,0.0,1.0);

    }


    TFile* fileOut = new TFile("outputMNewTherm/outSetDeltaExp.root","RECREATE");

    TCanvas* can[_N_CASES_][_N_PAIRS_];

    Double_t m0[_N_CASES_][_N_PAIRS_][_N_FIGURES_];
    Double_t gamma0[_N_CASES_][_N_PAIRS_][_N_FIGURES_];
    Double_t m0Err[_N_CASES_][_N_PAIRS_][_N_FIGURES_];
    Double_t gamma0Err[_N_CASES_][_N_PAIRS_][_N_FIGURES_];
    Double_t temp[_N_CASES_][_N_PAIRS_][_N_FIGURES_];
    Double_t tempErr[_N_CASES_][_N_PAIRS_][_N_FIGURES_];


    for(int k = 0; k < _N_CASES_; k++) {
        cout<<"\t"<<casesNames[k].Data()<<endl;
        for(int i = 0; i < _N_PAIRS_; i++) {
            cout<<"\t\t"<<pairsTitles[i].Data()<<endl;
            cout<<"\t\t\t"<<canvasesNames[2].Data()<<endl;
            can[k][i] = new TCanvas(Form("%s%s%s",casesNames[k].Data(),pairsTitles[i].Data(),canvasesNames[2].Data()),Form("%s%s%s",casesNames[k].Data(),pairsTitles[i].Data(),canvasesNames[2].Data()),1000,1000);

            for(int j = 0; j < _N_FIGURES_; j++) {
                cout<<"\t\t\t\t"<<figuresNames[j].Data()<<endl;
                dataHist[k][i]->Fit(fun[j],"M+","M+",fitMin,fitMax);
                
                fun[j]->SetRange(XMin,XMax);
                dataHist[k][i]->GetListOfFunctions()->Add(fun[j]);

                m0[k][i][j] = fun[j]->GetParameter(1);
                gamma0[k][i][j] = fun[j]->GetParameter(2);
                m0Err[k][i][j] = fun[j]->GetParError(1);
                gamma0Err[k][i][j] = fun[j]->GetParError(2);

                if (j != 2) {
                    temp[k][i][j] = fun[j]->GetParameter(3);
                    tempErr[k][i][j] = fun[j]->GetParError(3);
                }

            }

            can[k][i]->cd();
            dataHist[k][i]->Draw();


            makePaveText(can[k][i], canvasesNamesNice[2].Data(),0.45, 0.85, 0.7, 0.9, 0.04);
            makePaveText(can[k][i], "Monitz", 0.6, 0.75, 0.8, 0.8, 0.04, colors[0]);
            makePaveText(can[k][i], "Manley", 0.6, 0.7, 0.8, 0.75, 0.04, colors[1]);
            makePaveText(can[k][i], "Breit-Wigner", 0.6, 0.65, 0.8, 0.7 , 0.04, kRed);
            makePaveText(can[k][i], pairsNames[i].Data(),0.7, 0.5, 0.8, 0.55, 0.05);

            fileOut->cd();
            dataHist[k][i]->Write();
            can[k][i]->SaveAs(Form("outputMNewTherm/%s%s%s.png",casesNames[k].Data(),pairsTitles[i].Data(),canvasesNames[2].Data()));


        }
    }

    fileOut->Save();
    fileOut->Close();


    ofstream fileTxt;
    fileTxt.open(Form("outputMNewTherm/%sFitResults.txt",canvasesNames[2].Data()));



    for(int k = 0; k < _N_CASES_; k++) {
        fileTxt<<casesNames[k].Data()<<"\n";
        for(int i = 0; i < _N_PAIRS_; i++) {
            fileTxt<<"\t"<<pairsTitles[i].Data()<<"\n";
            for(int j = 0; j < _N_FIGURES_; j++) {
                fileTxt<<"\t\t"<<figuresNames[j].Data()<<"\n";
                if(j == 2)
                    fileTxt<<Form("\t\t\tM=%.3f (%.3f)\t#Gamma = %.3f (%.3f)",m0[k][i][j],m0Err[k][i][j],gamma0[k][i][j],gamma0Err[k][i][j])<<"\n";
                else
                    fileTxt<<Form("\t\t\tM=%.3f (%.3f)\t#Gamma = %.3f (%.3f)  \t#delta^{2} = %.2f \tT = %.3f (%.3f) MeV",m0[k][i][j],m0Err[k][i][j],gamma0[k][i][j],gamma0Err[k][i][j],delta2Mon,temp[k][i][j]*1e3,tempErr[k][i][j]*1e3)<<"\n";
            }

        }
        fileTxt<<"\n";
    }

    fileTxt.close();



}

void parDeltaExp(TH1D* dataHist[][_N_PAIRS_]) { // 0 - norm, 1 - m_delta, 2 - gamma_delta, 3 - delta, 4 - temp
    TF1* fun[_N_FIGURES_];
    fun[0] = new TF1(figuresNames[0].Data(),monitzParDeltaExp,XMin,XMax,5);
    fun[1] = new TF1(figuresNames[1].Data(),manleyParDeltaExp,XMin,XMax,5);
    fun[2] = new TF1(figuresNames[2].Data(),fitFunctionM,XMin,XMax,3);

    for(int i = 0; i < _N_FIGURES_; i++) {
        fun[i]->SetLineColor(colors[i]);
        fun[i]->SetLineStyle(lines[i]);
    }

    for(int j = 0; j < _N_FIGURES_; j++) {
        if (j == 2) {
            fun[j]->SetParameters(1.0,1.232,0.117);
        }
        else {
            fun[j]->SetParameters(10.0,1.232,0.117,delta2Mon,0.5);
            fun[j]->SetParLimits(3,0.01,1.0);
            fun[j]->SetParLimits(4,0.01,0.5);
        }
        fun[j]->SetParLimits(1,1.1,1.4);
        //fun[j]->SetParLimits(0,0.1,1e6);
        fun[j]->SetParLimits(2,0.0,2.0);

    }


    TFile* fileOut = new TFile("outputMNewTherm/outParDeltaExp.root","RECREATE");

    TCanvas* can[_N_CASES_][_N_PAIRS_];

    Double_t m0[_N_CASES_][_N_PAIRS_][_N_FIGURES_];
    Double_t gamma0[_N_CASES_][_N_PAIRS_][_N_FIGURES_];
    Double_t m0Err[_N_CASES_][_N_PAIRS_][_N_FIGURES_];
    Double_t gamma0Err[_N_CASES_][_N_PAIRS_][_N_FIGURES_];
    Double_t temp[_N_CASES_][_N_PAIRS_][_N_FIGURES_];
    Double_t tempErr[_N_CASES_][_N_PAIRS_][_N_FIGURES_];
    Double_t delta[_N_CASES_][_N_PAIRS_][_N_FIGURES_];
    Double_t deltaErr[_N_CASES_][_N_PAIRS_][_N_FIGURES_];


    for(int k = 0; k < _N_CASES_; k++) {
        cout<<"\t"<<casesNames[k].Data()<<endl;
        for(int i = 0; i < _N_PAIRS_; i++) {
            cout<<"\t\t"<<pairsTitles[i].Data()<<endl;
            cout<<"\t\t\t"<<canvasesNames[3].Data()<<endl;
            can[k][i] = new TCanvas(Form("%s%s%s",casesNames[k].Data(),pairsTitles[i].Data(),canvasesNames[3].Data()),Form("%s%s%s",casesNames[k].Data(),pairsTitles[k].Data(),canvasesNames[3].Data()),1000,1000);

            for(int j = 0; j < _N_FIGURES_; j++) {
                cout<<"\t\t\t\t"<<figuresNames[j].Data()<<endl;
                dataHist[k][i]->Fit(fun[j],"M+","M+",fitMin,fitMax);
                
                fun[j]->SetRange(XMin,XMax);
                dataHist[k][i]->GetListOfFunctions()->Add(fun[j]);

                m0[k][i][j] = fun[j]->GetParameter(1);
                gamma0[k][i][j] = fun[j]->GetParameter(2);
                m0Err[k][i][j] = fun[j]->GetParError(1);
                gamma0Err[k][i][j] = fun[j]->GetParError(2);

                if (j != 2) {
                    temp[k][i][j] = fun[j]->GetParameter(4);
                    tempErr[k][i][j] = fun[j]->GetParError(4);
                    delta[k][i][j] = fun[j]->GetParameter(3);
                    deltaErr[k][i][j] = fun[j]->GetParError(3);
                }

            }

            can[k][i]->cd();
            dataHist[k][i]->Draw();


            makePaveText(can[k][i], canvasesNamesNice[3].Data(),0.45, 0.85, 0.7, 0.9, 0.04);
            makePaveText(can[k][i], "Monitz", 0.6, 0.75, 0.8, 0.8, 0.04, colors[0]);
            makePaveText(can[k][i], "Manley", 0.6, 0.7, 0.8, 0.75, 0.04, colors[1]);
            makePaveText(can[k][i], "Breit-Wigner", 0.6, 0.65, 0.8, 0.7 , 0.04, kRed);
            makePaveText(can[k][i], pairsNames[i].Data(),0.7, 0.5, 0.8, 0.55, 0.05);

            fileOut->cd();
            dataHist[k][i]->Write();
            can[k][i]->SaveAs(Form("outputMNewTherm/%s%s%s.png",casesNames[k].Data(),pairsTitles[i].Data(),canvasesNames[3].Data()));


        }
    }


    fileOut->Save();
    fileOut->Close();

    ofstream fileTxt;
    fileTxt.open(Form("outputMNewTherm/%sFitResults.txt",canvasesNames[3].Data()));



    for(int k = 0; k < _N_CASES_; k++) {
        fileTxt<<casesNames[k].Data()<<"\n";
        for(int i = 0; i < _N_PAIRS_; i++) {
            fileTxt<<"\t"<<pairsTitles[i].Data()<<"\n";
            for(int j = 0; j < _N_FIGURES_; j++) {
                fileTxt<<"\t\t"<<figuresNames[j].Data()<<"\n";
                if(j == 2)
                    fileTxt<<Form("\t\t\tM=%.3f (%.3f)\t#Gamma = %.3f (%.3f)",m0[k][i][j],m0Err[k][i][j],gamma0[k][i][j],gamma0Err[k][i][j])<<"\n";
                else
                    fileTxt<<Form("\t\t\tM=%.3f (%.3f)\t#Gamma = %.3f (%.3f)  \t#delta^{2} = %.2f \tT = %.3f (%.3f) MeV",m0[k][i][j],m0Err[k][i][j],gamma0[k][i][j],gamma0Err[k][i][j],delta[k][i][j],temp[k][i][j]*1e3,tempErr[k][i][j]*1e3)<<"\n";
            }

        }
        fileTxt<<"\n";
    }

    fileTxt.close();

}


void compCanvases() {

    TH1D* hist[_N_CASES_][_N_PAIRS_][_N_CANVASES_];

    TFile* fileIn[_N_CANVASES_];
    for(int i = 0; i < _N_CANVASES_; i++) {
        fileIn[i] = new TFile(Form("outputMNewTherm/out%s.root",canvasesNames[i].Data()));
        for(int k = 0; k < _N_CASES_; k++) {
            for(int j = 0; j < _N_PAIRS_; j++) {
                hist[k][j][i] = (TH1D*)fileIn[i]->Get(Form("%s%sMFitHist",casesNames[k].Data(),pairsTitles[j].Data()));
                //hist[k][j][i]->SetTitle(canvasesNamesNice[i].Data());
                hist[k][j][i]->SetMarkerSize(0.9);
                hist[k][j][i]->GetXaxis()->SetTitleSize(0.05);
                hist[k][j][i]->GetYaxis()->SetTitleSize(0.05);
                hist[k][j][i]->GetXaxis()->SetLabelSize(0.05);
                hist[k][j][i]->GetYaxis()->SetLabelSize(0.05);

             }
        }
    }

    TCanvas* canOut[_N_CASES_][_N_PAIRS_];
    TFile* fileOut = new TFile("outputMNewTherm/outComp.root","RECREATE");
    for(int k = 0; k < _N_CASES_; k++) {
        for(int i = 0; i < _N_PAIRS_; i++) {
            canOut[k][i] = new TCanvas(Form("%s%sComp",casesNames[k].Data(),pairsTitles[i].Data()),Form("%s%sComp",casesNames[k].Data(),pairsTitles[i].Data()),1000,1000);
            canOut[k][i]->Divide(2,2);
            for(int j = 0; j < _N_CANVASES_; j++) {
                canOut[k][i]->cd(j+1);
                hist[k][i][j]->Draw();
                gStyle->SetPadLeftMargin(0.13);
                gStyle->SetPadBottomMargin(0.11);
                if (j == 1) {
                    makePaveText(canOut[k][i], canvasesNamesNice[j].Data(),0.69, 0.9, 0.9, 0.95, 0.025);
                    makePaveText(canOut[k][i], "Monitz", 0.8, 0.85, 0.9, 0.88, 0.025, colors[0]);
                    makePaveText(canOut[k][i], "Manley", 0.8, 0.82, 0.9, 0.85, 0.025, colors[1]);
                    makePaveText(canOut[k][i], "Breit-Wigner", 0.8, 0.79, 0.9, 0.82 , 0.025, kRed);
                    makePaveText(canOut[k][i], casesNames[k].Data(),0.85, 0.73, 0.9, 0.77, 0.025);
                    makePaveText(canOut[k][i], pairsNames[i].Data(),0.85, 0.69, 0.9, 0.73, 0.025);
                }
                else if(j == 0) 
                    makePaveText(canOut[k][i], canvasesNamesNice[j].Data(),0.25, 0.9, 0.5, 0.95, 0.025);
                else if(j == 3)
                    makePaveText(canOut[k][i], canvasesNamesNice[j].Data(),0.68, 0.4, 0.9, 0.45, 0.025);
                else
                    makePaveText(canOut[k][i], canvasesNamesNice[j].Data(),0.25, 0.4, 0.5, 0.45, 0.025);
            }

            fileOut->cd();
            canOut[k][i]->Write();
            canOut[k][i]->SaveAs(Form("outputMNewTherm/%s%sComp.png",casesNames[k].Data(),pairsTitles[i].Data()));

//canOut[k][i]->SaveAs(Form("outputMNewTherm/%s%sComp.eps",casesNames[k].Data(),pairsTitles[i].Data()));

        }

    }


        fileOut->Close();
        fileOut->Save();





}



void compMFitNewTherm() {
    
    //gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(20000);

    TH1D *dataHist[_N_CASES_][_N_PAIRS_];

    TFile* file[_N_CASES_];
    TString fileNames[_N_CASES_]= {"outputBasic/outCaseABasic.root","outputBasic/outCaseBBasic.root","outputBasic/outCaseCBasic.root"};
    Float_t events[_N_CASES_] = {9.97*10e6, 10e6, 10e6};
    Float_t dX = (XMax-XMin)/XBins;
    Int_t rebin = 10; 

    for(int i = 0; i < _N_CASES_; i++){
        file[i] = new TFile(fileNames[i].Data());
        for(int j = 0; j < _N_PAIRS_; j++){

                dataHist[i][j] = (TH1D*)file[i]->Get(Form("%sMHist",pairsTitles[j].Data()))->Clone(Form("%s%sMFitHist",casesNames[i].Data(),pairsTitles[j].Data()));       
                dataHist[i][j]->Rebin(rebin);
                dataHist[i][j]->Scale(1.0/events[i]/dX/rebin);
                dataHist[i][j]->GetXaxis()->SetRangeUser(1.05,XMax);
                dataHist[i][j]->SetMarkerSize(1.0);
                dataHist[i][j]->SetMarkerStyle(8);
            
        }
    }


    gStyle->SetOptStat(0);
    //setDeltaNoExp(dataHist);
   //parDeltaNoExp(dataHist);
    //setDeltaExp(dataHist);
   //parDeltaExp(dataHist);

    compCanvases();





}
