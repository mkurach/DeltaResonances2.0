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
#define _N_MRANGES_ 24
#define _N_FIGURES_ 2

int colors[]={1, 600, 629, 414, 802, 880, 819, 922,433,618}; // kBlack, kBlue, kRed, kGreen-2, kOrange+2, kGray, kViolet, kSpring
int markers[]= {20,21,24,25,28,34,47,43}; 



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

void multiplePlot(TH1D *histTab[], TCanvas* can, TString* entry, size_t n, bool logy = false, bool scaling = false, TString XLabel = "") {
    if(XLabel!="")
        histTab[0]->GetXaxis()->SetTitle(XLabel);

    Double_t max = 0;
    Double_t min = 0;

    for (int i = 0; i < n; i++) {
        if(scaling) {
            histTab[i]->Scale(pow(10,i));
            histTab[i]->GetFunction("fun")->SetParameter(1,histTab[i]->GetFunction("fun")->GetParameter(1)*pow(10,i));
            histTab[i]->GetFunction("fun")->SetLineColor(kBlack);

            //histTab[i]->Sumw2(kFALSE); //recalculate errors
            //histTab[i]->Sumw2();

            //entry[i] = TString::Format("(%i) %s ",i,entry[i].Data());
            //entry[i] += TString::Format(" 10^{%i}",i); //KRYPTYDA
        }
        if(histTab[i]->GetMaximum()>max)
            max = histTab[i]->GetMaximum();
        if(histTab[i]->GetMinimum()<min)
            min = histTab[i]->GetMinimum();
    }

    if (scaling) {
        histTab[0]->SetMaximum(1.15*10*max);
    }
    else {
        //histTab[0]->SetMaximum(2*max);
        histTab[0]->GetYaxis()->SetRangeUser(min,max*1.15);
    }
    histTab[0]->SetTitle("");

    TLegend* legend	= new TLegend(0.7,0.5,0.9,0.9, "",	"NDC");
    legend->SetFillColor(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.03);
    legend->SetLineWidth(0);
    legend->SetName("legend");

    can->cd();
    gPad->SetLeftMargin(0.3);
    for (int i = 0; i < n; i++) {
        
        /*histTab[i]->GetXaxis()->SetTitleSize(0.04);
        histTab[i]->GetYaxis()->SetTitleSize(0.04);
        histTab[i]->GetXaxis()->SetTitleFont(42);
        histTab[i]->GetYaxis()->SetTitleFont(42);
        histTab[i]->GetYaxis()->SetTitleOffset(1.2);

        histTab[i]->GetXaxis()->SetLabelFont(42);
        histTab[i]->GetYaxis()->SetLabelFont(42);
        histTab[i]->GetXaxis()->SetLabelSize(0.035);
        histTab[i]->GetYaxis()->SetLabelSize(0.035);*/

        //histTab[i]->SetMarkerStyle(markers[i]);
        histTab[i]->SetMarkerStyle(markers[0]);
        histTab[i]->SetMarkerSize(0.8);
        histTab[i]->SetMarkerColor(colors[i]);
        histTab[i]->SetLineColor(colors[i]);

        if(scaling && i!=0)
                legend->AddEntry(histTab[i],TString::Format("%s (x10^{%i})",entry[i].Data(),i));
        else
            legend->AddEntry(histTab[i],entry[i]);

        if(logy)
            gPad->SetLogy();
        if(scaling)
                    
        if (i == 0)
            histTab[i]->Draw();
        else {
            histTab[i]->Draw("same");
        }

        
    }
    legend->Draw("same");


}


Double_t getXRange(int i) {
    Double_t x;
    switch (i) {
        case 0:  x = (1.078+1.09)/2;   break;
        case 1:  x = (1.09+1.1)/2;     break;
        case 2:  x = (1.1+1.11)/2;    break;
        case 3:  x = (1.11+1.12)/2;    break;
        case 4:  x = (1.12+1.13)/2;     break;
        case 5:  x = (1.13+1.14)/2;     break;
        case 6:  x = (1.14+1.15)/2;     break;
        case 7:  x = (1.15+1.16)/2;     break;
        case 8:  x = (1.16+1.17)/2;     break;
        case 9:  x = (1.17+1.18)/2;     break;
        case 10:  x = (1.18+1.19)/2;    break;
        case 11:  x = (1.19+1.2)/2;     break;
        case 12:  x = (1.2+1.22)/2;     break;
        case 13:  x = (1.22+1.24)/2;     break;
        case 14:  x = (1.24+1.28)/2;    break;
        case 15:  x = (1.28+1.32)/2;     break;
        case 16:  x = (1.32+1.36)/2;    break;
        case 17:  x = (1.36+1.4)/2;     break;
        case 18:  x = (1.4+1.44)/2;     break;
        case 19:  x = (1.44+1.48)/2;     break;
        case 20:  x = (1.48+1.52)/2;     break;
        case 21:  x = (1.52+1.56)/2;     break;
        case 22:  x = (1.56+1.6)/2;     break;
        case 23:  x = (1.6+1.65)/2;     break;
        case 24:  x = (1.65+1.7)/2;     break;
    }

    return x;

}

Double_t fitFunctionMtMMidY(Double_t *x, Double_t *par) {
    if(x[0] >= 1.0 && x[0] <= 3.0)
        return par[1]*TMath::Exp(-x[0]/par[0]);
    else    
        return 0;
}

Double_t fitFunctionMtMFullY(Double_t *x, Double_t *par) {
    if(x[0] >= 1.0 && x[0] <= 3.0)
        return par[1]*TMath::Exp(-x[0]/par[0]);
    else    
        return 0;
}

Double_t fitFunctionTeffM(Double_t *x, Double_t *par) {
    if(x[0] >= 0 && x[0] <= 1.7)
        return TMath::Power(par[0],2)/2*x[0]+par[1];
    else    
        return 0;
}

void createMtMFit(int caseInt){
    
    TString Case;
    switch(caseInt) {
        case 0:
            Case = "CaseA";
            break;
        case 1:
            Case = "CaseB";
            break;
        case 2:
            Case = "CaseC";
            break;
    }
        
    //gROOT -> SetBatch(kTRUE); 

    TH1D* histMtM[_N_PAIRS_][_N_FIGURES_][_N_MRANGES_];
    TString pairsTitles[_N_PAIRS_] = {"PiPlusP","PiMinusP"};
    TString pairsNames[_N_PAIRS_] = {"#pi^{+}p","#pi^{-}p"};
    //TString MRanges[_N_MRANGES_] = {"1 #leq M < 1.1","1.1 #leq M < 1.15","1.15 #leq M < 1.2","1.2 #leq M < 1.25","1.25 #leq M < 1.3","1.3 #leq M < 1.35","1.35 #leq M < 1.4","1.4 #leq M < 1.45","1.45 #leq M < 1.5","1.5 #leq M < 1.6"};
    TString MRanges[_N_MRANGES_] = {"1.078 #leq M < 1.09","1.09 #leq M < 1.1","1.1 #leq M < 1.11","1.11 #leq M < 1.12","1.12 #leq M < 1.13",
                                "1.13 #leq M < 1.14","1.14 #leq M < 1.15","1.15 #leq M < 1.16","1.16 #leq M < 1.17","1.17 #leq M < 1.18",
                                "1.18 #leq M < 1.19","1.19 #leq M < 1.2","1.2 #leq M < 1.22","1.22 #leq M < 1.24","1.24 #leq M < 1.28",
                                "1.28 #leq M < 1.32","1.32 #leq M < 1.36","1.36 #leq M < 1.4","1.4 #leq M < 1.44","1.44 #leq M < 1.48",
                                "1.48 #leq M < 1.52","1.52 #leq M < 1.56","1.56 #leq M < 1.6","1.6 #leq M < 1.65"};//,"1.65 #leq M < 1.7"};
    TString histNames[_N_FIGURES_] = {"MidY","FullY"};
    
    Float_t events[_N_CASES_] = {9.97*10e6, 10e6, 10e6};
    Int_t XBins = 1000;
    Float_t XMin = 1;
	Float_t XMax = 3;
	Float_t dX = (XMax-XMin)/XBins;
    Float_t scale = 1.0/(events[caseInt]*dX);
    Int_t rebin = 10; 

    //TCanvas * can[_N_PAIRS_][_N_MRANGES_];
    TFile* file = new TFile("outputMt/out"+Case+"Mt.root");
    for(int j = 0; j < _N_PAIRS_; j++){
        for(int k = 0; k < _N_FIGURES_; k++){ 
            for(int l = 0; l < _N_MRANGES_; l++) {
                //can[j][k] = new TCanvas(Form("%s%sPtMFitting%i",Case.Data(),pairsTitles[j].Data(),k),Form("%s%sPtMFitting%i",Case.Data(),pairsTitles[j].Data(),k),1000,1000);
                histMtM[j][k][l] = (TH1D*)file->Get(Form("%sMtM%sHist%i",pairsTitles[j].Data(),histNames[k].Data(),l));

                //histMtM[j][k][l]->Sumw2();
                //histMtM[j][k]->Rebin(rebin);
                //histMtM[j][k]->Scale(scale/rebin);
                histMtM[j][k][l]->Scale(scale);
            }
        }
    }

    //FUNCTION
    Int_t nparams = 2;
    TF1 *fun[_N_FIGURES_];


    fun[0] = new TF1("funMidY",fitFunctionMtMMidY,1.0,3.0,nparams);
    fun[1]= new TF1("funFullY",fitFunctionMtMFullY,1.0,3.0,nparams);
    //0 - Teff,  1 - normalizacja

    Double_t teffParameter[_N_CASES_] = {0.15,0.1,0.1};
    Double_t normParameter[_N_CASES_] = {2e3,5e5,2e5};

    fun[0]->SetParameter(0,teffParameter[caseInt]);
    fun[0]->SetParameter(1,normParameter[caseInt]);

    fun[1]->SetParameter(0,teffParameter[caseInt]);
    fun[1]->SetParameter(1,normParameter[caseInt]);

    //FITTING

    gROOT -> SetBatch(kTRUE);

    Double_t fitMinOffset[_N_CASES_] = {0.3,0.15,0.15};
    Double_t fitMaxOffset[_N_CASES_] = {0.9,0.8,0.8};

    Double_t minRange = 1;
    Double_t maxRange = 3;

    TF1* funCop;

    Double_t tEff[_N_PAIRS_][_N_FIGURES_][_N_MRANGES_];
    Double_t tEffEr[_N_PAIRS_][_N_FIGURES_][_N_MRANGES_];
    Double_t mXValues[_N_PAIRS_][_N_FIGURES_][_N_MRANGES_];
    Double_t XErr[_N_PAIRS_][_N_FIGURES_][_N_MRANGES_];

    TFile* fileOut[_N_FIGURES_];
    fileOut[0] = new TFile(Form("outputMt/outFitMtMMid%s.root",Case.Data()),"RECREATE");
    fileOut[1] = new TFile(Form("outputMt/outFitMtMFull%s.root",Case.Data()),"RECREATE");


    for(int j = 0; j < _N_PAIRS_; j++) {
        for(int k = 0; k < _N_FIGURES_; k++) { 
            for(int l = 0; l < _N_MRANGES_ - 7; l++) {
                //cout<<Case<<"  "<<pairsTitles[j].Data()<<" numer range'u:"<<k<<endl;
                //histMtM[j][k][l]->Fit(Form("fun%s",histNames[l].Data()),"Q","Q",fitMin[caseInt],fitMax[caseInt]);
                histMtM[j][k][l]->Fit(Form("fun%s",histNames[k].Data()),"Q","Q",getXRange(l)+fitMinOffset[caseInt],getXRange(l)+fitMaxOffset[caseInt]);
                funCop = (TF1*)fun[k]->Clone("funCop");
                funCop->SetRange(minRange,maxRange);
                funCop->SetLineStyle(2);
                histMtM[j][k][l]->GetListOfFunctions()->Add(funCop);
                fileOut[k]->cd();
                histMtM[j][k][l]->Write();

                tEff[j][k][l] = fun[k]->GetParameter(0);
                tEffEr[j][k][l] = fun[k]->GetParError(0);
                mXValues[j][k][l] = getXRange(l);
                XErr[j][k][l] = 0;
            }

        }
    }

    fileOut[0]->Save();
    fileOut[0]->Close();

    fileOut[1]->Save();
    fileOut[1]->Close();


    gROOT -> SetBatch(kFALSE);

    //TEFF PRINTING

    /*for(int j = 0; j < _N_PAIRS_; j++) {
        for(int k = 0; k < _N_MRANGES_; k++) {
            cout<<Case<<"  "<<pairsTitles[j].Data()<<" numer range'u:"<<k<<endl;
            cout<<"\tTeff: "<<tEff[j][k]<<"\n"<<endl;
        }
        cout<<"******"<<endl;
    }*/
    
    //TEFF(M) PLOTS

    //TCanvas* can[_N_PAIRS_][_N_FIGURES_];
    TGraphErrors *grTeff[_N_PAIRS_][_N_FIGURES_];
    for(int i = 0; i < _N_PAIRS_; i++) {
        for(int j = 0; j < _N_FIGURES_; j++) {
            grTeff[i][j] = new TGraphErrors(_N_MRANGES_ -7,mXValues[i][j],tEff[i][j],XErr[i][j],tEffEr[i][j]);
            //can[i][j] = new TCanvas(Form("%i%i",i,j),Form("%i%i",i,j),1000,1000);
            grTeff[i][j]->GetYaxis()->SetTitle("T_{eff} (GeV)");
            grTeff[i][j]->GetXaxis()->SetTitle(Form("M_{%s} (GeV/c^{2})",pairsNames[i].Data()));

            //for(int j = 0; j < 3; j ++)
                //grTeff[i]->SetPoint(25+j,addPointsX[j],addPointsY[caseInt][j]);
            grTeff[i][j]->GetXaxis()->SetLimits(0,1.8);
            grTeff[i][j]->GetYaxis()->SetRangeUser(-0.05,0.3);
            grTeff[i][j]->GetYaxis()->SetTitleOffset(1.3);
            grTeff[i][j]->GetXaxis()->SetTitleOffset(1.3);
            grTeff[i][j]->SetTitle(Form("%s%s%sTeffM",Case.Data(),pairsTitles[i].Data(),histNames[j].Data()));
            grTeff[i][j]->SetMarkerStyle(20);
            grTeff[i][j]->SetMarkerSize(1);
            grTeff[i][j]->SetMarkerColor(kRed+2);
           // can[i][j]->cd();
            //grTeff[i][j]->Draw("ap");
            //fileOut->cd();
            //can[i][j]->SaveAs(Form("outputMt/TeffM%s%s%s.png",Case.Data(),pairsTitles[i].Data(),histNames[j].Data()));
        }
    }

    /*fileOut->Save();
	fileOut->Close();*/

    //PLOTTING nie działa na ten moment

    /*TCanvas* canAll[_N_PAIRS_];
    
    for(int i = 0; i < _N_PAIRS_; i ++) {
        canAll[i] = new TCanvas(Form("%s%sPtMFitted", Case.Data(),pairsTitles[i].Data()),Form("%s%sPtMFitted", Case.Data(),pairsTitles[i].Data()),1000,1000);
        multiplePlot(histMtM[i],canAll[i],MRanges,10,true,true,"M (GeV/c^{2})");
        makePaveText(canAll[i],Case.Data(),0.73,0.4,0.99,0.5,0.05);
        makePaveText(canAll[i],pairsNames[i].Data(),0.73,0.35,0.99,0.4,0.05);
        //canAll[i]->SaveAs(Form("outputFit/%s%sPtMFitted.png",Case.Data(),pairsTitles[i].Data()));
    }*/

    //TCanvas* canTeff[_N_PAIRS_];
    //canTeff[i] = new TCanvas(Form("%s%sTeff",Case.Data(),pairsTitles[i].Data()),Form("%s%sTeff",Case.Data(),pairsTitles[i].Data()),1000,1000);
    //singlePlot(histTeff[i],canTeff[i],"Rapidity", "T_{eff} (GeV)");

    
    
    //FUNCTION

    nparams = 2;
    TF1 *fun2 = new TF1("fun2",fitFunctionTeffM,0,1.7,nparams);
    //0 - beta, 1 - T_
    TF1* fun2Cop;

    //fun2->SetParLimits(0,0,1);
    fun2->SetParameter(0,0.5);

    //fun2->SetParLimits(1,-0.1,0.2);
    fun2->SetParameter(1,0.1);

    
    

    //ADDING POINTS FROM PARTICLES

    Double_t addPointsX[3] = {0.9395653,0.4936770,0.1395699};
    /*Double_t addPointsY[_N_CASES_][3] = {{0.145643,0.105242,0.0762062},
                                        {0.11814,0.0991819,0.0716104},
                                        {0.113079,0.0952766,0.0752646}};*/

    /*Double_t addPointsY[_N_CASES_][3][_N_FIGURES_] = {{{0.129938,0.125233},{0.101965,0.100151},{0.0679871,0.0688109}},
                                        {{0.112855,0.111528},{0.0940142,0.0950485},{0.0641781,0.065325}},
                                        {{0.108625,0.107424},{0.088614,0.0912045},{0.0667562,0.0683127}}};*/

    Double_t addPointsY[_N_CASES_][_N_FIGURES_][3] = {{{0.129938,0.101965,0.0679871},{0.125233,0.100151,0.0688109}},
                                                        {{0.112855,0.0940142,0.0641781},{0.111528,0.0950485,0.065325}},
                                                        {{0.108625,0.088614,0.0667562},{0.107424,0.0912045,0.0683127}}};


    TGraph* grPart[_N_FIGURES_];

    Double_t betaPar[_N_FIGURES_];
    Double_t tPar[_N_FIGURES_];
    Double_t betaParEr[_N_FIGURES_];
    Double_t tParEr[_N_FIGURES_];

    for(int j = 0; j < _N_FIGURES_; j++) {
        grPart[j]= new TGraph(3,addPointsX,addPointsY[caseInt][j]);
        grPart[j]->SetMarkerStyle(20);
        grPart[j]->SetMarkerSize(1);
        grPart[j]->SetMarkerColor(kBlue+2);
        grPart[j]->GetXaxis()->SetLimits(0,1.7);
        grPart[j]->GetYaxis()->SetRangeUser(-0.05,0.3);
        grPart[j]->GetYaxis()->SetTitle("T_{eff} (GeV)");
        grPart[j]->GetXaxis()->SetTitle("M (GeV/c^{2})");
        grPart[j]->SetTitle("");
        grPart[j]->GetYaxis()->SetTitleOffset(1.6);
        grPart[j]->GetXaxis()->SetTitleOffset(1.3);

        fun2->SetLineColor(kBlue);
        grPart[j]->Fit("fun2","","",0.1,1);
        betaPar[j] = fun2->GetParameter(0);
        tPar[j] = fun2->GetParameter(1);
        betaParEr[j] = fun2->GetParError(0);
        tParEr[j] = fun2->GetParError(1);

        fun2Cop = (TF1*)fun2->Clone("fun2Cop");
        fun2Cop->SetRange(0,1.6);
        fun2Cop->SetLineStyle(2);
        fun2Cop->SetLineColor(kBlue);
        grPart[j]->GetListOfFunctions()->Add(fun2Cop);
    }


    




    //TEFF M FIT


    Double_t beta[_N_PAIRS_][_N_FIGURES_];
    Double_t t[_N_PAIRS_][_N_FIGURES_];
    Double_t betaEr[_N_PAIRS_][_N_FIGURES_];
    Double_t tEr[_N_PAIRS_][_N_FIGURES_];

    Double_t betaTh[_N_CASES_] = {0.525672, 0.469742, 0.50362};
    Double_t tTh[_N_CASES_] = {49.6, 70.3, 63.1}; //(MeV)



    TCanvas* can[_N_PAIRS_][_N_FIGURES_];
    //TCanvas* canParticles[_N_PAIRS_];
   // TFile *fileOut2 = new TFile(Form("/u/mkurach/figures_with_data/moje/ladne/teffM%s.root",Case.Data()),"RECREATE");
    for(int j = 0; j < _N_PAIRS_; j++) {
        for(int k = 0; k < _N_FIGURES_; k++) {
        
            can[j][k] = new TCanvas(Form("%sTeffMFit",pairsTitles[j].Data()),Form("%sTeffMFit",pairsTitles[j].Data()),1015,1000);
            //canParticles[j] = new TCanvas(Form("%sTeffMFitParticles",pairsTitles[j].Data()),Form("%sTeffMFitParticles",pairsTitles[j].Data()),715,700);
            setBasicStyle();
            setCanvas(can[j][k]);
            //setCanvas(canParticles[j]);

            //cout<<pairsTitles[j].Data()<<endl;
            fun2->SetLineColor(kRed);
            grTeff[j][k]->Fit("fun2","Q","Q",1.07,1.3);
            beta[j][k] = fun2->GetParameter(0);
            t[j][k] = fun2->GetParameter(1);
            betaEr[j][k] = fun2->GetParError(0);
            tEr[j][k] = fun2->GetParError(1);

            fun2Cop = (TF1*)fun2->Clone(Form("fun2Cop%i",j));
            fun2Cop->SetRange(0,1.67);
            fun2Cop->SetLineStyle(2);
            fun2Cop->SetLineColor(kRed);
            grTeff[j][k]->GetListOfFunctions()->Add(fun2Cop);


            can[j][k]->cd();
            can[j][k]->SetLeftMargin(0.13);
            grTeff[j][k]->SetTitle("");
            grTeff[j][k]->GetXaxis()->SetRangeUser(0,1.7);
            //grTeff[j][k]->Draw("ap");
           // makePaveText(can[j][k],Case.Data(),0.6,0.6,0.99,0.99,0.04);
           // makePaveText(can[j][k],Form("<#beta> = %.3f",beta[j][k]),0.15,0.6,0.3,0.7,0.04,kRed+3);
            //makePaveText(can[j][k],Form("T = %.3f (MeV)",t[j][k]*1e3),0.15,0.55,0.3,0.65,0.04,kRed+3);
            //can[j][k]->SaveAs(Form("outputMt/%s%s%sTeffM.png",Case.Data(),pairsTitles[j].Data(),histNames[k].Data()));
            //can[j]->SaveAs(Form("/u/mkurach/figures_with_data/moje/ladne/%s%steffM.pdf",Case.Data(),pairsTitles[j].Data()));
            //can[j]->Write();

            //canParticles[j]->cd();
            
            grPart[k]->Draw("ap ");
            grTeff[j][k]->Draw("p same");
            
            makePaveText(can[j][k],Case.Data(),0.6,0.6,0.99,0.99,0.04);
            makePaveText(can[j][k],Form("<#beta> = %.2f (%.2f)",beta[j][k],betaEr[j][k]),0.55,0.3,0.7,0.4,0.035,kRed+3);
            makePaveText(can[j][k],Form("T = %.2f (%.2f) MeV",t[j][k]*1e3,tEr[j][k]*1e3),0.55,0.34,0.7,0.44,0.035,kRed+3);

            //makePaveText(can[j][k],Case.Data(),0.6,0.6,0.99,0.99,0.04);0.15,0.6,0.3,0.7,0.035
            makePaveText(can[j][k],Form(" <#beta> = %.2f (%.2f)",betaPar[k],betaParEr[k]),0.15,0.55,0.3,0.65,0.035,kBlue+3);
            makePaveText(can[j][k],Form("T = %.2f (%.2f) MeV",tPar[k]*1e3,tParEr[k]*1e3),0.15,0.59,0.3,0.69,0.035,kBlue+3);

            can[j][k]->SaveAs(Form("outputMt/%s%s%sTeffM.png",Case.Data(),pairsTitles[j].Data(),histNames[k].Data()));
            //canParticles[j]->SaveAs(Form("/u/mkurach/figures_with_data/moje/ladne/%s%steffMParticles.pdf",Case.Data(),pairsTitles[j].Data()));
            //canParticles[j]->Write();
        }

    }

    //fileOut2->Close();
    //fileOut2->Save();
    




    
    

    //PRINTING
    /*cout<<Case<<endl;
    for(int j = 0; j < _N_PAIRS_; j++) {
        cout<<"\t"<<pairsTitles[j].Data()<<":\n"<<endl;
        cout<<"\t\tBeta: "<<beta[j];
        cout<<"\tT: "<<t[j]*1e3<<" (MeV)\n"<<endl;
    }
    cout<<"\tTheoretical:\n"<<endl;
    cout<<"\t\tBeta: "<<betaTh[caseInt];
    cout<<"\tT: "<<tTh[caseInt]<<" (MeV)\n"<<endl;

    cout<<"\tFrom 3 particles:\n"<<endl;
    cout<<"\t\tBeta: "<<betaPar;
    cout<<"\tT: "<<tPar*1e3<<" (MeV)\n"<<endl;

    //cout<<"chi2: "<<fun2->GetChisquare()<<endl;
    //cout<<"ndf: "<<fun2->GetNDF()<<endl;*/

}

void compMtMFit() {
    //TGraphErrors grTeff[_N_CASES_][_N_PAIRS_];
   //for(int i = 0; i < _N_CASES_; i++)
    //  createMtMFit(i);
    //createMtMFit(0);    
    //createMtMFit(1);
    createMtMFit(2);
}