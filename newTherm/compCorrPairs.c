#define _N_CASES_ 3
#define _N_VERSIONS_ 4

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


void doublePlot(TH1D *gr1, TH1D *gr2, TCanvas* can, TString entry1, TString entry2, TString XLabel = "") {
gStyle->SetOptStat(0);
    if(XLabel!="")
        gr1->GetXaxis()->SetTitle(XLabel);

    Double_t max = 0;
    if(gr1->GetMaximum()>max)
        max = gr1->GetMaximum();
    if(gr2->GetMaximum()>max)
        max = gr2->GetMaximum();
    //gr2->SetMinimum(-0.5);
    //gr2->SetMaximum(max*1.1);
    //gr1->GetYaxis()->SetRangeUser(-0.5,max*1.1);
    gr2->Scale(gr1->Integral(gr1->FindBin(1.15),gr1->FindBin(2.0))/gr2->Integral(gr2->FindBin(1.15),gr2->FindBin(2.0)));
//H1D[i][1]->Integral(minRange,maxRange) / H1D[i][0]->Integral(minRange,maxRange
    gr1->GetXaxis()->SetTitleSize(0.04);
    gr1->GetYaxis()->SetTitleSize(0.04);
    gr1->GetXaxis()->SetTitleFont(42);
    gr1->GetYaxis()->SetTitleFont(42);
    gr1->GetYaxis()->SetTitleOffset(1.2);

    gr1->GetXaxis()->SetLabelFont(42);
    gr1->GetYaxis()->SetLabelFont(42);
    gr1->GetXaxis()->SetLabelSize(0.035);
    gr1->GetYaxis()->SetLabelSize(0.035);

    //gr2->GetXaxis()->SetRangeUser(1,1.7);
    //gr->GetXaxis()->SetRangeUser(1,1.7);
    //gr->GetYaxis()->SetLimits(-0.5, max*1.1);
    //gr2->SetMinimum(-0.5);
    

   //gr1->SetErrorY(0);

    //gr1->SetMarkerStyle(20);
    //gr1->SetMarkerSize(0.8);
    //gr1->SetMarkerColor(kBlue);
    //gr1->SetDrawOption("H");
    gr1->SetLineColor(kBlue+1);
    gr1->SetLineStyle(1);
    gr1->SetLineWidth(3);


    gr2->GetXaxis()->SetTitleSize(0.04);
    gr2->GetYaxis()->SetTitleSize(0.04);
    gr2->GetXaxis()->SetTitleFont(42);
    gr2->GetYaxis()->SetTitleFont(42);

    gr2 -> GetYaxis() -> SetTitleOffset(1.2);
    gr2->GetXaxis()->SetLabelFont(42);
    gr2->GetYaxis()->SetLabelFont(42);
    gr2->GetXaxis()->SetLabelSize(0.035);
    gr2->GetYaxis()->SetLabelSize(0.035);
    gr2 -> SetMarkerStyle(21);
    gr2 -> SetMarkerSize(0);
    gr2 -> SetMarkerColor(kRed+1);
    gr2 -> SetLineColor(kRed+1);
    //gr2 ->SetLineStyle(2);
    gr2 ->SetLineWidth(2);

    TLegend* legend	= new TLegend(0.6,0.72,0.9,0.82, "",	"NDC");
    legend->SetFillColor(0);
    legend->SetFillStyle(0);
    legend->SetLineWidth(0);
    legend->SetTextSize(0.03);
    legend->SetName("legend");
    
    legend->AddEntry(gr2,entry2);
    legend->AddEntry(gr1,entry1);


    can->cd();


    //gr->Draw("A4");
    gr1->Draw("C hist");
    gr2->Draw("C hist same");
    //gr2->Draw("C hist same");
    //gr1->Draw("C hist");
    //gr2->Draw("C hist same");
    //gr->Draw("a4 same");
    
    

    legend->Draw("same");
}

void compCorrPairs() {

    
    TString casesNames[_N_CASES_] = {"caseA","caseB","caseC"};
    
    
    TH1D* hist[_N_CASES_][_N_VERSIONS_];
    TH1D* histNorm[_N_CASES_];
    TFile* fileIn[_N_CASES_][_N_VERSIONS_];
    TFile* fileInNorm[_N_CASES_];
    TCanvas* can[_N_CASES_][_N_VERSIONS_];

    Double_t nEventsNorm[_N_CASES_] = {9.97*10e6, 10e6, 10e6};
    Double_t nEvents = 5e3 * 200;

    for(int i = 0; i < _N_CASES_; i++) {
        fileInNorm[i] = new TFile(Form("./outputTherm/outputCorrPairs/%snorm2.root",casesNames[i].Data()));
        histNorm[i] = (TH1D*)fileInNorm[i]->Get("PiMinusProt")->Clone(Form("PiMinusProt%snorm",casesNames[i].Data()));
       // histNorm[i]->Scale(1.0/nEventsNorm[i]);
        for(int j = 0; j < _N_VERSIONS_; j++){
            fileIn[i][j] = new TFile(Form("./outputTherm/outputCorrPairs/%s%i.root",casesNames[i].Data(),j+1));
            hist[i][j] = (TH1D*)fileIn[i][j]->Get("PiMinusProt")->Clone(Form("PiMinusProt%s%i",casesNames[i].Data(),j+1));
            //hist[i][j]->Scale(1.0/nEvents);
            can[i][j] = new TCanvas(Form("PiMinusProt%s%i",casesNames[i].Data(),j),Form("PiMinusProt%s%i",casesNames[i].Data(),j+1),1000,1000);
            doublePlot(hist[i][j],histNorm[i],can[i][j],"newTherm","oldTherm");
            can[i][j]->SaveAs(Form("./outputTherm/outputCorrPairs/%s%i.png",casesNames[i].Data(),j+1));
        }

    }







}