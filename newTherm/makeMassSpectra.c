#include <iostream>
//#include <TFile.h>
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
#include "TSystemDirectory.h"
#include "TList.h"
#include <vector>
#include "events2chain.C"
//#include "drawStyle.C"
#include <map>

#define _N_PARTICLES_ 15

TChain* addMultEvents2Chain(TString commaSeparatedList, ParticleCoor* aParticle, Int_t* aEvents){

    static StructEvent tStructEvents;
    TChain* tChainParts = new TChain(_PARTICLES_TREE_); //particle chain
    TChain* tChainEvent = new TChain(_EVENTS_TREE_);  //event chain

    tChainParts->SetBranchAddress(_PARTICLE_BRANCH_, aParticle); //linking branch from tree to chain
    tChainEvent->SetBranchAddress(_EVENTS_BRANCH_, &tStructEvents);
  

    TObjArray* arr = commaSeparatedList.Tokenize(",");
    Int_t n = arr->GetEntries(); //number of files loaded

    if(n == 0) {
        cout<<"No files found! Are u blind? Check paths bro!"<<endl;
        delete arr;
        return nullptr; 
    } 
    else {
        for(Int_t i = 0; i < n; i ++){ //for every file...

            TString name = ((TObjString*)arr->At(i))->GetString(); //get it's name...
            
            cerr << "Adding file: " << name << endl;
            if (!tChainParts->Add(name) || !tChainEvent->Add(name)){ //if any of TChains doesn't get filled properlu (Add() method returns 0) -> kill it.
                cout<<"cannot open file "<<name<<"! It doesn't exist or might be broken. TChain gets killed. Whatcha starin' @ m8? Go fix it!"<<endl;
                delete arr;
                return nullptr;
            }
        }
        //if all files added correctly -> print some stuff
        (*aEvents) = tChainEvent->GetEntries();

        cerr<<"Total number of events: "<<(*aEvents)<<endl;
        return tChainParts;
    }
    delete arr;
    return nullptr;

}

Int_t pid[_N_PARTICLES_] = {12112, 1214, 22112, 32214, 2122, 32212, 2216, 12216,5218, 2128,32124,22124, 9299,2218,12214};
std::map<Int_t,Int_t> pidMap;

void fillHistogram(ParticleCoor particle, TH1D *hist[]) {

    Int_t pdg = particle.pid;
    Int_t index = pidMap[TMath::Abs(pdg)];
    TLorentzVector v;
    v.SetPxPyPzE(particle.px,particle.py,particle.pz,particle.e);
    v.Boost(0.,0.,0.629);
     if(  TMath::Abs(pdg) == 12112||TMath::Abs(pdg) == 1214||TMath::Abs(pdg) == 22112||TMath::Abs(pdg) == 32214||TMath::Abs(pdg) == 2122||TMath::Abs(pdg) == 32212||TMath::Abs(pdg) == 2216||TMath::Abs(pdg) == 12216 ||TMath::Abs(pdg) == 5218||TMath::Abs(pdg) == 2128||TMath::Abs(pdg) ==32124||TMath::Abs(pdg) ==22124||TMath::Abs(pdg) ==9299||TMath::Abs(pdg) ==2218||TMath::Abs(pdg) ==12214) {
            if(particle.pid == particle.fatherpid)
                hist[index]->Fill(v.M());
            }

}
void createHistograms(TString aEventDir, TString aOutDir, Int_t aEventFiles, TString inputList) {

    //READING FILES 
    
    cout << "Reading events in folder: "<<aEventDir.Data()<< endl;
   	gSystem->MakeDirectory(aOutDir);
   	gROOT->SetBatch(kTRUE);
	
	static ParticleCoor Particle;
	static StructEvent tStructEvents;
    Int_t   Events;
	TChain* Chain;

    	if (aEventFiles == -1){ //use this to ruch on a batch farm (slurm)
        Chain = addMultEvents2Chain(inputList,&Particle,&Events);

        if(Chain==nullptr){
            cout<<"No TChain created! Killing analysis"<<endl;
            gROOT -> SetBatch(kFALSE);
            exit(1);
        }
    }
    else { //use for testing purposes, number of files you specify by changing maxFiles default value
        Int_t nFiles = 0;

        Chain = new TChain(_PARTICLES_TREE_);
        TChain* tChainEvent = new TChain(_EVENTS_TREE_); 

        Chain->SetBranchAddress(_PARTICLE_BRANCH_, &Particle); //linking branch from tree to chain
        tChainEvent->SetBranchAddress(_EVENTS_BRANCH_, &tStructEvents);

        TString inputFolder = aEventDir; //path to your data directory

        TSystemDirectory* inputDir = new TSystemDirectory("inputDir", inputFolder);
        TList* files = inputDir->GetListOfFiles(); //gets list of files in designated location

        cout<<"Files loaded:"<<endl;

        if(files==nullptr){
            cout<<"co jest kurwa"<<endl;
            exit(1);
        }

        for (Int_t i = 0; i <= files->LastIndex() && nFiles < aEventFiles; i++){           
            if (((TSystemFile*) files->At(i))->IsDirectory()) continue; //if one of the "files" is directory (directory inside directory) -> skip it.
           
            if(TString(((TSystemFile*)files->At(i))->GetName())){
                TString fileName = inputFolder + "/" + ((TSystemFile*) files->At(i))->GetName(); //get file name
                if(!fileName.Contains("root") || !fileName.Contains("event")) continue; //if non-root file or some other bullshit -> skip it

                Chain->Add(fileName);
                tChainEvent->Add(fileName);

                cout<<fileName<<endl;
                nFiles++;
            }
        }

        Events = tChainEvent->GetEntries();
        cerr << "Total number of events: " << Events << endl;
        delete tChainEvent;
    }
    //N*0(1440),N*0(1520),N*0(1535),Delta0(1600),Delta0(1620),N*0(1650),N*0(1675), N*0(1680), 
    //TString particlesTitles[_N_PARTICLES_] = {"N(1440)^{*0}"};
    TString particlesTitlesDiff[_N_PARTICLES_] = {"Ns1440zer","Ns1520zer","Ns1535zer","Dl1600zer","Dl1620zer","Ns1650zer","Ns1675zer","Ns1680zer","Ns2250zer","Ns2190zer","Ns1720zer","Ns1700zer","Dl2420zer","Dl1950zer","Dl1700zer"};

    TH1D* H1D[_N_PARTICLES_];
    Int_t XBins = 1000;
    Double_t XMin = 0.0;
    Double_t XMax = 3.0;



    for(int i = 0; i < _N_PARTICLES_; i++){
        H1D[i] = new TH1D(particlesTitlesDiff[i].Data(),"",XBins,XMin,XMax);
        H1D[i]->GetXaxis()->SetTitle("M (GeV/c^{2})");
        H1D[i]->GetYaxis()->SetTitle("dN/dM (GeV/c^{2})^{-1}");

        pidMap[pid[i]] = i;
    }


    Chain->GetEntry(0); 
	Int_t event_prev = Particle.eventid; //first event number
    Int_t event;

    int e = 0;

    for (Int_t i = 0; i < Chain->GetEntries() ; i++) { 
		Chain->GetEntry(i);
        event = Particle.eventid;

        fillHistogram(Particle,H1D); 

        if (event != event_prev){
            e++; 	
            event_prev = event;
            if (e%10000==0)
                cout<<"Event: "<<e<<endl;
        }
    }
    
    TFile *file = new TFile(aOutDir.Data(),"RECREATE");
	file->cd();
    for(int i = 0; i < _N_PARTICLES_; i++)
        H1D[i]->Write();
    file->Save();
	file->Close();

}


int makeMassSpectra(TString inputlist = "", TString outfile = "", Long64_t nDesEvents = -1, Int_t maxFiles = -1) {

    gROOT -> SetBatch(kTRUE);

    createHistograms("/u/mkurach/lustre/hades/user/kjedrzej/HubbDeltMotoMoto/H165E6D4femto",outfile,maxFiles,inputlist);

	gROOT -> SetBatch(kFALSE);
    return 0;

}