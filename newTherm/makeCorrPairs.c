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

//||Abs(pdg) == 5218||Abs(pdg) == 2128||Abs(pdg) ==32124||Abs(pdg) ==22124||Abs(pdg) == 9299||Abs(pdg) ==2218||Abs(pdg) ==12214
void seekParticles(vector<ParticleCoor> &particles, ParticleCoor &part, vector<ParticleCoor> &mothers) {
    if(((part.pid == 2212) && ( part.fatherpid == 2114 || part.fatherpid == 12112 || part.fatherpid == 1214 || part.fatherpid ==22112|| part.fatherpid ==32214|| part.fatherpid ==2122|| part.fatherpid ==32212|| part.fatherpid ==2216|| part.fatherpid ==12216|| part.fatherpid == 5218|| part.fatherpid ==2128|| part.fatherpid ==32124|| part.fatherpid ==22124|| part.fatherpid == 9299|| part.fatherpid ==2218|| part.fatherpid ==12214)) || ((part.pid == -211) && ( part.fatherpid == 2114 || part.fatherpid == 12112 || part.fatherpid == 1214 || part.fatherpid ==22112|| part.fatherpid ==32214|| part.fatherpid ==2122|| part.fatherpid ==32212|| part.fatherpid ==2216|| part.fatherpid ==12216|| part.fatherpid == 5218|| part.fatherpid ==2128|| part.fatherpid ==32124|| part.fatherpid ==22124|| part.fatherpid == 9299|| part.fatherpid ==2218|| part.fatherpid ==12214)))
        particles.push_back(part); 
    if ( (part.pid == 2114 && part.fatherpid == 2114) ||(part.pid == 12112 && part.fatherpid == 12112)|| (part.pid == 1214 && part.fatherpid == 1214)|| (part.pid == 22112 && part.fatherpid == 22112) || (part.pid ==32214 && part.fatherpid ==32214)|| (part.pid ==2122 && part.fatherpid ==2122)|| (part.pid ==32212 && part.fatherpid ==32212) || (part.pid ==2216 && part.fatherpid ==2216)|| (part.pid ==12216 && part.fatherpid ==12216)   || (part.pid == 5218 && part.fatherpid == 5218)|| (part.pid ==2128 && part.fatherpid ==2128)|| (part.pid ==32124 &&part.fatherpid ==32124)|| (part.pid ==22124 &&part.fatherpid ==22124)|| (part.pid == 9299 &&part.fatherpid == 9299)|| (part.pid ==2218 && part.fatherpid ==2218)|| (part.pid ==12214 && part.fatherpid ==12214))
        mothers.push_back(part);
}

void fillHistogram(vector<ParticleCoor> &particles, TH1D *hist, vector<ParticleCoor> &mothers) {

    TLorentzVector v1,v2;
    Int_t fatherEid;

    for(int i = 0; i < particles.size(); i++) {
        fatherEid = particles[i].fathereid;
        for(int j = 0; j < particles.size(); j++) {
            if((particles[j].fathereid == fatherEid) && (j != i)) {
                for(int k = 0; k < mothers.size();k++) {
                    if(mothers[k].eid == fatherEid) {
                        v1.SetPxPyPzE(particles[i].px,particles[i].py,particles[i].pz,particles[i].e);
                        v2.SetPxPyPzE(particles[j].px,particles[j].py,particles[j].pz,particles[j].e);
                        v1.Boost(0.,0.,0.629);
                        v2.Boost(0.,0.,0.629);
                        hist->Fill((v1+v2).M(), 0.5);
                    }
                }

                //particles.erase(particles.begin()+j);
                //break;

            }
        }
       // particles.erase(particles.begin()+i);


        
    }


    /*for(auto particle : particles) {
        cout<<"Particle pid: "<<particle.pid<<endl;
        cout<<"Particle father pid: "<<particle.fatherpid<<endl;
        cout<<"Particle root pid: "<<particle.rootpid<<endl;
        cout<<"Particle eid: "<<particle.eid<<endl;
        cout<<"Particle father eid: "<<particle.fathereid<<endl;
        cout<<"Particle event id : "<<particle.eventid<<endl;
        cout<<"Is decayed: "<<particle.decayed<<endl;
        cout<<endl;
    }
cout<<endl<<"matki"<<endl;
        for(auto particle : mothers) {
        cout<<"Particle pid: "<<particle.pid<<endl;
        cout<<"Particle father pid: "<<particle.fatherpid<<endl;
        cout<<"Particle root pid: "<<particle.rootpid<<endl;
        cout<<"Particle eid: "<<particle.eid<<endl;
        cout<<"Particle father eid: "<<particle.fathereid<<endl;
        cout<<"Particle event id : "<<particle.eventid<<endl;
        cout<<"Is decayed: "<<particle.decayed<<endl;
        cout<<endl;
    }*/

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

    Int_t XBins = 1000;
    Double_t XMin = 1.0;
    Double_t XMax = 1.8;

    TH1D* H1D = new TH1D("PiMinusProt","PiMinusProt",XBins,XMin,XMax);
    H1D->GetXaxis()->SetTitle("M (GeV/c^{2})");
    H1D->GetYaxis()->SetTitle("dN/dM (GeV/c^{2})^{-1}");



    Chain->GetEntry(0); 
	Int_t event_prev = Particle.eventid; //first event number
    Int_t event;
    int e = 0;
     vector<ParticleCoor> particles;
     vector<ParticleCoor> mothers;

    for (Int_t i = 0; i < Chain->GetEntries() ; i++) { 
		Chain->GetEntry(i);
        event = Particle.eventid;
        if(event == event_prev)
            seekParticles(particles,Particle,mothers);


        if (event != event_prev){
            e++;
            fillHistogram(particles,H1D,mothers);	
            event_prev = event;
            particles.clear();
            mothers.clear();
            Chain->GetEntry(i);
            seekParticles(particles,Particle,mothers);
            if (e%10000==0)
                cout<<endl<<"Event: "<<e<<endl;
        }
        //if(e==3)
           // break;
    }
    
    TFile *file = new TFile(aOutDir.Data(),"RECREATE");
	file->cd();
    H1D->Write();
    file->Save();
	file->Close();

}


int makeCorrPairs(TString inputlist = "", TString outfile = "", Long64_t nDesEvents = -1, Int_t maxFiles = -1) {

    gROOT -> SetBatch(kTRUE);

    createHistograms("/u/mkurach/lustre/hades/user/kjedrzej/HubbDeltMotoMoto/H165E6D4femto",outfile,maxFiles,inputlist);

	gROOT -> SetBatch(kFALSE);
    return 0;

}