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

#define _N_HIST_ 6
#define _N_PARTICLES_ 5

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

void seekParticles(std::vector<ParticleCoor> *particles[], ParticleCoor &part) {
    //int pdg = part.pid;
    switch(part.pid) {
        case -211: //piMinus
            particles[0]->push_back(part);
            break;
        case 111: //piZer
            particles[1]->push_back(part);
            break;
        case 211: //piPlus
            particles[2]->push_back(part);
            break;
        case 2212: //proton
            particles[3]->push_back(part);
            break;
        case 2112: //neutron
            particles[4]->push_back(part);
            break;
        default: 
            break;
    }
    

}

void fillHistogram(std::vector<ParticleCoor> *particles[], TH1D *hist[]) {


    TLorentzVector v1,v2;
    
    for(auto pion: (*particles[0])) {
        v1.SetPxPyPzE(pion.px,pion.py,pion.pz,pion.e);
        v1.Boost(0.,0.,0.629);
        for(auto prot: (*particles[3])){ //PiMinusProt
            if(pion.fathereid == prot.fathereid) {  
                v2.SetPxPyPzE(prot.px,prot.py,prot.pz,prot.e);
                v2.Boost(0.,0.,0.629);
                hist[0]->Fill((v1+v2).M());
            }
        }
        for(auto neut:  (*particles[4])){ //PiMinusNeutron
            if(pion.fathereid == neut.fathereid) {  
                v2.SetPxPyPzE(neut.px,neut.py,neut.pz,neut.e);
                v2.Boost(0.,0.,0.629);
                hist[1]->Fill((v1+v2).M());
            }
        }
    }

    for(auto pion: (*particles[1])) {
        v1.SetPxPyPzE(pion.px,pion.py,pion.pz,pion.e);
        v1.Boost(0.,0.,0.629);
        for(auto prot:  (*particles[3])){ //PiZeroProt
            if(pion.fathereid == prot.fathereid) {  
                v2.SetPxPyPzE(prot.px,prot.py,prot.pz,prot.e);
                v2.Boost(0.,0.,0.629);
                hist[2]->Fill((v1+v2).M());
            }
        }
        for(auto neut:  (*particles[4])){ //PiZeroNeutron
            if(pion.fathereid == neut.fathereid) {  
                v2.SetPxPyPzE(neut.px,neut.py,neut.pz,neut.e);
                v2.Boost(0.,0.,0.629);
                hist[3]->Fill((v1+v2).M());
            }
        }
    }

    for(auto pion: (*particles[2])) {
        v1.SetPxPyPzE(pion.px,pion.py,pion.pz,pion.e);
        v1.Boost(0.,0.,0.629);
        for(auto prot:  (*particles[3])){ //PiPlusProt
            if(pion.fathereid == prot.fathereid) {  
                v2.SetPxPyPzE(prot.px,prot.py,prot.pz,prot.e);
                v2.Boost(0.,0.,0.629);
                hist[4]->Fill((v1+v2).M());
            }
        }
        for(auto neut:  (*particles[4])){ //PiPlusNeutron
            if(pion.fathereid == neut.fathereid) {  
                v2.SetPxPyPzE(neut.px,neut.py,neut.pz,neut.e);
                v2.Boost(0.,0.,0.629);
                hist[5]->Fill((v1+v2).M());
            }
        }
    }










    /*TLorentzVector v1,v2;
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


        
    }*/


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


    TH1D* hist[_N_HIST_];
    TString titles[_N_HIST_] = {"PiMinusProt","PiMinusNeutron","PiZeroProton","PiZeroNeutron","PiPlusProton","PiPlusNeutron"};
    Int_t XBins = 1000;
    Double_t XMin = 1.0;
    Double_t XMax = 1.8;

    for(int i = 0; i < _N_HIST_; i++) {
        hist[i] = new TH1D(titles[i].Data(),titles[i].Data(),XBins,XMin,XMax);
        hist[i]->GetXaxis()->SetTitle("M (GeV/c^{2})");
        hist[i]->GetYaxis()->SetTitle("dN/dM (GeV/c^{2})^{-1}");

    }

    /*

    TH1D* H1D = new TH1D("PiMinusProt","PiMinusProt",XBins,XMin,XMax);
    H1D->GetXaxis()->SetTitle("M (GeV/c^{2})");
    H1D->GetYaxis()->SetTitle("dN/dM (GeV/c^{2})^{-1}");*/



    Chain->GetEntry(0); 
	Int_t event_prev = Particle.eventid; //first event number
    Int_t event;
    int e = 0;
    //vector<ParticleCoor> particles;
    std::vector<ParticleCoor> *particles[_N_PARTICLES_];
    for(int i = 0; i < _N_PARTICLES_; i++)
        particles[i] = new std::vector<ParticleCoor>;
    
    //cout<<Chain->GetEntries()<<endl;
    for (Int_t i = 0; i < Chain->GetEntries() ; i++) { 
		Chain->GetEntry(i);
        event = Particle.eventid;
        if(event == event_prev)
            seekParticles(particles,Particle);


        if (event != event_prev){
            e++;
            fillHistogram(particles,hist);	
            event_prev = event;
            for(int j = 0; j < _N_PARTICLES_; j++)
                particles[j]->clear();
            Chain->GetEntry(i);
            seekParticles(particles,Particle);
            if (e%10000==0)
                cout<<"Event: "<<e<<endl;
        }
        //if(e==3)
            //break;
    }
    
    TFile *file = new TFile(aOutDir.Data(),"RECREATE");
	file->cd();
    for(int i = 0; i < _N_HIST_; i++)
        hist[i]->Write();
    file->Save();
	file->Close();

}


int makeCorrPairsAll(TString inputlist = "", TString outfile = "", Long64_t nDesEvents = -1, Int_t maxFiles = -1) {

    gROOT -> SetBatch(kTRUE);

    createHistograms("/lustre/hades/user/mkurach/newTherm/caseB07/H225E0D4",outfile,maxFiles,inputlist);

	gROOT -> SetBatch(kFALSE);
    return 0;

}