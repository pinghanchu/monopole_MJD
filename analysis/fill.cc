#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TProof.h"
#include "TROOT.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TObject.h"
#include "TString.h"
#include <map>
#include <iostream>
#include <string>
#include <vector>
#include <stdio.h>
#include <getopt.h>

using namespace std;

void fill()
{

  ifstream fin("/global/homes/p/pchu/work/g4work/monopole_MJD/Detectorposition.txt");
  vector<string> Detectors;
  string input;
  while(getline(fin,input)){
    Detectors.push_back(input);
  }
  fin.close();

  const Int_t nChannel=58;
  Int_t Chan[nChannel]   = {111,112,113,114,121,122,123,124,131,132,133,134,141,142,143,144,145,151,152,153,154,161,162,163,164,171,172,173,174, 211,212,213,214,221,222,223,224,225,231,232,233,241,242,243,244,245,251,252,253,254,261,262,263,264,271,272,273,274};
  Int_t GoodBad[nChannel]= {1 , 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,1 , 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  map<int,int> detindex;
  for(int i=0;i<nChannel;i++){
    detindex[Chan[i]]=i;
  }

  double cm=0.01;
  double mm=0.001;
  string s;
  int j=0;
  int Detector_number;
  string Detector_name;
  string Detector_name_sol;
  string Detector_name_log;
  TVector3 Detector_position;
  int Detector_slices=0;
  vector<double> Detector_r;
  vector<double> Detector_z;
  map<string,int> detname;
  map<string,double> detpx;
  map<string,double> detpy;
  map<string,double> detpz;
  for (int i =0;i< (int)Detectors.size();i++){
    j = 0;
    Detector_number = 0;
    Detector_r.clear();
    Detector_z.clear();

    stringstream ss(Detectors[i]);
    while (ss >> s) {
      j++;
      if (j==1) Detector_number += atoi(s.c_str())*100;
      else if (j==2) Detector_number += atoi(s.c_str())*10;
      else if (j==3) Detector_number += atoi(s.c_str())*1;
      else if (j==4) Detector_position.SetX(atof(s.c_str())*mm);
      else if (j==5) Detector_position.SetY(atof(s.c_str())*mm);
      else if (j==6) Detector_position.SetZ(atof(s.c_str())*mm);
      else if (j==7) Detector_name = s;
      else if (j==8) Detector_slices = atoi(s.c_str());
      else if (j==9) Detector_r.push_back((double) atof(s.c_str())*mm);
      else if (j==10) Detector_z.push_back((double) atof(s.c_str())*mm);
      else if (j==11) Detector_r.push_back((double) atof(s.c_str())*mm);
      else if (j==12) Detector_z.push_back((double) atof(s.c_str())*mm);
      else if (j==13) Detector_r.push_back((double) atof(s.c_str())*mm);
      else if (j==14) Detector_z.push_back((double) atof(s.c_str())*mm);
    }
    if(Detector_slices==2){
      Detector_r.push_back(Detector_r[1]+0.00001*mm);
      Detector_z.push_back(Detector_z[1]+0.00001*mm);
    }
    Detector_position.SetZ(Detector_position.z() + Detector_z[2] - 20*cm);
    string detName = "Det_"+Detector_name;
    detname[detName]=Detector_number;
    detpx[detName]=Detector_position.X();
    detpy[detName]=Detector_position.Y();
    detpz[detName]=Detector_position.Z();
  }

  TChain* fTree = new TChain("MonopoleTree"); 
  fTree->Add("/global/homes/p/pchu/work/g4work/monopole_MJD/analysis/monopole.root");
  
  
  cout << "Added files" << endl;
  cout << "-------" << endl;
  Long64_t nentries = (Long64_t)fTree->GetEntries();
  cout << fTree->GetNtrees() << " trees with " << nentries << " entries " <<endl;
  cout << "-------" << endl;
  Double_t PID;
  Double_t Z;
  Double_t A;
  Double_t Charge;
  Double_t Energy;
  Double_t EDeposit;
  Double_t Time;
  Double_t Weight;
  Double_t PosX;
  Double_t PosY;
  Double_t PosZ;
  Int_t TrackID;
  Int_t ParentID;
  Double_t Velocity;
  Double_t U;
  Double_t V;
  Double_t W;
  Int_t EventID;
  Char_t VolumeName[25];

  fTree->SetBranchStatus("*",1);
  fTree->SetBranchAddress("PID",&PID);
  fTree->SetBranchAddress("Z",&Z);
  fTree->SetBranchAddress("A",&A);
  fTree->SetBranchAddress("Charge",&Charge);
  fTree->SetBranchAddress("Energy",&Energy);
  fTree->SetBranchAddress("EDeposit",&EDeposit);
  fTree->SetBranchAddress("Time",&Time);
  fTree->SetBranchAddress("Weight",&Weight);
  fTree->SetBranchAddress("PosX",&PosX);
  fTree->SetBranchAddress("PosY",&PosY);
  fTree->SetBranchAddress("PosZ",&PosZ);
  fTree->SetBranchAddress("TrackID",&TrackID);
  fTree->SetBranchAddress("ParentID",&ParentID);
  fTree->SetBranchAddress("Velocity",&Velocity);
  fTree->SetBranchAddress("U",&U);
  fTree->SetBranchAddress("V",&V);
  fTree->SetBranchAddress("W",&W);
  fTree->SetBranchAddress("EventID", &EventID);
  fTree->SetBranchAddress("VolumeName",VolumeName);
  Int_t eventid_old=0;
  Int_t Eventcount = 0;

  ofstream fout("./data.csv",ios::app);
  fout.precision(15);
  fout << "EventID,Time,Channel,Ekin,Edep,PosX,PosY,PosZ,Velocity,U,V,W,DetPX,DetPY,DetPZ"<< endl;
  for (Int_t i=0;i<nentries;i++){
    fTree->GetEntry(i);    
    Double_t pid=PID;
    Double_t z=Z;
    Double_t a=A;
    Double_t charge = Charge;
    Double_t energy = Energy;
    Double_t Edep = EDeposit;
    Double_t time = Time;
    Double_t weight = Weight;
    Double_t posx = PosX;
    Double_t posy = PosY;
    Double_t posz = PosZ;
    Int_t trackid = TrackID;
    Int_t parentid = ParentID;
    Double_t velocity = Velocity;
    Double_t u = U;
    Double_t v = V;
    Double_t w = W;
    Int_t eventid = EventID;
    string det_name = VolumeName;
    Int_t channel = detname[det_name];
    Double_t detPx = detpx[det_name];
    Double_t detPy = detpy[det_name];
    Double_t detPz = detpz[det_name];
    //cout << eventid << "," << det_name << endl;
    if(det_name!="world"){
      fout << eventid << "," << time << "," << channel << "," << energy <<","<<Edep << "," << posx << "," << posy << "," << posz << "," << velocity << "," << u << "," << v << ","<< w << ","<< detPx << "," << detPy << "," << detPz << endl;
    }
  }
  return 0;
}