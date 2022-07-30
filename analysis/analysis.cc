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
#include <numeric>
using namespace std;


void analysis()
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
  map<int,double> detpx1;
  map<int,double> detpy1;
  map<int,double> detpz1;

  for (int i =0;i< (int)Detectors.size();i++){
    j = 0;
    Detector_number = 0;
    Detector_r.clear();
    Detector_z.clear();
   
    stringstream ss(Detectors[i]);
    while (ss >> s) {
      j++;
      if (j==1) Detector_number += atoi(s.c_str())*100;    // C module
      else if (j==2) Detector_number += atoi(s.c_str())*10;// P string position
      else if (j==3) Detector_number += atoi(s.c_str())*1; // D detector position
      else if (j==4) Detector_position.SetX(atof(s.c_str())*mm); // x
      else if (j==5) Detector_position.SetY(atof(s.c_str())*mm); // y
      else if (j==6) Detector_position.SetZ(atof(s.c_str())*mm); // z
      else if (j==7) Detector_name = s;                          // name
      else if (j==8) Detector_slices = atoi(s.c_str());          //  2: BEGe (broad energy germanium) detector 3:Ortec
      else if (j==9) Detector_r.push_back((double) atof(s.c_str())*mm);  // 2: radius at bottom 3: radius at bottom ring
      else if (j==10) Detector_z.push_back((double) atof(s.c_str())*mm); // 2: bottom 3: bottom ring
      else if (j==11) Detector_r.push_back((double) atof(s.c_str())*mm); // 2: radius at top 3: raidus at middle ring
      else if (j==12) Detector_z.push_back((double) atof(s.c_str())*mm); // 2: top 3:middle ring
      else if (j==13) Detector_r.push_back((double) atof(s.c_str())*mm); // 3: radius at top ring
      else if (j==14) Detector_z.push_back((double) atof(s.c_str())*mm); // 3: top ring
    }
    if(Detector_slices==2){
      Detector_r.push_back(Detector_r[1]+0.00001*mm); // dead layer thickness=0.00001 mm
      Detector_z.push_back(Detector_z[1]+0.00001*mm); // dead layer thickness=0.00001 mm
    }
    Detector_position.SetZ(Detector_position.z() + Detector_z[2] - 20*cm); // the top position of the detector move to lower 20 cm position.
    string detName = "Det_"+Detector_name;
    detname[detName]=Detector_number;
    detpx[detName]=Detector_position.X();
    detpy[detName]=Detector_position.Y();
    detpz[detName]=Detector_position.Z();
    detpx1[Detector_number]=Detector_position.X();
    detpy1[Detector_number]=Detector_position.Y();
    detpz1[Detector_number]=Detector_position.Z();
  }

  TChain* fTree = new TChain("MonopoleTree"); 
  fTree->Add("/global/homes/p/pchu/work/g4work/monopole_MJD/analysis/monopole.root");
  ofstream fout("analysis.csv",ios::app); 
  fout.precision(12);
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

  Int_t TEventID;
  Double_t TEsum;
  vector<Double_t> TEhit;
  vector<Int_t> TChannel;
  vector<Double_t> TStartTime;
  vector<Double_t> TStartTime_sorted;
  vector<Double_t> TEndTime;
  vector<Double_t> TDetX;
  vector<Double_t> TDetY;
  vector<Double_t> TDetZ;
  
  /*
  TFile *hfile = TFile::Open("build.root","RECREATE");
  TTree *tree = new TTree("monopole","");
  tree->Branch("TEventID",&TEventID,"TEventID/I");
  tree->Branch("TEsum",&TEsum,"TEsum/D");
  tree->Branch("TEhit",&TEhit);
  tree->Branch("TChannel",&TChannel);
  tree->Branch("TStartTime",&TStartTime);
  tree->Branch("TEndTime",&TEndTime);
  */
  Double_t pid,z,a,charge,energy,Edep,timestamp,weight,posx,posy,posz,velocity,u,v,w,detPx,detPy,detPz;
  Int_t trackid,parentid,channel;
  string det_name;

  fTree->GetEntry(0);
  Int_t eventid=EventID;
  Int_t eventid_last = EventID;
  Int_t channel_idx;
  Double_t Esum = 0;
  Double_t starttime=Time;
  Double_t endtime = 0;
  vector<Double_t> energy_event;
  vector<Double_t> time_event;
  vector<Double_t> time_sorted_event;
  vector<Int_t> channel_event;

  Double_t energy_hit[58];
  Double_t starttime_hit[58];
  Double_t endtime_hit[58];

  for (Int_t j=0;j<nChannel;j++){
    energy_hit[j] = 0;
    starttime_hit[j] = 0;
    endtime_hit[j] = 0;
  }
  
  for (Int_t i=0;i<nentries;i++){
  //for (Int_t i=0;i<1000;i++){
    fTree->GetEntry(i);    
    pid=PID;
    z=Z;
    a=A;
    charge = Charge;
    energy = Energy;
    Edep = EDeposit;
    timestamp = Time;
    weight = Weight;
    posx = PosX;
    posy = PosY;
    posz = PosZ;
    trackid = TrackID;
    parentid = ParentID;
    velocity = Velocity;
    u = U;
    v = V;
    w = W;
    eventid = EventID;
    det_name = VolumeName;
    channel = detname[det_name];
    detPx = detpx[det_name];
    detPy = detpy[det_name];
    detPz = detpz[det_name];
    //cout << i << "," << det_name << endl;

    if (eventid!=eventid_last){
      //cout << eventid_last << "," << eventid << endl;
      if(energy_event.size()>0){

	// Initialize sum energy
	Esum = 0;

	// Sort time; get start/end time
	sort(time_sorted_event.begin(), time_sorted_event.end());
	starttime = time_sorted_event.at(0);
	endtime = time_sorted_event.at(time_sorted_event.size()-1);

	// Initialize hit energy and start/end time
	for (Int_t j=0;j<nChannel;j++){
	  energy_hit[j] = 0;
	  endtime_hit[j] = starttime;
	  starttime_hit[j] = endtime;
	}
	
	// Loop steps; sum hit energy
	for(size_t j = 0;j<energy_event.size();j++){
	  Esum = Esum + energy_event.at(j);	
	  channel_idx = detindex[channel_event.at(j)];
	  energy_hit[channel_idx] = energy_hit[channel_idx]+energy_event.at(j);
	  if (time_event.at(j)>endtime_hit[channel_idx]){
	    endtime_hit[channel_idx] = time_event.at(j);
	  }
	  if(time_event.at(j)<starttime_hit[channel_idx]){
	    starttime_hit[channel_idx] = time_event.at(j);
	  }
	}

	if(Esum>0){
          TEventID = eventid_last;
          TEsum = Esum;                                                                                                                                                                                   
	  //Save non-zero hits
          for(Int_t k=0;k<nChannel;k++){
	    if(energy_hit[k]>0){
	      TChannel.push_back(Chan[k]);
	      TEhit.push_back(energy_hit[k]);
	      TStartTime.push_back(starttime_hit[k]);
	      TStartTime_sorted.push_back(starttime_hit[k]);
	      TEndTime.push_back(endtime_hit[k]);
	      TDetX.push_back(detpx1[Chan[k]]);
	      TDetY.push_back(detpy1[Chan[k]]);
	      TDetZ.push_back(detpz1[Chan[k]]);
	    }
	  }


	  if(TEhit.size()>1){//Need multiplicity >=2 	    
	    //Calculate average hit energy
	    Double_t Ehit_sum = accumulate(TEhit.begin(),TEhit.end(),0.0);
            Double_t TEhit_mean = Ehit_sum/TEhit.size();
            Double_t Ehit_sq_sum = inner_product(TEhit.begin(),TEhit.end(),TEhit.begin(),0.0);
            Double_t TEhit_stdev = sqrt(Ehit_sq_sum/TEhit.size()-TEhit_mean*TEhit_mean);

	    //Calculate time difference
	    sort(TStartTime_sorted.begin(), TStartTime_sorted.end());
	    vector<Double_t> TimeDiff;
	    for(size_t j=0;j<TStartTime_sorted.size()-1;j++){
	      Double_t t1 = TStartTime_sorted.at(j);
	      Double_t t2 = TStartTime_sorted.at(j+1);
	      Double_t tdiff = t2-t1;
	      TimeDiff.push_back(tdiff);
	    }
	    
	    Double_t tdiff_sum = accumulate(TimeDiff.begin(),TimeDiff.end(),0.0);
	    Double_t TTdiff_mean = tdiff_sum/TimeDiff.size();
	    Double_t tdiff_sq_sum = inner_product(TimeDiff.begin(),TimeDiff.end(),TimeDiff.begin(),0.0);
	    Double_t TTdiff_stdev = sqrt(tdiff_sq_sum/TimeDiff.size()-TTdiff_mean*TTdiff_mean);
	    TimeDiff.clear();

	    //Calculate costheta and phi
	    vector<Double_t> CosTheta;
	    vector<Double_t> Phi;
	    for(size_t j = 0;j<TDetX.size();j++){	      
	      for(size_t k = j+1;k<TDetX.size();k++){
		Double_t x1 = TDetX.at(j);
		Double_t x2 = TDetX.at(k);
		Double_t y1 = TDetY.at(j);
		Double_t y2 = TDetY.at(k);
		Double_t z1 = TDetZ.at(j);
		Double_t z2 = TDetZ.at(k);
		Double_t costhetajk = (z1-z2)/sqrt(pow(x1-x2,2)+pow(y1-y2,2)+pow(z1-z2,2));
		CosTheta.push_back(costhetajk);
		Double_t phijk=0;
		if(x1!=x2){
		  phijk = atan((y1-y2)/(x1-x2));
		}else{
		  phijk=0;
		}
		Phi.push_back(phijk);
	      }
	    }
	    Double_t costheta_sum = accumulate(CosTheta.begin(), CosTheta.end(), 0.0);
	    Double_t TCostheta_mean = costheta_sum / CosTheta.size();

	    Double_t costheta_sq_sum = inner_product(CosTheta.begin(), CosTheta.end(), CosTheta.begin(), 0.0);
	    Double_t TCostheta_stdev = sqrt(costheta_sq_sum /CosTheta.size() - TCostheta_mean*TCostheta_mean);

	    Double_t phi_sum = accumulate(Phi.begin(), Phi.end(), 0.0);
            Double_t TPhi_mean = phi_sum/Phi.size();

            Double_t phi_sq_sum = inner_product(Phi.begin(), Phi.end(), Phi.begin(), 0.0);
            Double_t TPhi_stdev = sqrt(phi_sq_sum /Phi.size() - TPhi_mean * TPhi_mean);	    
	    CosTheta.clear();
	    Phi.clear();
	    // Print out
	    for(size_t j=0;j<TEhit.size();j++){
	      fout << eventid_last << "," << setprecision(5) << Esum << ","<< TEhit.size() << "," << setprecision(12) << TTdiff_mean << "," << TTdiff_stdev << "," << TEhit_mean <<"," << TEhit_stdev << "," << TCostheta_mean << "," << TCostheta_stdev << "," << TPhi_mean << "," << TPhi_stdev << "," << TChannel.at(j) << "," << TEhit.at(j) << "," << setprecision(12) << TStartTime.at(j) << "," << TEndTime.at(j)<<"," << TDetX.at(j) << "," << TDetY.at(j) << "," << TDetZ.at(j) << endl; 
	    }
	  }
	}
      }
      energy_event.clear();
      time_event.clear();
      time_sorted_event.clear();
      channel_event.clear();

      TChannel.clear();
      TEhit.clear();      
      TStartTime.clear();
      TStartTime_sorted.clear();
      TEndTime.clear();
      TDetX.clear();
      TDetY.clear();
      TDetZ.clear();
      if(det_name!="world"){
	energy_event.push_back(Edep);
	time_event.push_back(timestamp);
	time_sorted_event.push_back(timestamp);
	channel_event.push_back(channel);
      }
    }else if(eventid==eventid_last){
      if(det_name!="world"){
	energy_event.push_back(Edep);
	time_event.push_back(timestamp);
	time_sorted_event.push_back(timestamp);
	channel_event.push_back(channel);
      }
    }     
    eventid_last=eventid;
  }
  //  hfile->Write();

  return 0;

}
