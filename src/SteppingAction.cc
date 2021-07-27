//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file exoticphysics/monopole/src/SteppingAction.cc
/// \brief Implementation of the SteppingAction class
//
// $Id$
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"
#include "G4Step.hh"
#include "G4RunManager.hh"
#include "Run.hh"
#include "EventAction.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(EventAction* evt)
  : G4UserSteppingAction(),fEventAction(evt)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{


  G4int eventID = fEventAction->GetEventID();

  G4double edep = aStep->GetTotalEnergyDeposit();
  if (edep <= 0.) { return; }

  Run* run = static_cast<Run*>(
             G4RunManager::GetRunManager()->GetNonConstCurrentRun());
 
  
  //G4Event *event=run->GetCurrentEvent();
  //G4int eventId = event->GetEventID();

  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  //Bragg curve
  G4double x = aStep->GetPreStepPoint()->GetPosition().z();
  G4double dx = aStep->GetPostStepPoint()->GetPosition().z() - x;
  x += dx*G4UniformRand() - run->GetOffsetZ();
  run->FillHisto(1, x, edep);

  G4VPhysicalVolume* volume = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume();

  G4String volumename = volume->GetName();
  //G4String particleType = aStep->GetTrack()->GetDefinition()->GetParticleType();
  //G4String particleName = aStep->GetTrack()->GetDefinition()->GetParticleName();
  //G4ThreeVector position = aStep->GetTrack()->GetPosition();
  //G4double p_time = aStep->GetTrack()->GetGlobalTime();


  const G4ParticleDefinition* particle = aStep->GetTrack()->GetParticleDefinition();
  G4String name   = particle->GetParticleName();
  G4int pid       = particle->GetPDGEncoding();
  G4int Z         = particle->GetAtomicNumber();
  G4int A         = particle->GetAtomicMass();
  G4double charge = particle->GetPDGCharge();
  G4double energy = aStep->GetTrack()->GetKineticEnergy();
  G4double time   = aStep->GetTrack()->GetGlobalTime();
  G4double weight = aStep->GetTrack()->GetWeight();
  G4ThreeVector position = aStep->GetTrack()->GetPosition();
  G4double PosX = position.getX();
  G4double PosY = position.getY();
  G4double PosZ = position.getZ();
  G4int trackID = aStep->GetTrack()->GetTrackID();
  G4int parentID = aStep->GetTrack()->GetParentID();
  G4ThreeVector momentumdirection = aStep->GetTrack()->GetMomentumDirection();
  G4double MomX = momentumdirection.getX();
  G4double MomY = momentumdirection.getY();
  G4double MomZ = momentumdirection.getZ();
  G4double velocity = aStep->GetTrack()->GetVelocity();

  //run->ParticleCount(name,energy,iVol);                                                                    
  //G4int id =0;
  analysisManager->FillNtupleDColumn(0, double(pid));
  analysisManager->FillNtupleDColumn(1, double(Z));
  analysisManager->FillNtupleDColumn(2, double(A));
  analysisManager->FillNtupleDColumn(3, charge);
  analysisManager->FillNtupleDColumn(4, energy);
  analysisManager->FillNtupleDColumn(5, edep);
  analysisManager->FillNtupleDColumn(6, time);
  analysisManager->FillNtupleDColumn(7, weight);
  analysisManager->FillNtupleDColumn(8, PosX);
  analysisManager->FillNtupleDColumn(9, PosY);
  analysisManager->FillNtupleDColumn(10,PosZ);
  analysisManager->FillNtupleIColumn(11,trackID);
  analysisManager->FillNtupleIColumn(12,parentID);
  analysisManager->FillNtupleDColumn(13,velocity);
  analysisManager->FillNtupleDColumn(14,MomX);
  analysisManager->FillNtupleDColumn(15,MomY);
  analysisManager->FillNtupleDColumn(16,MomZ);
  analysisManager->FillNtupleIColumn(17,eventID);
  analysisManager->AddNtupleRow();

  
  //G4cout << "Volume Name:" << volumename.c_str() << "; Particle:" <<particleType.c_str() << "," << particleName.c_str() << "; Energy Deposit:" << edep << "; Position:"<< position << "; Time:" << p_time << G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
