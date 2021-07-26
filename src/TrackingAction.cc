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
/// \file exoticphysics/monopole/src/TrackingAction.cc
/// \brief Implementation of the TrackingAction class
//
// $Id$
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "TrackingAction.hh"
#include "G4RunManager.hh"
#include "Run.hh"

#include "G4Track.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackingAction::TrackingAction()
 : G4UserTrackingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackingAction::~TrackingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TrackingAction::PreUserTrackingAction(const G4Track* track)
{
  G4cout << track->GetTrackID() << G4endl;
  //Run* run = static_cast<Run*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());

  //G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  //which volume ?                         
  //G4LogicalVolume* lVolume = track->GetVolume()->GetLogicalVolume();
  /*
  const G4ParticleDefinition* particle = track->GetParticleDefinition();
  G4String name   = particle->GetParticleName();
  G4int pid       = particle->GetPDGEncoding();
  G4int Z         = particle->GetAtomicNumber();
  G4int A         = particle->GetAtomicMass();
  G4double charge = particle->GetPDGCharge();
  G4double energy = track->GetKineticEnergy();
  G4double time   = track->GetGlobalTime();
  G4double weight = track->GetWeight();
  G4ThreeVector position = track->GetPosition();
  G4double PosX = position.getX();
  G4double PosY = position.getY();
  G4double PosZ = position.getZ(); 
  //run->ParticleCount(name,energy,iVol);
  G4int id =0;
  analysisManager->FillNtupleDColumn(id,0, double(pid));
  analysisManager->FillNtupleDColumn(id,1, double(Z));
  analysisManager->FillNtupleDColumn(id,2, double(A));
  analysisManager->FillNtupleDColumn(id,3, energy);
  analysisManager->FillNtupleDColumn(id,4, time);
  analysisManager->FillNtupleDColumn(id,5, weight);
  analysisManager->FillNtupleDColumn(id,6, PosX);
  analysisManager->FillNtupleDColumn(id,7, PosY);
  analysisManager->FillNtupleDColumn(id,8,PosZ);
  analysisManager->FillNtupleDColumn(id,9,charge);
  analysisManager->AddNtupleRow(id);
  */
}

void TrackingAction::PostUserTrackingAction(const G4Track* aTrack)
{
  G4cout << aTrack->GetTrackID() << G4endl;
  // extract Projected Range of primary particle
  //if (aTrack->GetTrackID() == 1) {
  //  Run* run = static_cast<Run*>(
  //            G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  //  G4double x = aTrack->GetPosition().x() - run->GetOffsetX();
  //  run->AddProjRange(x);
  // }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

