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
/// \file exoticphysics/monopole/src/PrimaryGeneratorAction.cc
/// \brief Implementation of the PrimaryGeneratorAction class
//
// $Id$
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleDefinition.hh"
#include "G4Event.hh"
#include "G4SystemOfUnits.hh"
#include "g4root.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction(DetectorConstruction* det)
  :G4VUserPrimaryGeneratorAction(),fParticleGun(0),fDetector(det),
   bPrimPositionDefined(false)
{
  fParticleGun  = new G4ParticleGun(1);
  fParticleGun->SetParticleEnergy(100*GeV);
  //fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., -1.));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  /*  
  //this function is called at the begining of event
  //  if(0 == anEvent->GetEventID() || !bPrimPositionDefined) {
  G4double z0 = 0.5*(fDetector->GetWorldSizeZ()) + 1*um;
  G4double x0 = 2*(G4UniformRand()-0.5)*m;
  G4double y0 = 2*(G4UniformRand()-0.5)*m;
  G4cout << "X Y " << x0 << " " << y0 << G4endl << G4endl;
  fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
  bPrimPositionDefined = true;
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., -1.));   
  //  }  
 
  */
  
  G4double radius = 0.5*(fDetector->GetWorldSizeZ()) + 1*um;
  G4double theta = 0.5*(G4UniformRand())*3.1415926;
  G4double phi = 2.*(G4UniformRand())*3.1415926;

  G4double z0 = radius*cos(theta);
  G4double x0 = radius*sin(theta)*cos(phi);
  G4double y0 = radius*sin(theta)*sin(phi);

  fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
  bPrimPositionDefined = true;
  //fParticleGun->SetParticleMomentumDirection(G4ThreeVector(-x0/radius,-y0/radius,-z0/radius)); 
  
  //G4double rho0 = 5*m;
  //G4double phi0 = 2.*(G4UniformRand())*3.1415926;
  
  //G4double z1 = -5*m;
  //G4double x1 = rho0*cos(phi0);
  //G4double y1 = rho0*sin(phi0);

  //G4double x1 = (G4UniformRand()-0.5)*2*m;
  //G4double y1 = (G4UniformRand()-0.5)*2*m;
  //G4double z1 = (G4UniformRand()-0.5)*2*m;

  G4double x1 = (G4UniformRand()-0.5)*1*m;
  G4double y1 = (G4UniformRand()-0.5)*1*m;
  G4double z1 = -G4UniformRand()*1*m;    
  
  //G4double x1 = 0*m;
  //G4double y1 = 0*m;
  //G4double z1 = 0*m;
  G4double vx = x1-x0;
  G4double vy = y1-y0;
  G4double vz = z1-z0;
  G4double v = sqrt(pow(vx,2)+pow(vy,2)+pow(vz,2));
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(vx/v,vy/v,vz/v));     
  
  /*
  // particle must go through the floor of the hemisphere
  double radius = 30.0*m;
  double theLowerX = radius-2*radius*G4UniformRand();
  double theLowerY = radius-2*radius*G4UniformRand();
  double theLowerZ = -5.0*m;
  double outsideBound = std::pow(theLowerX,2) + std::pow(theLowerY,2); //eq for the circle inside box
  while(outsideBound > std::pow(radius,2)) { //if the point falls outside the circle, select a new one                               
    theLowerX = radius-2*radius*G4UniformRand();
    theLowerY = radius-2*radius*G4UniformRand();
    outsideBound = std::pow(theLowerX,2) + std::pow(theLowerY,2);
  }

  // create particle on surface of hemisphere of size radius, that travels through the bottom surface                                  // pick bottom point inside a circle                                                                                                 ////starting position on the hemisphere                                                                                             
  double theUpperPhi = 3.14159-2*3.1415*G4UniformRand();
  double theUpperTheta = acos(G4UniformRand());
  double theUpperX = radius * cos(theUpperPhi) * sin(theUpperTheta);
  double theUpperY = radius * sin(theUpperPhi) * sin(theUpperTheta);
  double theUpperZ = radius * cos(theUpperTheta)+theLowerZ;
  //G4cout << "UpperX= " << theUpperX << " " << "UpperY= " << theUpperY << " "<< "UpperZ= " << theUpperZ << " " << G4endl;            
  fParticleGun->SetParticlePosition(G4ThreeVector(theUpperX, theUpperY, theUpperZ));

  //compute track direction (note left handed coord system)                                                                           
  double Xtrk = theLowerX-theUpperX;
  double Ytrk = theLowerY-theUpperY;
  double Ztrk = theLowerZ-theUpperZ;
  double Rtrk = std::sqrt(Xtrk*Xtrk + Ytrk*Ytrk + Ztrk*Ztrk);
  double ThetaTrk=acos(Ztrk/Rtrk);
  double PhiTrk=0.;
  if (Ytrk == 0){
    if (Xtrk > 0) PhiTrk = 3.14159/2.;
    if (Xtrk < 0) PhiTrk = -3.14159/2.;
  }
  else{
    PhiTrk=atan(std::abs(Xtrk)/std::abs(Ytrk));
    if ((Ytrk < 0) && (Xtrk > 0)) PhiTrk = 3.14159/2.+ PhiTrk;
    if ((Ytrk < 0) && (Xtrk <= 0)) PhiTrk = -3.14159/2. - PhiTrk;
    if ((Ytrk > 0) && (Xtrk < 0)) PhiTrk = - PhiTrk;
  }
  //G4cout << "ThetaTrk= " << ThetaTrk << " PhiTrk= " << PhiTrk << G4endl << G4endl;                                                  
  //set direction of particle (inward)                                                                                                
  G4ThreeVector dir(std::sin(ThetaTrk)*std::cos(PhiTrk),std::sin(ThetaTrk)*std::sin(PhiTrk),std::cos(ThetaTrk));
  fParticleGun->SetParticleMomentumDirection(dir);
  */  
  fParticleGun->GeneratePrimaryVertex(anEvent);
  fEventID = anEvent->GetEventID();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

