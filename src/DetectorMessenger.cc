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
/// \file exoticphysics/monopole/src/DetectorMessenger.cc
/// \brief Implementation of the DetectorMessenger class
//
// $Id$
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorMessenger.hh"

#include "DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger(DetectorConstruction * det)
  :G4UImessenger(),
   fDetector(det),
   fDetDir(0),
   fMaterCmd(0),
   fSizeXCmd(0),
   fSizeYCmd(0),
   fSizeZCmd(0),
   fStepSizeCmd(0),
   fMagFieldCmd(0)
    
{ 
  fDetDir = new G4UIdirectory("/testex/det/");
  fDetDir->SetGuidance("detector construction commands");
      
  fMaterCmd = new G4UIcmdWithAString("/testex/det/setMat",this);
  fMaterCmd->SetGuidance("Select material of the box.");
  fMaterCmd->SetParameterName("choice",false);
  fMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fSizeXCmd = new G4UIcmdWithADoubleAndUnit("/testex/det/setSizeX",this);
  fSizeXCmd->SetGuidance("Set sizeX of the absorber");
  fSizeXCmd->SetParameterName("SizeX",false);
  fSizeXCmd->SetRange("SizeX>0.");
  fSizeXCmd->SetUnitCategory("Length");
  fSizeXCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fSizeYCmd = new G4UIcmdWithADoubleAndUnit("/testex/det/setSizeY",this);
  fSizeYCmd->SetGuidance("Set sizeY of the absorber");
  fSizeYCmd->SetParameterName("SizeY",false);
  fSizeYCmd->SetRange("SizeY>0.");
  fSizeYCmd->SetUnitCategory("Length");
  fSizeYCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fSizeZCmd = new G4UIcmdWithADoubleAndUnit("/testex/det/setSizeZ",this);
  fSizeZCmd->SetGuidance("Set sizeZ of the absorber");
  fSizeZCmd->SetParameterName("SizeZ",false);
  fSizeZCmd->SetRange("SizeZ>0.");
  fSizeZCmd->SetUnitCategory("Length");
  fSizeZCmd->AvailableForStates(G4State_PreInit,G4State_Idle);


  fStepSizeCmd = new G4UIcmdWithADoubleAndUnit("/testex/det/setStepSize",this);
  fStepSizeCmd->SetGuidance("Set maxStepSize in the absorber");
  fStepSizeCmd->SetParameterName("StepSize",false);
  fStepSizeCmd->SetRange("StepSize>0.");
  fStepSizeCmd->SetUnitCategory("Length");
  fStepSizeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);


  fMagFieldCmd = new G4UIcmdWithADoubleAndUnit("/testex/det/setField",this);
  fMagFieldCmd->SetGuidance("Define magnetic field.");
  fMagFieldCmd->SetGuidance("Magnetic field will be in Z direction.");
  fMagFieldCmd->SetParameterName("Bz",false);
  fMagFieldCmd->SetUnitCategory("Magnetic flux density");
  fMagFieldCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fUpdateCmd = new G4UIcmdWithoutParameter("/testex/det/update",this);
  fUpdateCmd->SetGuidance("Update calorimeter geometry.");
  fUpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  fUpdateCmd->SetGuidance("if you changed geometrical value(s).");
  fUpdateCmd->AvailableForStates(G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
  delete fMaterCmd;
  delete fSizeXCmd;
  delete fSizeYCmd;
  delete fSizeZCmd;
  delete fStepSizeCmd;
  delete fMagFieldCmd;
  delete fUpdateCmd;
  delete fDetDir;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if( command == fMaterCmd )
   { fDetector->SetMaterial(newValue);}
   
  if( command == fSizeXCmd )
   { fDetector->SetSizeX(fSizeXCmd->GetNewDoubleValue(newValue));}
   
  if( command == fSizeYCmd )
   { fDetector->SetSizeY(fSizeYCmd->GetNewDoubleValue(newValue));}
      
  if( command == fSizeZCmd )
    { fDetector->SetSizeZ(fSizeZCmd->GetNewDoubleValue(newValue));}

  if( command == fMagFieldCmd )
   { fDetector->SetMagField(fMagFieldCmd->GetNewDoubleValue(newValue));}

  if( command == fStepSizeCmd )
   { fDetector->SetMaxStepSize(fStepSizeCmd->GetNewDoubleValue(newValue));}
              
  if( command == fUpdateCmd )
   { fDetector->UpdateGeometry();}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
