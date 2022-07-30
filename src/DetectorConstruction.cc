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
/// \file exoticphysics/monopole/src/DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
// $Id$
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Polycone.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4UniformMagField.hh"

#include "G4GeometryManager.hh"
#include "G4UserLimits.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4UnitsTable.hh"
#include "G4NistManager.hh"

#include "G4MonopoleFieldSetup.hh"
//#include "G4FieldManager.hh"
//#include "G4TransportationManager.hh"
#include "G4ThreeVector.hh"
#include "G4RunManager.hh" 
#include "G4SystemOfUnits.hh"

#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
 : G4VUserDetectorConstruction(),
   fWorldMaterial(0),           
   fAbsorMaterial(0),
   fLogAbsor(0),
   fMonFieldSetup(0),
   fZMagFieldValue(0.),
   fDetectorMessenger(0)
{
  // default parameter values
  
  fAbsorSizeX = fAbsorSizeY = fAbsorSizeZ = 10.*cm; // dummy parameters; replaced by the real MJD detectors
  fWorldSizeX = fWorldSizeY = 3000.*m;
  fWorldSizeZ = 3000.*m;
  fMaxStepSize = 5 * mm;
  G4cout << "World size X,Y:"<< fWorldSizeX <<"mm"<< G4endl;

  G4cout << "World size Z:"<< fWorldSizeZ << "mm"<< G4endl;

  //  fMonFieldSetup = G4MonopoleFieldSetup::GetMonopoleFieldSetup();
  fMonFieldSetup = new G4MonopoleFieldSetup();

  SetMaterial("G4_Ge");
  
  G4NistManager* nist = G4NistManager::Instance();
  fWorldMaterial = nist->FindOrBuildMaterial("G4_SILICON_DIOXIDE");

  //G4Material* mat_Silicon = nist->FindOrBuildMaterial("G4_Si");
  //G4Material* mat_Oxygen = nist->FindOrBuildMaterial("G4_O");
  //G4Material* mat_Rock= new G4Material("Rock", 2.86 * g / cm3, 5);
  //G4double sum_minors = 0.22/1000000. + 0.33/1000000. + 0.96/100.;
  //mat_Rock->AddMaterial(mat_Silicon, 1./3.*(1-sum_minors));
  //mat_Rock->AddMaterial(mat_Oxygen, 2./3.*(1-sum_minors));

  //fWorldMaterial=  G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");  

  // create commands for interactive definition of the detector
  fDetectorMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
  delete fDetectorMessenger;
  //  delete fMonFieldSetup;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{ 
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  /****************************    World   *****************************/
  G4Box * sWorld = new G4Box("world",                                
               fWorldSizeX, fWorldSizeY, fWorldSizeZ);

  //G4Material* mat_Rock= new G4Material("Rock", 2.86 * g / cm3, 5);

  G4LogicalVolume * lWorld = new G4LogicalVolume(sWorld,
                                                 fWorldMaterial,
                                                 "world");        

  G4VPhysicalVolume * pWorld = new G4PVPlacement(0,      //no rotation
                             G4ThreeVector(),            //at (0,0,0)
                           lWorld,                       //logical volume
                           "world",                      //name
                           0,                            //mother  volume
                           false,                        //no boolean operation
                           0);                           //copy number


  /**************************    Absorber    ***************************/
  /*G4Box * sAbsor = new G4Box("Absorber",                   
  fAbsorSizeX / 2, fAbsorSizeY / 2, fAbsorSizeZ / 2);        

  //Ge density: 5.33 g/cm^3
  //If each detector is about a kg, volume=~188 cm^3
  //Cylinder 4.5 cm radius, 3 cm half-z has a volume of ~191 cm^3
  //Use this if the other implementation is too much
  //G4Tubs *sAbsor = new G4Tubs("GeDet",0,4.5*cm,3*cm,0,2*M_PI);

    fLogAbsor = new G4LogicalVolume(sAbsor, fAbsorMaterial, "Absorber");
  fLogAbsor->SetUserLimits(new G4UserLimits(fMaxStepSize));
  
  new G4PVPlacement(0,                                //no rotation
                    G4ThreeVector(),                //at (0,0,0)
                    fLogAbsor,                        //logical volume
                   "Absorber",                        //name
  lWorld,                               //mother  volume
  false,                        //no boolean operation
  0);                                //copy number
  */
  //PrintParameters();

  //read detector file
  using namespace std;
  vector<string> Detectors;
  string input;
  G4String s1;

  G4int Detector_number;
  G4String Detector_name;
  G4String Detector_name_sol;
  G4String Detector_name_log;
  G4ThreeVector Detector_position;
  G4int Detector_slices=0;
  vector<G4double> Detector_r;
  vector<G4double> Detector_z;
  G4RotationMatrix* Detrotation = new G4RotationMatrix();
  Detrotation->rotateX(180*deg);

  int j =0;

  ifstream in ("Detectorposition.txt");
  while(getline(in,input)){
    Detectors.push_back(input);
  }
  in.close();
  G4cout << "-------" << Detectors.size() <<G4endl;

  for (int i =0;i< (int)Detectors.size();i++){
    G4cout << i << " " << Detectors[i] << endl;
    j = 0;
    Detector_number = 0;
    Detector_r.clear();
    Detector_z.clear();
    stringstream ss(Detectors[i]);

    while (ss >> s1) {
      j++;
      if (j==1) Detector_number += atoi(s1.c_str())*100;
      else if (j==2) Detector_number += atoi(s1.c_str())*10;
      else if (j==3) Detector_number += atoi(s1.c_str())*1;
      else if (j==4) Detector_position.setX(atof(s1.c_str())*mm);
      else if (j==5) Detector_position.setY(atof(s1.c_str())*mm);
      else if (j==6) Detector_position.setZ(atof(s1.c_str())*mm);
      else if (j==7) Detector_name = "Det_" + s1;
      else if (j==8) Detector_slices = atoi(s1.c_str());
      else if (j==9) Detector_r.push_back((G4double) atof(s1.c_str())*mm);
      else if (j==10) Detector_z.push_back((G4double) atof(s1.c_str())*mm);
      else if (j==11) Detector_r.push_back((G4double) atof(s1.c_str())*mm);
      else if (j==12) Detector_z.push_back((G4double) atof(s1.c_str())*mm);
      else if (j==13) Detector_r.push_back((G4double) atof(s1.c_str())*mm);
      else if (j==14) Detector_z.push_back((G4double) atof(s1.c_str())*mm);
    }

    if(Detector_slices==2){
      Detector_r.push_back(Detector_r[1]+0.00001*mm);
      Detector_z.push_back(Detector_z[1]+0.00001*mm);
    }

    Detector_position.setZ(Detector_position.z() + Detector_z[2] - 20*cm);

    G4cout << j << " --- " << Detector_name << " " << Detector_number << " " << Detector_position <<  G4endl;

    Detector_name_sol = Detector_name + "_sol";
    Detector_name_log = Detector_name + "_log";
    const G4double r_i[] = {0,0,0};
    const G4double r[] = {Detector_r[0],Detector_r[1],Detector_r[2]};
    const G4double z[] = {Detector_z[0],Detector_z[1],Detector_z[2]};

    G4Polycone* Det_solid = new G4Polycone(Detector_name_sol,0,2*M_PI,3,z,r_i,r);
    //G4LogicalVolume* 
    fLogAbsor = new G4LogicalVolume(Det_solid,fAbsorMaterial,Detector_name_log);
    fLogAbsor->SetUserLimits(new G4UserLimits(fMaxStepSize));

    new G4PVPlacement (Detrotation,Detector_position,fLogAbsor,Detector_name,lWorld,false,Detector_number,0);
  }




  /************     always return the World volume     *****************/
  return pWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PrintParameters()
{
  G4cout << "\n---------------------------------------------------------\n";
  //G4cout << "---> The Absorber is " << G4BestUnit(fAbsorSizeX, "Length")
  //       << " of " << fAbsorMaterial->GetName() << G4endl;
  G4cout << "\n---------------------------------------------------------\n";

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetSizeX(G4double value)
{
  if(value > 0.0) {
    fAbsorSizeX = value; 
    //fWorldSizeX = 1.2 * fAbsorSizeX;
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
  }
}
  
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetSizeY(G4double value)
{
  if(value > 0.0) {
    fAbsorSizeY = value;
    //fWorldSizeY = 1.2 * fAbsorSizeY;
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
  }
}  

void DetectorConstruction::SetSizeZ(G4double value)
{
  if(value > 0.0) {
    fAbsorSizeZ = value;
    //fWorldSizeZ = 1.2 *fAbsorSizeZ;
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
  }
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetMaterial(const G4String& namemat)
{
  // search the material by its name   
  G4Material* mat = G4NistManager::Instance()->FindOrBuildMaterial(namemat);
  if(!mat) {
    G4cout << "!!! DetectorConstruction::SetMaterial: WARNING Material <"
           << namemat << "> does not exist in DB" << G4endl;
    return;
  }
  // new material is found out
  if (mat != fAbsorMaterial) {
    fAbsorMaterial = mat;
    if(fLogAbsor) { fLogAbsor->SetMaterial(mat); }
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


// void DetectorConstruction::SetMagField(G4double fieldValue)
// {
//   fMonFieldSetup->SetMagField(fieldValue);

//   //apply a global uniform magnetic field along Z axis
//   G4FieldManager * fieldMgr = 
//     G4TransportationManager::GetTransportationManager()->GetFieldManager();
    
//   if (fMagField) { delete fMagField; }    //delete the existing magn field

//   if (fieldValue != 0.)                   // create a new one if non nul
//     {
//       fMagField = new G4UniformMagField(G4ThreeVector(0., 0., fieldValue));
//       fieldMgr->SetDetectorField(fMagField);
//       fieldMgr->CreateChordFinder(fMagField);
//     }
//    else
//     {
//       fMagField = 0;
//       fieldMgr->SetDetectorField(fMagField);
//     }
// }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void DetectorConstruction::ConstructSDandField()
{

  // Define magnetic field
  bool bNewFieldValue = false;
  if ( fFieldMessenger.Get() != 0 ) {
    G4ThreeVector fieldSet =  fFieldMessenger.Get()->GetFieldValue();
    if(fieldSet.z()!=fZMagFieldValue) bNewFieldValue = true;
  }
  else bNewFieldValue = true;

  // Monopole particule specific magnetic field
  if(bNewFieldValue&&fZMagFieldValue!=0.)
    fMonFieldSetup->SetMagField(fZMagFieldValue, true);

  if ( bNewFieldValue ) { 
    // Create global magnetic field messenger.
    // Uniform magnetic field is then created automatically if
    // the field value is not zero.

    if(fZMagFieldValue!=0.)
      {
        G4ThreeVector fieldValue = G4ThreeVector(0.,0.,fZMagFieldValue);
        G4GlobalMagFieldMessenger* msg =  
                              new G4GlobalMagFieldMessenger(fieldValue);
        msg->SetVerboseLevel(1);
        G4AutoDelete::Register(msg);
        fFieldMessenger.Put( msg );        
      }
  }
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetMaxStepSize(G4double step)
{
  fMaxStepSize = step;
  if(fLogAbsor) { fLogAbsor->SetUserLimits(new G4UserLimits(fMaxStepSize)); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
