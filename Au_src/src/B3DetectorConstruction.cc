//
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
//
/// \file B3DetectorConstruction.cc
/// \brief Implementation of the B3DetectorConstruction class

#include "B3DetectorConstruction.hh"

#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PSDoseDeposit.hh"
#include "G4VisAttributes.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
//oOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B3DetectorConstruction::B3DetectorConstruction()
: G4VUserDetectorConstruction(),
  fCheckOverlaps(true)
{
  DefineMaterials();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B3DetectorConstruction::~B3DetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B3DetectorConstruction::DefineMaterials()
{ }
G4Material* get_material(const std::string& name) {
  
  auto* nist = G4NistManager::Instance();
  auto* mat = nist->FindOrBuildMaterial(name);
  if (mat!=0) return mat;

  
  std::vector<int>    natom;
  std::vector<double> wfrac;
    if (name=="Vacuum"){
      mat = nist->ConstructNewMaterial(name, {"N","O"}, wfrac={0.7,0.3}, 1e-8*gram/cm3);
    }
    return mat;
  }
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B3DetectorConstruction::Construct()
{
  //  detector Parameters
  //
  G4double det_dx = 1.*cm, det_dy = 0.5*cm, det_dz = 4.*cm;
  G4double theta = 5.;
  G4double r = 30.*cm;
    // target parameters
    //
  G4double target_xy = 10.*mm;
  G4double target_z = 0.0025*mm;
    //
    //
    G4double acryl_r = 35.*cm;
    G4double acryl_sizeZ  = 20.*cm;
    //
  G4NistManager* nist = G4NistManager::Instance();
  G4Material* target_mat1 = nist->FindOrBuildMaterial("G4_Au");
  G4Material* target_mat2 = nist->FindOrBuildMaterial("G4_Al");
  G4Material* target_mat3 = nist->FindOrBuildMaterial("G4_Cu");
  G4Material* det_mat   = nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
  //
  // World
  //
  G4double world_size = 1.2*acryl_r;
  
  G4Box* solidWorld =
    new G4Box("World",                       //its name
       0.5*world_size, 0.5*world_size, 0.5*world_size); //its size
      
  G4LogicalVolume* logicWorld =
    new G4LogicalVolume(solidWorld,          //its solid
                        get_material("Vacuum"),         //its material
                        "World");            //its name
                                   
  G4VPhysicalVolume* physWorld =
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      fCheckOverlaps);       // checking overlaps
  //
  //acryl tub
  //
    G4Tubs* solidacryl =
      new G4Tubs("acryl",                       //its name
         0., 0.5*acryl_r, 0.5*acryl_sizeZ, 0., twopi); //its size
        
    G4LogicalVolume* logicacryl =
      new G4LogicalVolume(solidacryl,          //its solid
                          get_material("Vacuum"),         //its material
                          "World");            //its name
                                     
    G4VPhysicalVolume* physacryl =
      new G4PVPlacement(0,                     //no rotation
                        G4ThreeVector(),       //at (0,0,0)
                        logicacryl,            //its logical volume
                        "World",               //its name
                        logicWorld,                     //its mother  volume
                        false,                 //no boolean operation
                        0,                     //copy number
                        fCheckOverlaps);       // checking overlaps
  //
  // target
  //
    G4Box* solidtarget =
    new G4Box("target",
              0.5*target_xy, 0.5*target_z, 0.5*target_xy);
      
  G4LogicalVolume* logictarget =
    new G4LogicalVolume(solidtarget,           //its solid
                        target_mat1,         //its material
                        "target");             //its name
   
    G4VPhysicalVolume* phystarget =
      new G4PVPlacement(0,                     //no rotation
                        G4ThreeVector(),       //at (0,0,0)
                        logictarget,            //its logical volume
                        "target",               //its name
                        logicacryl,                     //its mother  volume
                        false,                 //no boolean operation
                        0,                     //copy number
                        fCheckOverlaps);       // checking overlaps
                 
  //
  // detectors
  //
  G4Box* soliddet =
    new G4Box("detector",
              det_dx/2, det_dy/2, det_dz/2);
 
    G4LogicalVolume* logicdet =
        new G4LogicalVolume(soliddet,
                            det_mat,
                            "detectorLV");
    
G4double maxID = 360/theta  ;
  for (G4int detID = 0; detID <= maxID ; detID++) 
  {
    G4double phi1 = detID*theta*twopi/360;
    G4RotationMatrix rotm1  = G4RotationMatrix();
    rotm1.rotateZ(-1*phi1);
    G4ThreeVector uz = G4ThreeVector(std::sin(phi1), std::cos(phi1),  0.);
    G4ThreeVector position1 = (0.5*r + 0.5*det_dy)*uz;
    G4Transform3D transform1 = G4Transform3D(rotm1,position1);
    
    G4VPhysicalVolume* physdet =
     new G4PVPlacement(transform1,
                       logicdet,
                       "detector",
                       logicacryl,
                       false,
                       detID,
                       fCheckOverlaps);
 
  }

/*
G4double maxID = 3;
for (G4int detID = 0; detID <= maxID ; detID++) 
  {

    if (detID == 0)
    {
      G4double phi1 = 0*twopi/360;
      G4RotationMatrix rotm1  = G4RotationMatrix();
      rotm1.rotateZ(phi1);
      G4ThreeVector uz = G4ThreeVector(std::sin(phi1), std::cos(phi1),  0.);
      G4ThreeVector position1 = (0.5*r + 0.5*det_dy)*uz;
      G4Transform3D transform1 = G4Transform3D(rotm1,position1);

      G4VPhysicalVolume* physdet =
     new G4PVPlacement(transform1,
                       logicdet,
                       "detector",
                       logicacryl,
                       false,
                       detID,
                       fCheckOverlaps);
    }

    if (detID == 1)
    {
      G4double phi1 = 10*twopi/360;
      G4RotationMatrix rotm1  = G4RotationMatrix();
      rotm1.rotateZ(phi1);
      G4ThreeVector uz = G4ThreeVector(std::sin(-1*phi1), std::cos(-1*phi1),  0.);
      G4ThreeVector position1 = (0.5*r + 0.5*det_dy)*uz;
      G4Transform3D transform1 = G4Transform3D(rotm1,position1);

      G4VPhysicalVolume* physdet =
     new G4PVPlacement(transform1,
                       logicdet,
                       "detector",
                       logicacryl,
                       false,
                       detID,
                       fCheckOverlaps);
    }

    if (detID == 2)
    {
      G4double phi1 = 15*twopi/360;
      G4RotationMatrix rotm1  = G4RotationMatrix();
      rotm1.rotateZ(-1*phi1);
      G4ThreeVector uz = G4ThreeVector(std::sin(phi1), std::cos(phi1),  0.);
      G4ThreeVector position1 = (0.5*r + 0.5*det_dy)*uz;
      G4Transform3D transform1 = G4Transform3D(rotm1,position1);

      G4VPhysicalVolume* physdet =
     new G4PVPlacement(transform1,
                       logicdet,
                       "detector",
                       logicacryl,
                       false,
                       detID,
                       fCheckOverlaps);
    }
 
  }
  */
  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl; 

  //always return the physical World
  //
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

     void B3DetectorConstruction::ConstructSDandField()
     {
      G4SDManager::GetSDMpointer()->SetVerboseLevel(1);
    
         G4MultiFunctionalDetector* det = new G4MultiFunctionalDetector("detector");
         G4SDManager::GetSDMpointer()->AddNewDetector(det);
         G4VPrimitiveScorer* primitiv1 = new G4PSEnergyDeposit("edep");
         det->RegisterPrimitive(primitiv1);
         SetSensitiveDetector("detectorLV",det);
     }


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
