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
/// \file B3PrimaryGeneratorAction.cc
/// \brief Implementation of the B3PrimaryGeneratorAction class

#include "B3PrimaryGeneratorAction.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ChargedGeantino.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B3PrimaryGeneratorAction::B3PrimaryGeneratorAction()
 : G4VUserPrimaryGeneratorAction(),
   fParticleGun(0)
{
  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);

  // default particle kinematic

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* Particle
                    = particleTable->FindParticle("alpha");
  fParticleGun->SetParticleDefinition(Particle);
  fParticleGun->SetParticlePosition(G4ThreeVector(0.,-1.,0.));
  fParticleGun->SetParticleEnergy(5.4*MeV);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,1.,0.));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B3PrimaryGeneratorAction::~B3PrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B3PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  G4ParticleDefinition* particle = fParticleGun->GetParticleDefinition();
  

  // randomized position
  //
  ///G4double x0  = 0*cm, y0  = 0*cm, z0  = 0*cm;
  ///G4double dx0 = 0*cm, dy0 = 0*cm, dz0 = 0*cm;
  
  G4double x0  = 0*cm, y0  = -1*cm, z0  = 0*cm;
  G4double dx0 = 1*mm, dy0 = 0*mm, dz0 = 1*mm;
  /*
  x0 += dx0*(G4UniformRand()-0.5);
  y0 += dy0*(G4UniformRand()-0.5);
  z0 += dz0*(G4UniformRand()-0.5);
  */
  fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
  
   
    // randomized momentum
    //
    
    G4double coli_r = 1*mm;
    G4double y = 10*mm;
    G4double theta = 4*CLHEP::pi*(G4UniformRand()-0.5);
    G4double x = std::cos(theta)*coli_r*2*(G4UniformRand()-0.5);
    G4double z = std::sin(theta)*coli_r*2*(G4UniformRand()-0.5);
    
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(x,y,z));
   
    
    
  //create vertex
  //
  fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

