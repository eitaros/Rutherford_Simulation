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
/// \file B3aEventAction.cc
/// \brief Implementation of the B3aEventAction class

#include "B3aEventAction.hh"
#include "B3aRunAction.hh"
#include "B3Analysis.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"

#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4THitsMap.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B3aEventAction::B3aEventAction(B3aRunAction* runAction)
 : G4UserEventAction(), 
   fRunAction(runAction),
   fCollID_det(-1),
   fcopyNb(-1)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B3aEventAction::~B3aEventAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B3aEventAction::BeginOfEventAction(const G4Event* /*evt*/)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B3aEventAction::EndOfEventAction(const G4Event* evt )
{
  auto analysisManager = G4AnalysisManager::Instance();
   //Hits collections
  //  
  G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
  if(!HCE) return;
               
   // Get hits collections IDs
  if (fCollID_det < 0) {
    G4SDManager* SDMan = G4SDManager::GetSDMpointer();  
    fCollID_det   = SDMan->GetCollectionID("detector/edep");
  }
    //入射したらfcopyNbを引数にしてCountEventを起動
    //
    const G4double eThreshold = 0.*keV;
    G4int nbOfFired = 0;
    
    G4THitsMap<G4double>* evtMap =
                       (G4THitsMap<G4double>*)(HCE->GetHC(fCollID_det));
                 
    std::map<G4int,G4double*>::iterator itr;
    for (itr = evtMap->GetMap()->begin(); itr != evtMap->GetMap()->end(); itr++) {
      fcopyNb  = (itr->first);
      G4double edep = *(itr->second);
      if (edep > eThreshold) nbOfFired++;
      //G4cout << "det" << fcopyNb << ": " << edep/MeV << " MeV " << G4endl;
    }
    if (fcopyNb > 36)
    {
      fcopyNb = 72 - fcopyNb;
    }
    if ( nbOfFired == 1 ){fRunAction->CountEvent(fcopyNb);}
      
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
