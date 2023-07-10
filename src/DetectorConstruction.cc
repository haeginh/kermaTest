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
/// \file B1/src/DetectorConstruction.cc
/// \brief Implementation of the B1::DetectorConstruction class

#include "DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4PSEnergyDeposit.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();

  // Envelope parameters
  //
  G4Material* water = nist->FindOrBuildMaterial("G4_WATER");
  G4Material* tissue = nist->FindOrBuildMaterial("G4_TISSUE_SOFT_ICRU-4");
  G4Material* air = nist->FindOrBuildMaterial("G4_AIR");

  //
  // World
  //
  
  auto solidWorld = new G4Box("World", 5*m, 5*m, 5*m);                              
  auto logicWorld = new G4LogicalVolume(solidWorld, air, "world");
  auto physWorld = new G4PVPlacement(nullptr,  // no rotation
    G4ThreeVector(),                           // at (0,0,0)
    logicWorld,                                // its logical volume
    "World",                                   // its name
    nullptr,                                   // its mother  volume
    false,                                     // no boolean operation
    0,                                         // copy number
    false);                            // overlaps checking

  //
  // tissue
  //
  auto solidThin = new G4Tubs("Thin", 0, 10*cm, 1*mm, 0, 2*M_PI);  
  auto solidThick = new G4Tubs("Thick", 0, 10*cm, 0.5*cm, 0, 2*M_PI);  
  fScoringVolume = new G4LogicalVolume(solidThin, tissue, "Thin");
  auto logicThick = new G4LogicalVolume(solidThick, tissue, "Thick");
  new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 1.1*cm), fScoringVolume, "Thin", logicWorld, false, 0, false);
  new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 0.5*cm), logicThick, "Thick", logicWorld, false, 0, false);
  
  return physWorld;
}

void DetectorConstruction::ConstructSDandField()
{
  G4SDManager* sdMan = G4SDManager::GetSDMpointer();
  G4MultiFunctionalDetector* mfd = new G4MultiFunctionalDetector("det");
  sdMan->AddNewDetector(mfd);
  mfd->RegisterPrimitive(new G4PSEnergyDeposit("edep"));
  SetSensitiveDetector(fScoringVolume, mfd);
}

