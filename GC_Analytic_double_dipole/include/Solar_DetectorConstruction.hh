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

#ifndef Solar_DetectorConstruction_h
#define Solar_DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4Cache.hh"
#include "G4ThreeVector.hh"


class G4VPhysicalVolume;
class G4LogicalVolume;

class G4Material;
class G4UserLimits;
class G4GlobalMagFieldMessenger;

class G4Box;
class G4Region;

//class Solar_CalorimeterSD;
class Solar_FieldSetup;

/// Detector construction class to define materials and geometry.

class Solar_DetectorConstruction : public G4VUserDetectorConstruction
{
  public:

    Solar_DetectorConstruction();
    virtual ~Solar_DetectorConstruction();

  public:

    virtual G4VPhysicalVolume* Construct();
    // Gives the magnetic field B at a given position defined by yTrack
    //void GetFieldValue(const G4double yTrack[], G4double B[]) const;
    //G4ThreeVector GetFieldValue(const G4ThreeVector position) const;

    virtual void ConstructSDandField();

    G4double newBigDist;

    G4bool  fCheckOverlaps; // option to activate checking of volumes overlaps

    G4Box* WorldBox;
    G4Box* CoronaBox;
    G4Box* AtmBox;

    G4LogicalVolume*WorldLog;
    G4LogicalVolume*AtmLog;
    G4LogicalVolume*CoronaLog;

    G4VPhysicalVolume*WorldPhys;
    G4VPhysicalVolume*AtmPhys;
    G4VPhysicalVolume*CoronaPhys;


    G4int ncomponents;
    G4double fractionmass, density, densityCorona;

    G4double x, y, z;

    G4Region* atmRegion;
    G4Region* corRegion;
    G4Region* sphRegion;

  private:

    G4LogicalVolume* fMagneticLogical;
    G4Cache<Solar_FieldSetup*> fEmFieldSetup;
    G4UserLimits* SolarUserLimits;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
