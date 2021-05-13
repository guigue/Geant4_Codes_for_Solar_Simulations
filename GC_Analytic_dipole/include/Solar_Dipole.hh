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


// class Solar_Dipole
//
// Class description:
//
// Class to implement a dipole magnetic field

// -------------------------------------------------------------------

#ifndef SOLAR_DIPOLE_HH
#define SOLAR_DIPOLE_HH

#include "G4MagneticField.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4SystemOfUnits.hh"

class Solar_Dipole : public G4MagneticField
{
  public: // with description

    Solar_Dipole();

   ~Solar_Dipole();

    void GetFieldValue(const G4double yTrack[],
                            G4double B[]) const;
    G4Field* Clone() const;

  private:

    G4double B0 = 0.1* tesla;
    G4double Rs = 6.96e8*m;
    G4double Rf = 0.010*Rs;
    G4double lz = 0.0267632*Rs;  // loop top z-coordinate (meter)
    G4double lh = 0.0152632*Rs;  // loop height (meter)
    G4double dz = (lz - lh);  // dipoloe depth (meter)

};
#endif
