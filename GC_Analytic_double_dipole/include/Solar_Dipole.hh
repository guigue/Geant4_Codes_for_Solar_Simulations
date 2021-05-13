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


    G4double B0F = 0.1* tesla; // value for the magnetic foot
    G4double B0S = 0.1* tesla; // value for thr magnetic foot for the second dipole
    G4double Rs = 6.96e8*m;
    G4double Rf = 0.010*Rs;  // distance from center to magnetic foot
    G4double Rf2 = 0.010*Rs;  // distance from center to magnetic foot for the second dipole
    G4double lz = 0.0267632*Rs;  // loop top z-coordinate (meter)
    G4double lh = 0.0152632*Rs;  // loop height (meter)
    //G4double displacement = 500*km;

    G4double dzF = (lz - lh);  // dipoloe depth (meter) for First dipole
    G4double dzS = (lz - lh);  // dipoloe depth (meter) for Second dipole

    G4double dxF = Rf+(500*km); // Displacement in x axis for first dipole
    G4double dyF = 0; // Displacement in y axis for first dipole

    G4double dxS = -Rf2-(500*km); // Displacement in x axis for first dipole
    G4double dyS = 0*km; // Displacement in y axis for first dipole

};
#endif
