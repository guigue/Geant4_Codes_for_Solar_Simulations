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


 // Guiding Center equation of motion. The guiding center code
 //depends of x0, y0, z0, vp0, v, dt, 'the field' (writed in other module)
 //particle mass (implicit in other module) and the value of the charge q.

 //In Geant's User Guide Developer the advising is start writing a new equation
 //of motion from the class G4Mag_EqRhs.

#ifndef G4MAG_GC
#define G4MAG_GC

#include "G4Mag_EqRhs.hh"
#include "G4ChargeState.hh"
//#include "Solar_TabulatedField.hh"

class Solar_Dipole;

class Solar_GC : public G4Mag_EqRhs
{
    public:

        Solar_GC( Solar_Dipole* MagField);
        virtual ~Solar_GC();

        // Constructor and Destructor. No actions.

    void EvaluateRhsGivenB( const G4double y[],
                            const G4double B[],
                                  G4double dydx[] ) const;

        // Given the value of Magnetic Field B, this function
        // calculates the value of the derivate dydx.

    virtual void SetChargeMomentumMass ( G4ChargeState particleCharge,
                                         G4double MomentumXc,
                                         G4double mass);

    void GetFieldValue_forGrads(const double point[4],
              double *Bfield) const;


    private:
    //G4ThreeVector p_par;
    //G4ThreeVector p_perp;
    //G4ThreeVector xi;

    //G4double p_par_magnitude;
    //G4double p_perp_magnitude;

    G4double charge, mass;
    G4double beta, gamma;
    G4double MDM;

    //G4ThreeVector xi;
    //G4ThreeVector zeta;

    //G4ThreeVector deltaXi;
    //G4ThreeVector epsilonZeta;


};

#endif
