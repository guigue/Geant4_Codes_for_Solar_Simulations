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

#include "Solar_Dipole.hh"
#include "G4RotationMatrix.hh"
#include "G4SystemOfUnits.hh"

static G4RotationMatrix IdentityMatrix;

Solar_Dipole::Solar_Dipole()
{
}


G4Field* Solar_Dipole::Clone() const
{
    return new Solar_Dipole(/*fGradient, fOrigin, fpMatrix*/);
}

/////////////////////////////////////////////////////////////////////////

Solar_Dipole::~Solar_Dipole()
{
}

// Should decalere function as G4ThreeVector
void Solar_Dipole::GetFieldValue(const G4double y[7],
                                       G4double B[12]) const
{



  // Magnetic Momentum

  G4double Mf = B0*pow(Rf*Rf+dz*dz,5.0/2.0)/sqrt(4.0*Rf*Rf*Rf*Rf + 5.0*(Rf*Rf)*(dz*dz) + dz*dz*dz*dz);

  //G4double fac =  Mf/pow(x*x+y*y+z*z, 5.0/2.0);

  G4double fac =  Mf/pow(y[0]*y[0]+y[1]*y[1]+(y[2]+dz)*(y[2]+dz), 5.0/2.0);
  G4double fac2 = Mf/pow(y[0]*y[0] + y[1]*y[1] + (y[2] + dz)*(y[2] + dz),7.0/2.0);

  //G4double Bx1 = (2.0*y[0]*y[0] -y[1]*y[1]- (y[2]+dz)*(y[2]+dz))*fac;
  //G4double By1 = 3.0*y[0]*y[1]*fac;
  //G4double Bz1 = 3.0*y[0]*(y[2]+dz)*fac;

  //G4double Bx2 = ((2.0*(y[0]-x0)*(y[0]-x0)) -((y[1]-y0)*(y[1]-y0))- (y[2]+dz)*(y[2]+dz))*fac;
  //G4double By2 = 3.0*(y[0]-x0)*(y[1]-y0)*fac;
  //G4double Bz2 = 3.0*(y[0]-x0)*(y[2]+dz)*fac;

  //Bx, By, Bz;
  B[0] = (2.0*y[0]*y[0] -y[1]*y[1]- (y[2]+dz)*(y[2]+dz))*fac;
  B[1] = 3.0*y[0]*y[1]*fac;
  B[2] = 3.0*y[0]*(y[2]+dz)*fac;
  // Gradients
  B[3] = -3.0*y[0]*(-3.0*dz*dz + 2.0*y[0]*y[0] - 6.0*dz*y[2] - 3.0*(y[1]*y[1] + y[2]*y[2]))*fac2; //dBxdx
  B[4] =  3.0*y[1]*(dz*dz - 4.0*y[0]*y[0] + y[1]*y[1] + 2.0*dz*y[2] + y[2]*y[2])*fac2; //dBxdy
  B[5] =  3.0*(y[2] + dz)*(dz*dz- 4.0*y[0]*y[0] + y[1]*y[1] + 2.0*dz*y[2] + y[2]*y[2])*fac2; //dBxdz
  B[6] =  3.0*y[1]*(dz*dz - 4.0*y[0]*y[0] + y[1]*y[1] + 2.0*dz*y[2] + y[2]*y[2])*fac2; //dBydx
  B[7] =  3.0*y[0]*(dz*dz + y[0]*y[0] -4.0*y[1]*y[1] + 2.0*dz*y[2] + y[2]*y[2])*fac2; //dBydy
  B[8] = -15.0*y[0]*y[1]*(y[2] + dz)*fac2; //dBydz
  B[9] =  3.0*(y[2] + dz)*(dz*dz - 4.0*y[0]*y[0] + y[1]*y[1] + 2.0*dz*y[2] + y[2]*y[2])*fac2; //dBzdx
  B[10] = -15.0*y[0]*y[1]*(y[2] + dz)*fac2; //dBzdy
  B[11] = 3.0*y[0]*(-4.0*dz*dz + y[0]*y[0] + y[1]*y[1] - 8.0*dz*y[2] - 4.0*y[2]*y[2])*fac2; //dBzdz

}
