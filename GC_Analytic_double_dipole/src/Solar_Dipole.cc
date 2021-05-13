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
  G4double Mf = B0F*pow(Rf*Rf+dzF*dzF,5.0/2.0)/sqrt(4.0*Rf*Rf*Rf*Rf + 5.0*(Rf*Rf)*(dzF*dzF) + dzF*dzF*dzF*dzF);
  G4double Ms = B0S*pow(Rf2*Rf2+dzS*dzS,5.0/2.0)/sqrt(4.0*Rf2*Rf2*Rf2*Rf2 + 5.0*(Rf2*Rf2)*(dzS*dzS) + dzS*dzS*dzS*dzS);

  //G4double fac =  Mf/pow(x*x+y*y+z*z, 5.0/2.0);

  // factor for the first dipole and its derivate
  G4double facF = Mf/pow((y[0]+dxF)*(y[0]+dxF)+(y[1]+dyF)*(y[1]+dyF)+(y[2]+dzF)*(y[2]+dzF), 5.0/2.0);
  G4double facF2 = Mf/pow((dxF+y[0])*(dxF+y[0])+(dyF+y[1])*(dyF+y[1])+(dzF+y[2])*(dzF+y[2]), 7.0/2.0);

  // factor for the second dipole
  G4double facS = -Ms/pow((y[0]+dxS)*(y[0]+dxS)+(y[1]+dyS)*(y[1]+dyS)+(y[2]+dzS)*(y[2]+dzS), 5.0/2.0);
  G4double facS2 = -Ms/pow((dxS+y[0])*(dxS+y[0])+(dyS+y[1])*(dyS+y[1])+(dzS+y[2])*(dzS+y[2]), 7.0/2.0);

  // Bx, By and Bz values for the First and Second dipoles
  G4double Bx1 = (2*(y[0]+dxF)*(y[0]+dxF)-(y[1]+dyF)*(y[1]+dyF)-(y[2]+dzF)*(y[2]+dzF))*facF;
  G4double By1 = 3*(y[0]+dxF)*(y[1]+dyF)*facF;
  G4double Bz1 = 3*(y[0]+dxF)*(y[2]+dzF)*facF;

  G4double Bx2 = (2*(y[0]+dxS)*(y[0]+dxS)-(y[1]+dyS)*(y[1]+dyS)-(y[2]+dzS)*(y[2]+dzS))*facS;
  G4double By2 = 3*(y[0]+dxS)*(y[1]+dyS)*facS;
  G4double Bz2 = 3*(y[0]+dxS)*(y[2]+dzS)*facS;

  // Derivatives for the First and Second dipoles

G4double dBxdx1 = -((3.0*(dxF+y[0])*2.0*dxF*dxF-3.0*dyF*dyF-3.0*dzF*dzF+4.0*dxF*y[0]+2.0*y[0]*y[0]-6.0*dyF*y[1]-3.0*y[1]*y[1]-6.0*dzF*y[2]-3.0*y[2]*y[2])*facF2); //dBxdxF
G4double dBxdy1 =  3.0*(dyF+y[1])*(-4.0*dxF*dxF+dyF*dyF+dzF*dzF-8.0*dxF*y[0]-4*y[0]*y[0]+2.0*dyF*y[1]+y[1]*y[1]+2.0*dzF*y[2]+y[2]*y[2])*facF2; //dBxdyF
G4double dBxdz1 =  3.0*(dzF+y[2])*(-4.0*dxF*dxF+dyF*dyF+dzF*dzF-8.0*dxF*y[0]- 4.0*y[0]*y[0] + 2.0*dyF*y[1]+ y[1]*y[1] + 2.0*dzF*y[2] + y[2]*y[2])*facF2; //dBxdzF
G4double dBydx1 =  3.0*(dyF+y[1])*(-4.0*(dxF+y[0])*(dxF+y[0])+(dyF+y[1])*(dyF+y[1])+(dzF+y[2])*(dzF+y[2]))*facF2; //dBydxF
G4double dBydy1 =  3.0*(dxF+y[0])*((dxF+y[0])*(dxF+y[0])-4.0*(dyF+y[1])*(dyF+y[1])+(dzF+y[2])*(dzF+y[2]))*facF2; //dBydyF
G4double dBydz1 = -(15.0*(dxF+y[0])*(dyF+y[1])*(dzF + y[2])*facF2); //dBydzF
G4double dBzdx1 =  3.0*(y[2]+dzF)*(-4.0*(dxF+y[0])*(dxF+y[0])+(dyF+y[1])*(dyF+y[1])+(dzF+y[2])*(dzF+y[2]))*facF2; //dBzdxF
G4double dBzdy1  = -15.0*(dxF+y[0])*(dyF+y[1])*(dzF + y[2])*facF2; //dBzdyF
G4double dBzdz1  = 3.0*(dxF+y[0])*((dxF+y[0])*(dxF+y[0])+(dyF+y[1])*(dyF+y[1])-4.0*(dzF+y[2])*(dzF+y[2]))*facF2; //dBzdzF


G4double dBxdx2 = -((3.0*(dxS+y[0])*2.0*dxS*dxS-3.0*dyS*dyS-3.0*dzS*dzS+4.0*dxS*y[0]+2.0*y[0]*y[0]-6.0*dyS*y[1]-3.0*y[1]*y[1]-6.0*dzS*y[2]-3.0*y[2]*y[2])*facS2); //dBxdxS
G4double dBxdy2 =  3.0*(dyS+y[1])*(-4.0*dxS*dxS+dyS*dyS+dzS*dzS-8.0*dxS*y[0]-4*y[0]*y[0]+2.0*dyS*y[1]+y[1]*y[1]+2.0*dzS*y[2]+y[2]*y[2])*facS2; //dBxdyS
G4double dBxdz2 =  3.0*(dzS+y[2])*(-4.0*dxS*dxS+dyS*dyS+dzS*dzS-8.0*dxS*y[0]- 4.0*y[0]*y[0] + 2.0*dyS*y[1]+ y[1]*y[1] + 2.0*dzS*y[2] + y[2]*y[2])*facS2; //dBxdzS
G4double dBydx2 =  3.0*(dyS+y[1])*(-4.0*(dxS+y[0])*(dxS+y[0])+(dyS+y[1])*(dyS+y[1])+(dzS+y[2])*(dzS+y[2]))*facS2; //dBydxS
G4double dBydy2 =  3.0*(dxS+y[0])*((dxS+y[0])*(dxS+y[0])-4.0*(dyS+y[1])*(dyS+y[1])+(dzS+y[2])*(dzS+y[2]))*facS2; //dBydyS
G4double dBydz2 = -(15.0*(dxS+y[0])*(dyS+y[1])*(dzS + y[2])*facS2); //dBydzS
G4double dBzdx2 =  3.0*(y[2]+dzS)*(-4.0*(dxS+y[0])*(dxS+y[0])+(dyS+y[1])*(dyS+y[1])+(dzS+y[2])*(dzS+y[2]))*facS2; //dBzdxS
G4double dBzdy2  = -15.0*(dxS+y[0])*(dyS+y[1])*(dzS + y[2])*facS2; //dBzdyS
G4double dBzdz2  = 3.0*(dxS+y[0])*((dxS+y[0])*(dxS+y[0])+(dyS+y[1])*(dyS+y[1])-4.0*(dzS+y[2])*(dzS+y[2]))*facS2; //dBzdzS

// Values of double dipole to pass at Solar_GC class
  //Bx, By, Bz;
  B[0] = Bx1+Bx2;
  B[1] = By1+By2;
  B[2] = Bz1+Bz2;
  // Gradients
  B[3] = dBxdx1+dBxdx2; //dBxdx
  B[4] = dBxdy1+dBxdy2; //dBxdy
  B[5] = dBxdz1+dBxdz2;//dBxdz
  B[6] = dBydx1+dBydx2;//dBydx
  B[7] = dBydy1+dBydy2;//dBydy
  B[8] = dBydz1+dBydz2;//dBydz
  B[9] = dBzdx1+dBzdx2; //dBzdx
  B[10] = dBzdy1+dBzdy2; //dBzdy
  B[11] = dBzdz1+dBzdz2;//dBzdz

}
