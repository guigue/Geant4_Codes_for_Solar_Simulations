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


#include "Solar_GC.hh"
#include "Solar_Dipole.hh"
#include "Solar_FieldSetup.hh"
#include "G4MagneticField.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "vector"
#include "fstream"
#include "globals.hh"    // For DBL_MAX
#include <cmath>

Solar_GC::Solar_GC( Solar_Dipole* MagField )
  : G4Mag_EqRhs( MagField ), charge(0.), mass(0.)
  {}

Solar_GC::~Solar_GC() {}

void
Solar_GC::
 SetChargeMomentumMass( G4ChargeState particleCharge,G4double MomentumXc,G4double particleMass )

{
   G4Mag_EqRhs::SetChargeMomentumMass( particleCharge, MomentumXc, mass );

   charge = particleCharge.GetCharge();
   mass = particleMass;
}

void
Solar_GC::EvaluateRhsGivenB( const G4double y[],const G4double B[],G4double dydx[] ) const
{

  G4ThreeVector unitx(1,0,0);
  G4ThreeVector unity(0,1,0);
  G4ThreeVector unitz(0,0,1);

  // Magnitude of B

  G4double B_magnitude = std::sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]);

  // Unit Vector of B

  G4ThreeVector unitb ((B[0]/B_magnitude), (B[1]/B_magnitude), (B[2]/B_magnitude));

  // Momentum vector

  G4ThreeVector p(y[3],y[4],y[5]);                                              // p[] -> y[] = p[]*c, i.e. units of y[] is MeV

  // Momentum components

  G4ThreeVector p_par  = ((p.dot(unitb))*unitb);                                // paralel component (vector)
  G4ThreeVector p_perp = (p-p_par);                                             // perpendicular component (vector)

  // Momentum magnitude

  G4double p_magnitude = std::sqrt(y[3]*y[3] + y[4]*y[4] + y[5]*y[5]);

  // Inverse of momentum magnitude

  G4double inv_momentum_magnitude = 1.0/p_magnitude;

  // Magnitudes of momentum components

  G4double p_par_magnitude  = std::sqrt(p_par[0]*p_par[0]+p_par[1]*p_par[1]+p_par[2]*p_par[2]);
  G4double p_perp_magnitude = std::sqrt(p_perp[0]*p_perp[0]+p_perp[1]*p_perp[1]+p_perp[2]*p_perp[2]);

  // Total energy, gamma and cof

  G4double Energy = std::sqrt(sqr(p_magnitude)+sqr(mass));
  G4double xgamma = Energy/(mass);
  G4double cof = FCof()*inv_momentum_magnitude;   // FCof() = 299.72 mm/ns -> 0.29972 MeV/(tesla*mm)

  // Velocity vector

  G4ThreeVector v = p/(mass*xgamma);                                            // v[] -> beta =v[]/c

  // Velocity components from momentum components

  G4ThreeVector v_par  = p_par/(mass*xgamma);
  G4ThreeVector v_perp = p_perp/(mass*xgamma);

  // Velocity magnitude

  G4double v_magnitude  = std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);

  // Magnitudes of velocity components

  //G4double v_par_magnitude  = std::sqrt(v_par[0]*v_par[0]+v_par[1]*v_par[1]+v_par[2]*v_par[2]);
  //G4double v_perp_magnitude = std::sqrt(v_perp[0]*v_perp[0]+v_perp[1]*v_perp[1]+v_perp[2]*v_perp[2]);

  // Unit vectors xi and zeta

  if ( unitb[2] ==  1.0 ){unitb[2] =  0.999999;}
  if ( unitb[2] == -1.0 ){unitb[2] = -0.999999;}

  G4ThreeVector xi   = (1.0/std::sqrt(1.0-unitb[2]*unitb[2]))*(unitb.cross(unitz)); // unit vector
  G4ThreeVector zeta = xi.cross(unitb);                                         // unit vector

  // Delta and Epsilon

  G4double delta   = p_perp.dot(xi)/p_perp_magnitude;                           // scalar
  G4double epsilon = p_perp.dot(zeta)/p_perp_magnitude;                         // scalar

  // Operations with delta and xi

  G4ThreeVector deltaXi = delta*xi;                                             // vector
  G4ThreeVector epsilonZeta = epsilon*zeta;                                     // vector
  G4ThreeVector deltaXiPLUSepsilonZeta = deltaXi + epsilonZeta;                 // vector

  // Gradient of B

  //--------------------------------------------------------------------
  //Check if it's correct B[0]+d or should be expressed as B[0+d] in C++
  //---------------------------------------------------------------------
  // The gradB is broken and declared as G4double gradBx, gradBy and gradBz instead
  // of G4ThreeVector to adapt it to what G4 wanted: dy/dx[3], dy/dx[4], and dy/dx[5].

  G4ThreeVector gradientB((unitb[0]*B[3] + unitb[1]*B[6] + unitb[2]*B[9]),
                          (unitb[0]*B[4] + unitb[1]*B[7] + unitb[2]*B[10]),
                          (unitb[0]*B[5] + unitb[1]*B[8] + unitb[2]*B[11]));

  // derivatives of b and Gradient of b

  G4double bx_compx = (B[3]/B_magnitude) - (unitb[0]*gradientB[0]/B_magnitude);
  G4double bx_compy = (B[4]/B_magnitude) - (unitb[0]*gradientB[1]/B_magnitude);
  G4double bx_compz = (B[5]/B_magnitude) - (unitb[0]*gradientB[2]/B_magnitude);

  G4double by_compx = (B[6]/B_magnitude) - (unitb[1]*gradientB[0]/B_magnitude);
  G4double by_compy = (B[7]/B_magnitude) - (unitb[1]*gradientB[1]/B_magnitude);
  G4double by_compz = (B[8]/B_magnitude) - (unitb[1]*gradientB[2]/B_magnitude);

  G4double bz_compx = (B[9]/B_magnitude)  - (unitb[2]*gradientB[0]/B_magnitude);
  G4double bz_compy = (B[10]/B_magnitude) - (unitb[2]*gradientB[1]/B_magnitude);
  G4double bz_compz = (B[11]/B_magnitude) - (unitb[2]*gradientB[2]/B_magnitude);

  G4ThreeVector b_gradientUnitb((unitb[0]*bx_compx + unitb[1]*bx_compy + unitb[2]*bx_compz),
                                (unitb[0]*by_compx + unitb[1]*by_compy + unitb[2]*by_compz),
                                (unitb[0]*bz_compx + unitb[1]*bz_compy + unitb[2]*bz_compz));

  // derivatives of Xi and Gradient of Xi

  G4double xiX_compx = ((1.0/(std::sqrt(1.0-unitb[2]*unitb[2]))*by_compx) + (unitb[1]*unitb[2]/(std::pow(1.0-unitb[2]*unitb[2],3.0/2.0)))*bz_compx);
  G4double xiX_compy = ((1.0/(std::sqrt(1.0-unitb[2]*unitb[2]))*by_compy) + (unitb[1]*unitb[2]/(std::pow(1.0-unitb[2]*unitb[2],3.0/2.0)))*bz_compy);
  G4double xiX_compz = ((1.0/(std::sqrt(1.0-unitb[2]*unitb[2]))*by_compz) + (unitb[1]*unitb[2]/(std::pow(1.0-unitb[2]*unitb[2],3.0/2.0)))*bz_compz);

  G4double xiY_compx = ((-1.0/(std::sqrt(1.0-unitb[2]*unitb[2]))*bx_compx) - (unitb[0]*unitb[2]/(std::pow(1.0-unitb[2]*unitb[2],3.0/2.0)))*bz_compx);
  G4double xiY_compy = ((-1.0/(std::sqrt(1.0-unitb[2]*unitb[2]))*bx_compy) - (unitb[0]*unitb[2]/(std::pow(1.0-unitb[2]*unitb[2],3.0/2.0)))*bz_compy);
  G4double xiY_compz = ((-1.0/(std::sqrt(1.0-unitb[2]*unitb[2]))*bx_compz) - (unitb[0]*unitb[2]/(std::pow(1.0-unitb[2]*unitb[2],3.0/2.0)))*bz_compz);


  G4ThreeVector b_gradientXi((xiX_compx*unitb[0] + xiX_compy*unitb[1] + xiX_compz*unitb[2]),(xiY_compx*unitb[0] + xiY_compy*unitb[1] + xiY_compz*unitb[2]),0.0);

  // derivatives of Zeta and Gradient of Zeta

  G4double zetaX_compx = (-1.0/(std::sqrt(1.0-unitb[2]*unitb[2]))*(unitb[0]*bz_compx+unitb[2]*bx_compx)-((unitb[0]*unitb[2]*unitb[2])/(std::pow((1.0-unitb[2]*unitb[2]),(3.0/2.0)))*bz_compx));
  G4double zetaX_compy = (-1.0/(std::sqrt(1.0-unitb[2]*unitb[2]))*(unitb[0]*bz_compy+unitb[2]*bx_compy)-((unitb[0]*unitb[2]*unitb[2])/(std::pow((1.0-unitb[2]*unitb[2]),(3.0/2.0)))*bz_compy));
  G4double zetaX_compz = (-1.0/(std::sqrt(1.0-unitb[2]*unitb[2]))*(unitb[0]*bz_compz+unitb[2]*bx_compz)-((unitb[0]*unitb[2]*unitb[2])/(std::pow((1.0-unitb[2]*unitb[2]),(3.0/2.0)))*bz_compz));

  G4double zetaY_compx = (-1.0/(std::sqrt(1.0-unitb[2]*unitb[2]))*(unitb[1]*bz_compx+unitb[2]*by_compx)-((unitb[1]*unitb[2]*unitb[2])/(std::pow((1.0-unitb[2]*unitb[2]),(3.0/2.0)))*bz_compx));
  G4double zetaY_compy = (-1.0/(std::sqrt(1.0-unitb[2]*unitb[2]))*(unitb[1]*bz_compy+unitb[2]*by_compy)-((unitb[1]*unitb[2]*unitb[2])/(std::pow((1.0-unitb[2]*unitb[2]),(3.0/2.0)))*bz_compy));
  G4double zetaY_compz = (-1.0/(std::sqrt(1.0-unitb[2]*unitb[2]))*(unitb[1]*bz_compz+unitb[2]*by_compz)-((unitb[1]*unitb[2]*unitb[2])/(std::pow((1.0-unitb[2]*unitb[2]),(3.0/2.0)))*bz_compz));

  G4double zetaZ_compx = (-unitb[2]/(std::sqrt(1.0-unitb[2]*unitb[2]))*bz_compx);
  G4double zetaZ_compy = (-unitb[2]/(std::sqrt(1.0-unitb[2]*unitb[2]))*bz_compy);
  G4double zetaZ_compz = (-unitb[2]/(std::sqrt(1.0-unitb[2]*unitb[2]))*bz_compz);

  G4ThreeVector b_gradientZeta((zetaX_compx*unitb[0] + zetaX_compy*unitb[1] + zetaX_compz*unitb[2]),
                               (zetaY_compx*unitb[0] + zetaY_compy*unitb[1] + zetaY_compz*unitb[2]),
                               (zetaZ_compx*unitb[0] + zetaZ_compy*unitb[1] + zetaZ_compz*unitb[2]));

  // ***********************************************************************

  // GC equations of motion neglecting drift

  // ***********************************************************************

/*  G4double parcf1 = p_perp_magnitude*p_perp_magnitude*unitb.dot(gradientB)/(2.0*p.dot(unitb)*B_magnitude);
  G4double parcf2 = p.dot(unitb);

  G4double coef1  = p_perp_magnitude*unitb.dot(gradientB)/(2.0*B_magnitude);
  G4double coef2  = p_perp_magnitude;

  G4double dparXdt = (-parcf1*unitb[0]) + (parcf2*b_gradientUnitb[0]);
  G4double dparYdt = (-parcf1*unitb[1]) + (parcf2*b_gradientUnitb[1]);
  G4double dparZdt = (-parcf1*unitb[2]) + (parcf2*b_gradientUnitb[2]);

  dydx[0] = unitb[0]*v.dot(unitb)/v_magnitude;
  dydx[1] = unitb[1]*v.dot(unitb)/v_magnitude;
  dydx[2] = unitb[2]*v.dot(unitb)/v_magnitude;

  // ***********************************************************************

  dydx[3] = charge*(coef1*(deltaXiPLUSepsilonZeta[0]) + coef2*(delta*b_gradientXi[0]+epsilon*b_gradientZeta[0]) + (dparXdt))*v.dot(unitb)/v_magnitude;
  dydx[4] = charge*(coef1*(deltaXiPLUSepsilonZeta[1]) + coef2*(delta*b_gradientXi[1]+epsilon*b_gradientZeta[1]) + (dparYdt))*v.dot(unitb)/v_magnitude;
  dydx[5] = charge*(coef1*(deltaXiPLUSepsilonZeta[2]) + coef2*(delta*b_gradientXi[2]+epsilon*b_gradientZeta[2]) + (dparZdt))*v.dot(unitb)/v_magnitude;
*/
  // ***********************************************************************

  // GC equations of motion including drift

  // ***********************************************************************

  G4ThreeVector bcrossgradB = unitb.cross(gradientB);

  G4double drcoef = (p_magnitude*p_magnitude + p_par_magnitude*p_par_magnitude)/(2.0*p.dot(unitb)*B_magnitude*B_magnitude*FCof());

  G4double parcf1 = p_perp_magnitude*p_perp_magnitude*unitb.dot(gradientB)/(2.0*p.dot(unitb)*B_magnitude);
  G4double parcf2 = p.dot(unitb);

  G4double coef1  = p_perp_magnitude*unitb.dot(gradientB)/(2.0*B_magnitude);
  G4double coef2  = p_perp_magnitude;

  G4double dparXdt = (-parcf1*unitb[0]) + (parcf2*b_gradientUnitb[0]);
  G4double dparYdt = (-parcf1*unitb[1]) + (parcf2*b_gradientUnitb[1]);
  G4double dparZdt = (-parcf1*unitb[2]) + (parcf2*b_gradientUnitb[2]);

  dydx[0] = unitb[0]*v.dot(unitb)/v_magnitude + drcoef*bcrossgradB[0]*v.dot(unitb)/v_magnitude;
  dydx[1] = unitb[1]*v.dot(unitb)/v_magnitude + drcoef*bcrossgradB[1]*v.dot(unitb)/v_magnitude;
  dydx[2] = unitb[2]*v.dot(unitb)/v_magnitude + drcoef*bcrossgradB[2]*v.dot(unitb)/v_magnitude;

// ***********************************************************************

  dydx[3] = (coef1*(deltaXiPLUSepsilonZeta[0]) + coef2*(delta*b_gradientXi[0]+epsilon*b_gradientZeta[0]) + (dparXdt))*v.dot(unitb)/v_magnitude;
  dydx[4] = (coef1*(deltaXiPLUSepsilonZeta[1]) + coef2*(delta*b_gradientXi[1]+epsilon*b_gradientZeta[1]) + (dparYdt))*v.dot(unitb)/v_magnitude;
  dydx[5] = (coef1*(deltaXiPLUSepsilonZeta[2]) + coef2*(delta*b_gradientXi[2]+epsilon*b_gradientZeta[2]) + (dparZdt))*v.dot(unitb)/v_magnitude;

//   ***********************************************************************


  // NL equations of motion

//   dydx[0] = y[3]*inv_momentum_magnitude;       // (d/ds)x = Vx/V
//   dydx[1] = y[4]*inv_momentum_magnitude;       // (d/ds)y = Vy/V
//   dydx[2] = y[5]*inv_momentum_magnitude;       // (d/ds)z = Vz/V

//  ***********************************************************************

//   dydx[3] = cof*(y[4]*B[2] - y[5]*B[1]);       // Ax = a*(Vy*Bz - Vz*By)
//   dydx[4] = cof*(y[5]*B[0] - y[3]*B[2]);       // Ay = a*(Vz*Bx - Vx*Bz)
//   dydx[5] = cof*(y[3]*B[1] - y[4]*B[0]);       // Az = a*(Vx*By - Vy*Bx)

  // ************************************************************************

  return ;
}
