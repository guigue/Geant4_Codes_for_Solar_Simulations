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


#ifndef Solar_FieldSetup_h
#define Solar_FieldSetup_h 1

#include "G4MagneticField.hh"
#include "G4UniformMagField.hh"
#include "G4PropagatorInField.hh"
#include "G4MagIntegratorStepper.hh"
#include "Solar_GC.hh"


//#include "G4ChordFinder.hh"

class G4FieldManager;
class G4ChordFinder;
class G4Mag_UsualEqRhs;
class G4Mag_EqRhs;
class G4EquationOfMotion;
class G4MagIntegratorStepper;
class Solar_FieldMessenger;
class BLineEquation;
class Solar_FieldSetup
{
public:
  Solar_FieldSetup(/*G4double*/);  //  The value of the field
  //Solar_FieldSetup();               //  A zero field

  virtual ~Solar_FieldSetup();

  void SetStepperType( G4int i )
     { fStepperType = i; CreateStepperAndChordFinder(); }

  void SetStepper();

  /// Set the Minimum Step length
  void SetMinStep(G4double s) { fMinStep = s; }

  /// Set the delta chord length
  void SetDeltaChord(G4double dcr) { fDeltaChord = dcr; }

  /// Set the delta one step length
  void SetDeltaOneStep(G4double stp) { fDeltaOneStep = stp; }

  /// Set the delta intersection length
  void SetDeltaIntersection(G4double its) { fDeltaIntersection = its; }

  /// Set the minimum eps length
  void SetEpsMin(G4double eps) { fEpsMin = eps; }

  /// Set the maximum eps length
  void SetEpsMax(G4double eps) { fEpsMax = eps; }

  /// Set Large Acceptable Step
  void SetLargeAccStep(G4double laraccstep) { fLarAccStep = laraccstep; }

  void InitialiseAll();    //  Set parameters and call method below
  void CreateStepperAndChordFinder();

  void SetFieldValue(G4double FootValue);
  //void SetFieldValue(G4double      fieldValue);
  //G4ThreeVector GetFieldValue();

  //void SetDipoleB0(G4double aVal);

protected:

  //magnetic dipole in coordinate
  // Find the global Field Manager

  G4FieldManager*         GetGlobalFieldManager();

  G4FieldManager*         fFieldManager;
  G4ChordFinder*          fChordFinder;
  Solar_GC*               fEquation;
  //G4Mag_UsualEqRhs*       fEquation;
  //G4Mag_EqRhs*            fEquation;
  //G4MagneticField*        fMagneticField;
  Solar_TabulatedField*   fMagneticField;
  G4PropagatorInField*    fFieldPropagator;


  G4MagIntegratorStepper* fStepper;
  G4int                   fStepperType;

  G4double                fMinStep;
  G4double                fDeltaChord;
  G4double                fDeltaOneStep;
  G4double                fDeltaIntersection;
  G4double                fEpsMin;
  G4double                fEpsMax;
  G4double                fLarAccStep;
  G4double                AR;

  Solar_FieldMessenger*      fFieldMessenger;

};

#endif
