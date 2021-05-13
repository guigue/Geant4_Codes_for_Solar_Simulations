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

#include "Solar_FieldSetup.hh"
#include "Solar_FieldMessenger.hh"

#include "Solar_DetectorConstruction.hh"
#include "Solar_Dipole.hh"

#include "G4MagneticField.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4Mag_EqRhs.hh"
#include "Solar_GC.hh"
//#include "BLineEquation.hh"
#include "G4EquationOfMotion.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4ChordFinder.hh"

#include "G4ExplicitEuler.hh"
#include "G4ImplicitEuler.hh"
#include "G4SimpleRunge.hh"
#include "G4SimpleHeum.hh"
#include "G4ClassicalRK4.hh"
#include "G4HelixExplicitEuler.hh"
#include "G4HelixImplicitEuler.hh"
#include "G4HelixSimpleRunge.hh"
#include "G4CashKarpRKF45.hh"
#include "G4RKG3_Stepper.hh"
#include "G4ConstRK4.hh"
#include "G4NystromRK4.hh"
#include "G4HelixMixedStepper.hh"
#include "G4ExactHelixStepper.hh"

// Newest steppers - from Release 10.3-beta (June 2013)
#include "G4BogackiShampine23.hh"
#include "G4BogackiShampine45.hh"
#include "G4DormandPrince745.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//  Constructors:

Solar_FieldSetup::Solar_FieldSetup(/*G4double pGradient*/)
 : fFieldManager(0),
   fChordFinder(0),
   fEquation(0),
   //fEquationType(0),
   fMagneticField(new Solar_Dipole(/*pGradient*/)),
   fStepper(0),
   fStepperType(0),
   fMinStep(0.),
   fFieldMessenger(0),
   fDeltaChord(0),
   fDeltaOneStep(0),
   fDeltaIntersection(0),
   fEpsMin(0),
   fEpsMax(0),
   fLarAccStep(0),
   fFieldPropagator(0)

   {
      //G4cout << " Solar_FieldSetup: magnetic field set to Uniform ( "
      //        << pGradient << " ) " << G4endl;
      InitialiseAll();
   }

/*Solar_FieldSetup::Solar_FieldSetup()
 : fFieldManager(0),
   fChordFinder(0),
   fEquation(0),
   fMagneticField(new Solar_Dipole(G4double())),
   fStepper(0),
   fStepperType(0),
   fMinStep(0.),
   fFieldMessenger(0),
   fDeltaChord(0),
   fDeltaOneStep(0),
   fDeltaIntersection(0),
   fEpsMin(0),
   fEpsMax(0),
   fLarAccStep(0),
   fFieldPropagator(0)

   {
      G4cout << " Solar_FieldSetup: magnetic field set to Uniform (0.0, 0.0) "
             << G4endl;
      InitialiseAll();
   }*/

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Solar_FieldSetup::InitialiseAll()
{
  fFieldMessenger = new Solar_FieldMessenger(this);

  //fEquation = new BLineEquation(fMagneticField);
  //fEquation = new G4Mag_UsualEqRhs(fMagneticField);
  fEquation = new Solar_GC(fMagneticField);

  fMinStep     = 1.0*mm; // minimal step of 1 mm is default

  fStepperType = 4;      // ClassicalRK4 is default stepper

  fDeltaChord = 1.0*mm;

  fDeltaOneStep = 1.0*mm;

  fDeltaIntersection = 1.0*mm;

  fEpsMin = 1.0e-9*mm;

  fEpsMax = 1.0e-8*mm;

  //fLarAccStep = 10* km;


  fFieldManager = G4TransportationManager::GetTransportationManager()
                    ->GetFieldManager();
  CreateStepperAndChordFinder();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Solar_FieldSetup::~Solar_FieldSetup()
{
  delete fMagneticField;
  delete fChordFinder;
  delete fStepper;
  //delete fEquation;
  delete fFieldMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Solar_FieldSetup::CreateStepperAndChordFinder()
{

  //  Get transportation, field, and propagator managers
  G4TransportationManager* transportManager =
         G4TransportationManager::GetTransportationManager();

  // Update field

  G4cout<< " Solar_FieldSetup::CreateStepperAndChordFinder() called "
        << " to reset Stepper."  << G4endl;

  SetStepper();
  //SetEquation();
  //G4cout<<"The minimal step is equal to "<<fMinStep/mm<<" mm"<<G4endl;

  fFieldManager->SetDetectorField(fMagneticField );

  fFieldPropagator = transportManager->GetPropagatorInField();


  if (fChordFinder) delete fChordFinder;

  fChordFinder = new G4ChordFinder( fMagneticField, fMinStep,fStepper );

  //Set Accuracy parameters

 fChordFinder->SetDeltaChord( fDeltaChord);

 fFieldManager->SetAccuraciesWithDeltaOneStep(fDeltaOneStep);
 fFieldManager->SetDeltaIntersection(fDeltaIntersection);
 fFieldManager->SetChordFinder( fChordFinder );

 fFieldPropagator->SetMinimumEpsilonStep(fEpsMin);
 fFieldPropagator->SetMaximumEpsilonStep(fEpsMax);
 //fFieldPropagator->SetLargestAcceptableStep(fLarAccStep);

 G4cout         <<  "##################################################" << G4endl;
 G4cout         <<  "###"<<" Accuracy Parameters: " << G4endl;
 G4cout         <<  "###"<<"           MinStep= " << G4BestUnit(fMinStep,"Length")  << G4endl;
 G4cout         <<  "###"<<"           DeltaChord= " << G4BestUnit(fDeltaChord,"Length") << G4endl;
 G4cout         <<  "###"<<"           DeltaOneStep= " << G4BestUnit(fDeltaOneStep,"Length") << G4endl;
 G4cout         <<  "###"<<"                    " << G4endl;
 G4cout         <<  "###"<<"           DeltaIntersection= " << G4BestUnit(fDeltaIntersection,"Length")  << G4endl;
 G4cout         <<  "###"<<"           EpsMin= " << fEpsMin <<  G4endl;
 G4cout         <<  "###"<<"           EpsMax= " << fEpsMax <<  G4endl;
 G4cout         <<  "###"<<"           LargeAcceptableStep= " << G4BestUnit( fLarAccStep,"Length") << G4endl;
 G4cout         <<  "####################################################" << G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Solar_FieldSetup::SetStepper()
{
// Set stepper according to the stepper type

  if (fStepper) delete fStepper;

  switch ( fStepperType )
  {
    case 0:
      fStepper = new G4ExplicitEuler( fEquation );
      G4cout<<"G4ExplicitEuler is chosen."<<G4endl;
      break;
    case 1:
      fStepper = new G4ImplicitEuler( fEquation );
      G4cout<<"G4ImplicitEuler is chosen"<<G4endl;
      break;
    case 2:
      fStepper = new G4SimpleRunge( fEquation );
      G4cout<<"G4SimpleRunge is chosen"<<G4endl;
      break;
    case 3:
      fStepper = new G4SimpleHeum( fEquation );
      G4cout<<"G4SimpleHeum is chosen"<<G4endl;
      break;
    case 4:
      fStepper = new G4ClassicalRK4( fEquation );
      G4cout<<"G4ClassicalRK4 (default) is chosen"<<G4endl;
      break;
    case 5:
      fStepper = new G4HelixExplicitEuler( fEquation );
      G4cout<<"G4HelixExplicitEuler is chosen"<<G4endl;
      break;
    case 6:
      fStepper = new G4HelixImplicitEuler( fEquation );
      G4cout<<"G4HelixImplicitEuler is chosen"<<G4endl;
      break;
    case 7:
      fStepper = new G4HelixSimpleRunge( fEquation );
      G4cout<<"G4HelixSimpleRunge is chosen"<<G4endl;
      break;
    case 8:
      fStepper = new G4CashKarpRKF45( fEquation );
      G4cout<<"G4CashKarpRKF45 is chosen"<<G4endl;
      break;
    case 9:
      fStepper = new G4RKG3_Stepper( fEquation );
      G4cout<<"G4RKG3_Stepper is chosen"<<G4endl;
      break;
    case 10:
       fStepper = new G4ExactHelixStepper( fEquation );
       G4cout<<"G4ExactHelixStepper is chosen"<<G4endl;
       break;
    case 11:
       fStepper = new G4HelixMixedStepper( fEquation );
       G4cout<<"G4HelixMixedStepper is chosen"<<G4endl;
       break;
    case 12:
       fStepper = new G4ConstRK4( fEquation );
       G4cout<<"G4ConstRK4 Stepper is chosen"<<G4endl;
       break;
    case 13:
      fStepper = new G4NystromRK4( fEquation );
      G4cout<<" G4NystromRK4 Stepper is chosen"<<G4endl;
      break;
    case 14:
    case 23:
      fStepper = new G4BogackiShampine23( fEquation );
      G4cout<<"G4BogackiShampine23 Stepper is chosen"<<G4endl;
      break;
    case 15:
    case 45:
      fStepper = new G4BogackiShampine45( fEquation );
      G4cout<<"G4BogackiShampine45 Stepper is chosen"<<G4endl;
      break;
    case 457:
    case 745:
      fStepper = new G4DormandPrince745( fEquation );
      G4cout<<"G4DormandPrince745 Stepper is chosen"<<G4endl;
      break;
    default:
      fStepper = new G4ClassicalRK4( fEquation );
      G4cout<<"G4ClassicalRK4 Stepper (default) is chosen"<<G4endl;
      // fStepper = new G4DormandPrince745( fEquation );
      // G4cout<<"G4DormandPrince745 (default) Stepper is chosen"<<G4endl;
      break;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
