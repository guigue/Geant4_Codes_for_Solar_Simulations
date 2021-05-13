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


#include "Solar_FieldMessenger.hh"

#include "Solar_FieldSetup.hh"
//#include "Solar_Dipole.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Solar_FieldMessenger::Solar_FieldMessenger(Solar_FieldSetup* fieldSetup)
 : G4UImessenger(),
   fEMfieldSetup(fieldSetup),
   fFieldDir(0),
   fStepperCmd(0),
   fMagFieldCmd(0),
   fMinStepCmd(0),
   fEpsMinCmd(0),
   fEpsMaxCmd(0),
   fDeltaChordCmd(0),
   fDeltaOneStepCmd(0),
   fDeltaIntersectionCmd(0),
   fLarAccStepCmd(0),
   fUpdateCmd(0)
{
  fFieldDir = new G4UIdirectory("/field/");
  fFieldDir->SetGuidance("Solar field tracking control.");

  fStepperCmd = new G4UIcmdWithAnInteger("/field/setStepperType",this);
  fStepperCmd->SetGuidance("Select stepper type for magnetic field");
  fStepperCmd->SetParameterName("choice",true);
  fStepperCmd->SetDefaultValue(4);
  fStepperCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fUpdateCmd = new G4UIcmdWithoutParameter("/field/update",this);
  fUpdateCmd->SetGuidance("Update field values.");
  fUpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  fUpdateCmd->SetGuidance("if you changed chord finder or magnetic field value(s).");
  fUpdateCmd->AvailableForStates(G4State_Idle);

  fMagFieldCmd = new G4UIcmdWithADoubleAndUnit("/field/setFieldZ",this);
  fMagFieldCmd->SetGuidance("Define magnetic field.");
  fMagFieldCmd->SetGuidance("Magnetic field will be in Z direction.");
  fMagFieldCmd->SetParameterName("Bz",false,false);
  fMagFieldCmd->SetDefaultUnit("tesla");
  fMagFieldCmd->AvailableForStates(G4State_Idle);

  fMinStepCmd = new G4UIcmdWithADoubleAndUnit("/field/setMinStep",this);
  fMinStepCmd->SetGuidance("Define minimal step");
  fMinStepCmd->SetGuidance("Magnetic field will be in Z direction.");
  fMinStepCmd->SetParameterName("min step",false,false);
  fMinStepCmd->SetDefaultUnit("mm");
  fMinStepCmd->AvailableForStates(G4State_Idle);

  fEpsMinCmd = new G4UIcmdWithADoubleAndUnit("/field/setEpsilonMin",this);
  fEpsMinCmd->SetGuidance("Define epsilon min");
  fEpsMinCmd->SetParameterName("epsilon min",false,false);
  fEpsMinCmd->SetDefaultUnit("mm");
  fEpsMinCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fEpsMaxCmd = new G4UIcmdWithADoubleAndUnit("/field/setEpsilonMax",this);
  fEpsMaxCmd->SetGuidance("Define epsilon max");
  fEpsMaxCmd->SetParameterName("epsilon max",false,false);
  fEpsMaxCmd->SetDefaultUnit("mm");
  fEpsMaxCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fDeltaChordCmd = new G4UIcmdWithADoubleAndUnit("/field/setDeltaChord",this);
  fDeltaChordCmd->SetGuidance("Define delta chord");
  fDeltaChordCmd->SetParameterName("delta chord",false,false);
  fDeltaChordCmd->SetDefaultUnit("mm");
  fDeltaChordCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fDeltaOneStepCmd = new G4UIcmdWithADoubleAndUnit("/field/setDeltaOneStep",this);
  fDeltaOneStepCmd->SetGuidance("Define delta one step");
  fDeltaOneStepCmd->SetParameterName("delta one step",false,false);
  fDeltaOneStepCmd->SetDefaultUnit("mm");
  fDeltaOneStepCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fDeltaIntersectionCmd = new G4UIcmdWithADoubleAndUnit("/field/setDeltaIntersection",this);
  fDeltaIntersectionCmd->SetGuidance("Define delta Intersection");
  fDeltaIntersectionCmd->SetParameterName("delta intersection",false,false);
  fDeltaIntersectionCmd->SetDefaultUnit("mm");
  fDeltaIntersectionCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fLarAccStepCmd = new G4UIcmdWithADoubleAndUnit("/field/setLargeAcceptableStep",this);
  fLarAccStepCmd->SetGuidance("Define Large Acceptable Step");
  fLarAccStepCmd->SetParameterName("Large Acceptable Step",false,false);
  fLarAccStepCmd->SetDefaultUnit("mm");
  fLarAccStepCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Solar_FieldMessenger::~Solar_FieldMessenger()
{
  delete fStepperCmd;
  delete fMagFieldCmd;
  delete fMinStepCmd;
  delete fEpsMinCmd;
  delete fEpsMaxCmd;
  delete fFieldDir;
  delete fUpdateCmd;
  delete fDeltaChordCmd;
  delete fDeltaOneStepCmd;
  delete fDeltaIntersectionCmd;
  delete fLarAccStepCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Solar_FieldMessenger::SetNewValue( G4UIcommand* command, G4String newValue)
{
  if( command == fStepperCmd )
    fEMfieldSetup->SetStepperType(fStepperCmd->GetNewIntValue(newValue));
  if( command == fUpdateCmd )
    fEMfieldSetup->CreateStepperAndChordFinder();
  if( command == fMinStepCmd )
    fEMfieldSetup->SetMinStep(fMinStepCmd->GetNewDoubleValue(newValue));
  if( command == fEpsMinCmd )
    fEMfieldSetup->SetEpsMin(fEpsMinCmd->GetNewDoubleValue(newValue));
  if( command == fEpsMaxCmd )
    fEMfieldSetup->SetEpsMax(fEpsMaxCmd->GetNewDoubleValue(newValue));
  if( command == fDeltaChordCmd )
    fEMfieldSetup->SetDeltaChord(fDeltaChordCmd->GetNewDoubleValue(newValue));
  if( command == fDeltaOneStepCmd )
    fEMfieldSetup->SetDeltaOneStep(fDeltaOneStepCmd->GetNewDoubleValue(newValue));
  if( command == fDeltaIntersectionCmd)
    fEMfieldSetup->SetDeltaIntersection(fDeltaIntersectionCmd->GetNewDoubleValue(newValue));
  if( command == fLarAccStepCmd)
    fEMfieldSetup->SetLargeAccStep(fLarAccStepCmd->GetNewDoubleValue(newValue));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
