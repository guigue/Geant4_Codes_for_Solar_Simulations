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
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef Solar_FieldMessenger_h
#define Solar_FieldMessenger_h 1

#include "G4UImessenger.hh"

class Solar_FieldSetup;
class G4UIdirectory;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Solar_FieldMessenger: public G4UImessenger
{
  public:
    Solar_FieldMessenger(Solar_FieldSetup* );
    virtual ~Solar_FieldMessenger();

    virtual void SetNewValue(G4UIcommand*, G4String);

  private:

    Solar_FieldSetup*          fEMfieldSetup;

    G4UIdirectory*             fFieldDir;

    G4UIcmdWithAnInteger*      fStepperCmd;

    G4UIcmdWithADoubleAndUnit* fMagFieldCmd;
    G4UIcmdWithADoubleAndUnit* fMinStepCmd;
    G4UIcmdWithADoubleAndUnit* fEpsMinCmd;
    G4UIcmdWithADoubleAndUnit* fEpsMaxCmd;
    G4UIcmdWithADoubleAndUnit* fDeltaChordCmd;
    G4UIcmdWithADoubleAndUnit* fDeltaOneStepCmd;
    G4UIcmdWithADoubleAndUnit* fDeltaIntersectionCmd;
    G4UIcmdWithADoubleAndUnit* fLarAccStepCmd;

    G4UIcmdWithoutParameter*   fUpdateCmd;
};

#endif
