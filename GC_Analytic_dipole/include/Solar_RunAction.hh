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

#ifndef Solar_RunAction_h
#define Solar_RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
class G4ParticleDefinition;
class G4Transportation;
class G4CoupledTransportation;

class G4Run;

class Solar_RunAction : public G4UserRunAction
{
  public:
    Solar_RunAction();
    virtual ~Solar_RunAction();

  public:

    virtual void BeginOfRunAction(const G4Run*);
    virtual void   EndOfRunAction(const G4Run*);

    // Helper method to change the Transportation's 'looper' parameters
  void ChangeLooperParameters(const G4ParticleDefinition* particleDef );

  // Helper method to find the Transportation process for a particle type
  std::pair<G4Transportation*, G4CoupledTransportation*>
     findTransportation( const G4ParticleDefinition * particleDef,
                         bool reportError= true );

public:
  void     SetNumberOfTrials( G4int val )   { theNumberOfTrials  =  val; }
  void     SetWarningEnergy( double val )   { theWarningEnergy = val; }
  void     SetImportantEnergy( double val ) { theImportantEnergy = val; }
  G4int    GetNumberOfTrials() { return theNumberOfTrials; }
  G4double GetWarningEnergy()  { return theWarningEnergy; }
  G4double GetImportantEnergy() { return theImportantEnergy; }

private:

  // Values for initialising 'loopers' parameters of Transport process
  G4int    theNumberOfTrials  =  0;    //  0  Default will not overwrite
  G4double theWarningEnergy   =  0;  // -1  Default values - non operational
  G4double theImportantEnergy =  0;  // -1  Default - will not overwrite

  int    theVerboseLevel = 0;
};

#endif
