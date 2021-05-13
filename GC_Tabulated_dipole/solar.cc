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
// $Id: exampleSolar_.cc 86065 2014-11-07 08:51:15Z gcosmo $
//
/// \file exampleSolar_.cc
/// \brief Main program of the Solar_ example

#include "Solar_DetectorConstruction.hh"
#include "Solar_ActionInitialization.hh"
#include "Solar_RunAction.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4PhysListFactory.hh"
#include "G4UImanager.hh"
#include "QGSP_BERT.hh"
#include "FTFP_BERT.hh"
#include "QGSP_BIC_HP.hh"
#include "QGSP_BERT_HP.hh"
#include "G4StepLimiterPhysics.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "G4EmParameters.hh"
#include "G4HadronicProcessStore.hh"

#include "Randomize.hh"
#include "time.h"

// For Printing statistic from Transporation process(es)
//#include "G4Electron.hh"
//#include "G4Transportation.hh"
//#include "G4CoupledTransportation.hh"


namespace {
  void PrintUsage() {
    G4cerr << " Usage: " << G4endl;
    G4cerr << " exampleED [-m macro ] [-p physList ] [-u UIsession]  [-t nThreads]" << G4endl;
    G4cerr << "   note: -t option is available only for multi-threaded mode." << G4endl;
    G4cerr << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
  // Evaluate arguments
  //
  if ( argc > 9 ) {
    PrintUsage();
    return 1;
  }

  G4String macro;
  G4String session;
  G4String physicsListName;
  G4int nofThreads = 0;
  for ( G4int i=1; i<argc; i=i+2 ) {
    if      ( G4String(argv[i]) == "-m" ) macro = argv[i+1];
    else if ( G4String(argv[i]) == "-u" ) session = argv[i+1];
    else if ( G4String(argv[i]) == "-p" ) physicsListName = argv[i+1];
#ifdef G4MULTITHREADED
    else if ( G4String(argv[i]) == "-t" ) {
      nofThreads = G4UIcommand::ConvertToInt(argv[i+1]);
    }
#endif
    else {
      PrintUsage();
      return 1;
    }
  }

  // Detect interactive mode (if no arguments) and define UI session
  //
  G4UIExecutive* ui = 0;
  if ( ! macro.size() ) {
    ui = new G4UIExecutive(argc, argv, session);
  }

  //choose the Random engine
 CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine());
 //set random seed with system time
 G4long seed = time(NULL);
 CLHEP::HepRandom::setTheSeed(seed);

  // Construct the run manager
  //
#ifdef G4MULTITHREADED
  G4MTRunManager * runManager = new G4MTRunManager;
  if ( nofThreads > 0 ) {
    runManager->SetNumberOfThreads(nofThreads);
  }
#else
  G4RunManager * runManager = new G4RunManager;
#endif

  // Set mandatory initialization classes
  //
  // Detector construction
  auto detConstruction = new Solar_DetectorConstruction();
  runManager->SetUserInitialization(detConstruction);

  // Configure the use of low thresholds for looping particles
  //  ( appropriate for typical applications using low-energy physics. )
  //auto plHelper = G4PhysicsListHelper::GetPhysicsListHelper();
  //plHelper->UseLowLooperThresholds();
  // Request a set of pre-selected values of the parameters for looping
  //  particles

  auto physicsList = new QGSP_BERT_HP;
  physicsList->RegisterPhysics(new G4StepLimiterPhysics());
  runManager->SetUserInitialization(physicsList);




  // User action initialization
  auto actionInitialization = new Solar_ActionInitialization();
  runManager->SetUserInitialization(actionInitialization);

  // Fine grained control of thresholds for looping particles
  //auto runAction= new Solar_RunAction();
  //runAction->SetWarningEnergy(   10.0 * keV );
              // Looping particles with E < 10 keV will be killed after 1 step
              //   with warning.
              // Looping particles with E > 10 keV will generate a warning.
  //runAction->SetImportantEnergy( 100.0  * keV );
  //runAction->SetNumberOfTrials( 1000000000);
              // Looping particles with E > 0.1 MeV will survive for up to
              //  30 'tracking' steps, and only be killed if they still loop.
  // Note: this mechanism overwrites the thresholds established by
  //       the call to UseLowLooperThresholds() above.

  // Suppress large verbosity from EM & hadronic processes
  G4EmParameters::Instance()->SetVerbose(0);
  G4HadronicProcessStore::Instance()->SetVerbose(0);

   // Initialize G4 kernel
  //
  runManager->Initialize();


  // Initialize visualization
  //
  G4VisManager* visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();

  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  if ( macro.size() ) {
    // batch mode
    G4String command = "/control/execute ";
    UImanager->ApplyCommand(command+macro);
  }
  else {
    // interactive mode : define UI session
    UImanager->ApplyCommand("/control/execute init_vis.mac");
    ui->SessionStart();
    delete ui;
  }

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted
  // in the main() program !

  delete visManager;
  delete runManager;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
