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

#include "Solar_RunAction.hh"
#include "Solar_Analysis.hh"
#include "Solar_SD.hh"

#include "G4Run.hh"
#include "G4ParticleDefinition.hh"
#include "G4Electron.hh"
#include "G4ProcessManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Transportation.hh"
#include "G4CoupledTransportation.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Solar_RunAction::Solar_RunAction()
 : G4UserRunAction()
{

  theWarningEnergy   = 1.0 * keV;  // Arbitrary
  theImportantEnergy = 10.0 * keV;  // Arbitrary
  theNumberOfTrials  = 1000000000;  // Arbitrary

  // For GC Applications should set theNumberOfTrial up to 100 at maximum

  // Applications should determine these thresholds according to
  //  - physics requirements, and
  //  - the computing cost of continuing integration for looping tracks


  // Create analysis manager
// The choice of analysis technology is done via selectin of a namespace
// in B4Analysis.hh
G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
analysisManager->SetVerboseLevel(1);
G4cout << "Using " << analysisManager->GetType()
       << " analysis manager" << G4endl;


// ntuple id = 0
  analysisManager->CreateNtuple("Spectrums", "AtmToCorona");
  analysisManager->CreateNtupleIColumn("TrackID");            // column id = 0
  analysisManager->CreateNtupleIColumn("ParentID");           // column id = 1
  analysisManager->CreateNtupleIColumn("PDG");                // column id = 2
  analysisManager->CreateNtupleDColumn("Ekin");               // column id = 3
  analysisManager->CreateNtupleDColumn("Xpos");               // column id = 4
  analysisManager->CreateNtupleDColumn("Ypos");               // column id = 5
  analysisManager->CreateNtupleDColumn("Zpos");               // column id = 6
  analysisManager->CreateNtupleDColumn("Momentum x");         // column id = 7
  analysisManager->CreateNtupleDColumn("Momentum y");         // column id = 8
  analysisManager->CreateNtupleDColumn("Momentum z");         // column id = 9
  analysisManager->CreateNtupleDColumn("time");               // column id = 10
  analysisManager->CreateNtupleDColumn("weight");             // column id = 11
  analysisManager->CreateNtupleDColumn("Polarization x");     // column id = 12
  analysisManager->CreateNtupleDColumn("Polarization y");     // column id = 13
  analysisManager->CreateNtupleDColumn("Polarization z");     // column id = 14
  analysisManager->CreateNtupleDColumn("VPosition x");        // column id = 15
  analysisManager->CreateNtupleDColumn("VPosition y");        // column id = 16
  analysisManager->CreateNtupleDColumn("VPosition z");        // column id = 17
  analysisManager->CreateNtupleDColumn("VMomentum x");        // column id = 18
  analysisManager->CreateNtupleDColumn("VMomentum y");        // column id = 19
  analysisManager->CreateNtupleDColumn("VMomentum z");        // column id = 20
  analysisManager->CreateNtupleDColumn("VKE");                // column id = 21
  analysisManager->CreateNtupleSColumn("CPName");             // column id = 22*/
  analysisManager->CreateNtupleSColumn("SD1");
  //analysisManager->CreateNtupleSColumn("SD2");
  analysisManager->FinishNtuple(0);

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Solar_RunAction::~Solar_RunAction()
{
  delete G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Solar_RunAction::BeginOfRunAction(const G4Run* aRun)
{
  // Open an output file
  //
  G4String fileName = "analysis/spectrums.g4";
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->OpenFile(fileName);

  //G4RunManager::GetRunManager()->SetRandomNumberStore(true);

  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;

  G4cout << " Calling Solar_RunAction::ChangeLooperParameters() " << G4endl;
  ChangeLooperParameters(G4Electron::Definition() );
}

void Solar_RunAction::
ChangeLooperParameters(const G4ParticleDefinition* particleDef )
{
  //if( particleDef == nullptr )
  //   particleDef = G4Electron::Definition();
  auto transportPair= findTransportation( particleDef );
  auto transport = transportPair.first;
  auto coupledTransport = transportPair.second;

  if( transport != nullptr )
  {
     // Change the values of the looping particle parameters of Transportation
     if( theWarningEnergy   >= 0.0 )
        transport->SetThresholdWarningEnergy(  theWarningEnergy );
     if( theImportantEnergy >= 0.0 )
        transport->SetThresholdImportantEnergy(  theImportantEnergy );
     if( theNumberOfTrials  > 0 )
        transport->SetThresholdTrials( theNumberOfTrials );
  }
  else if( coupledTransport != nullptr )
  {
     // Change the values for Coupled Transport
     if( theWarningEnergy   >= 0.0 )
        coupledTransport->SetThresholdWarningEnergy(  theWarningEnergy );
     if( theImportantEnergy >= 0.0 )
        coupledTransport->SetThresholdImportantEnergy(  theImportantEnergy );
     if( theNumberOfTrials  > 0 )
        coupledTransport->SetThresholdTrials( theNumberOfTrials );
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Solar_RunAction::EndOfRunAction(const G4Run* /*run*/)
{
  // save histograms
  //
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();


  analysisManager->Write();
  analysisManager->CloseFile();

  if( theVerboseLevel > 1 )
     G4cout << G4endl << G4endl
            << " ###########  Track Statistics for Transportation process(es) "
            << " ########### " << G4endl
            << " ############################################## "
            << " ####################### " << G4endl << G4endl;

  auto transportPair= findTransportation( G4Electron::Definition() );
  auto transport = transportPair.first;
  auto coupledTransport = transportPair.second;
  if( transport)             {        transport->PrintStatistics(G4cout); }
  else if( coupledTransport) { coupledTransport->PrintStatistics(G4cout); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

std::pair<G4Transportation*, G4CoupledTransportation*>
Solar_RunAction::findTransportation( const G4ParticleDefinition* particleDef,
                                  bool reportError )
{
  const auto *partPM=  particleDef->GetProcessManager();

  G4VProcess* partTransport = partPM->GetProcess("Transportation");
  auto transport= dynamic_cast<G4Transportation*>(partTransport);

  partTransport = partPM->GetProcess("CoupledTransportation");
  auto coupledTransport=
     dynamic_cast<G4CoupledTransportation*>(partTransport);

  if( reportError && !transport && !coupledTransport )
  {
     G4cerr << "Unable to find Transportation process for particle type "
            << particleDef->GetParticleName()
            << "  ( PDG code = " << particleDef->GetPDGEncoding() << " ) "
            << G4endl;
  }

  return
     std::make_pair( transport, coupledTransport );
         // <G4Transportation*, G4CoupledTransportation*>


}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
