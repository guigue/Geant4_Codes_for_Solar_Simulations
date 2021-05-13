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


#include "Solar_SD.hh"
#include "Solar_Analysis.hh"
#include "Solar_RunAction.hh"


#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4VTouchable.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4ios.hh"
#include "G4VProcess.hh"
#include "G4SystemOfUnits.hh"
#include "G4TouchableHistory.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Solar_SD::Solar_SD(const G4String& name, G4int ntupleId)
 : G4VSensitiveDetector(name),
   fNtupleId(0)
   //fPProtonsKE(0),
   //fGammasKE(0),
   //fElectronsKE(0),
   //fPositronsKE(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Solar_SD::~Solar_SD()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Solar_SD::Initialize(G4HCofThisEvent* /*hce*/)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


G4bool Solar_SD::ProcessHits(G4Step* step, G4TouchableHistory*)
{

//Current Track
G4Track*aTrack = step->GetTrack();
//Track ID
G4int TrackID = aTrack->GetTrackID();
//Parent ID
G4int ParentID = aTrack->GetParentID();
//code PDG
G4int pdgCode = aTrack-> GetDefinition()->GetPDGEncoding();
// preStepPoint
G4StepPoint* preStepPoint = step->GetPreStepPoint();
//G4StepPoint* postStepPoint = step->GetPostStepPoint();
//if(!(postStepPoint->GetStepStatus() == fGeomBoundary)) return true;

const G4String particle = step->GetTrack()->GetDefinition()->GetParticleName();

const G4VProcess* CurrentProcess=preStepPoint->GetProcessDefinedStep();

// Obtain local coordinates:
  // const G4VTouchable* touchable = preStepPoint->GetTouchable();
  // G4ThreeVector globalPosition = preStepPoint->GetPosition();
  // G4ThreeVector localPosition
  //   = touchable->GetHistory()->GetTopTransform().TransformPoint(globalPosition);

G4ThreeVector position = preStepPoint->GetPosition();
G4double time   = preStepPoint->GetGlobalTime();
//G4double localTime = preStepPoint->GetLocalTime();
  // Weight:
G4double weight = preStepPoint->GetWeight();
  // Kenergy
//G4double KinE = step->GetTrack()->GetKineticEnergy();
G4double KinE = preStepPoint->GetKineticEnergy();

G4ThreeVector pol = preStepPoint->GetPolarization();

G4ThreeVector Momentum = preStepPoint->GetMomentumDirection();

  // Vertex informations
    //Position
const  G4ThreeVector VPos = step->GetTrack()->GetVertexPosition();
G4double VPosx = VPos.x();
G4double VPosy = VPos.y();
G4double VPosz = VPos.z();
  //Momentum
const G4ThreeVector VMomentum = step->GetTrack()->GetVertexMomentumDirection();
    //Energy
const G4double VKinE = step->GetTrack()->GetVertexKineticEnergy();
    //Process
G4String CPName;
 if(step->GetTrack()->GetCreatorProcess()!=0)
    CPName = step->GetTrack()->GetCreatorProcess()->GetProcessName();

G4String SDet1 = preStepPoint->GetSensitiveDetector()->GetName();
//G4String SDet2 = postStepPoint->GetSensitiveDetector()->GetName();


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%TESTING%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//

G4bool condition = step->IsFirstStepInVolume();
//G4bool condition2 = step->IsLastStepInVolume();


if (CurrentProcess != 0)
{
  const G4String StepProcessName = CurrentProcess->GetProcessName();
  


  if (SDet1 == "CoronaSD" && condition == true &&  particle == "gamma")
  {
  // Store hit in the ntuple
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  analysisManager->FillNtupleIColumn(0, 0, TrackID);
  analysisManager->FillNtupleIColumn(0, 1, ParentID);
  analysisManager->FillNtupleIColumn(0, 2, pdgCode);
  analysisManager->FillNtupleDColumn(0, 3, KinE/MeV);
  analysisManager->FillNtupleDColumn(0, 4, position.x()/cm);
  analysisManager->FillNtupleDColumn(0, 5, position.y()/cm);
  analysisManager->FillNtupleDColumn(0, 6, position.z()/cm);
  analysisManager->FillNtupleDColumn(0, 7, Momentum.x());
  analysisManager->FillNtupleDColumn(0, 8, Momentum.y());
  analysisManager->FillNtupleDColumn(0, 9, Momentum.z());
  analysisManager->FillNtupleDColumn(0, 10, time/s);
  analysisManager->FillNtupleDColumn(0, 11, weight);
  analysisManager->FillNtupleDColumn(0, 12, pol.x());
  analysisManager->FillNtupleDColumn(0, 13, pol.y());
  analysisManager->FillNtupleDColumn(0, 14, pol.z());
  analysisManager->FillNtupleDColumn(0, 15, VPosx/cm);
  analysisManager->FillNtupleDColumn(0, 16, VPosy/cm);
  analysisManager->FillNtupleDColumn(0, 17, VPosz/cm);
  analysisManager->FillNtupleDColumn(0, 18, VMomentum.x());
  analysisManager->FillNtupleDColumn(0, 19, VMomentum.y());
  analysisManager->FillNtupleDColumn(0, 20, VMomentum.z());
  analysisManager->FillNtupleDColumn(0, 21, VKinE/MeV);
  analysisManager->FillNtupleSColumn(0, 22, CPName);
  analysisManager->FillNtupleSColumn(0, 23, SDet1);
  //analysisManager->FillNtupleSColumn(0, 24, SDet2);

  analysisManager->AddNtupleRow(0);
  }



//if(StepProcessName=="Transportation" && SD1 == "CoronaSD" && particle == "gamma"){
if (SDet1 == "CoronaSD" && condition == true && particle == "e-")
{
// Store hit in the ntuple
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  analysisManager->FillNtupleIColumn(0, 0, TrackID);
  analysisManager->FillNtupleIColumn(0, 1, ParentID);
  analysisManager->FillNtupleIColumn(0, 2, pdgCode);
  analysisManager->FillNtupleDColumn(0, 3, KinE/MeV);
  analysisManager->FillNtupleDColumn(0, 4, position.x()/cm);
  analysisManager->FillNtupleDColumn(0, 5, position.y()/cm);
  analysisManager->FillNtupleDColumn(0, 6, position.z()/cm);
  analysisManager->FillNtupleDColumn(0, 7, Momentum.x());
  analysisManager->FillNtupleDColumn(0, 8, Momentum.y());
  analysisManager->FillNtupleDColumn(0, 9, Momentum.z());
  analysisManager->FillNtupleDColumn(0, 10, time/s);
  analysisManager->FillNtupleDColumn(0, 11, weight);
  analysisManager->FillNtupleDColumn(0, 12, pol.x());
  analysisManager->FillNtupleDColumn(0, 13, pol.y());
  analysisManager->FillNtupleDColumn(0, 14, pol.z());
  analysisManager->FillNtupleDColumn(0, 15, VPosx/cm);
  analysisManager->FillNtupleDColumn(0, 16, VPosy/cm);
  analysisManager->FillNtupleDColumn(0, 17, VPosz/cm);
  analysisManager->FillNtupleDColumn(0, 18, VMomentum.x());
  analysisManager->FillNtupleDColumn(0, 19, VMomentum.y());
  analysisManager->FillNtupleDColumn(0, 20, VMomentum.z());
  analysisManager->FillNtupleDColumn(0, 21, VKinE/MeV);
  analysisManager->FillNtupleSColumn(0, 22, CPName);
  analysisManager->FillNtupleSColumn(0, 23, SDet1);
  //analysisManager->FillNtupleSColumn(0, 24, SD2);

  analysisManager->AddNtupleRow(0);


}


if(SDet1 == "CoronaSD" && condition == true && particle == "e+"){

G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

analysisManager->FillNtupleIColumn(0, 0, TrackID);
analysisManager->FillNtupleIColumn(0, 1, ParentID);
analysisManager->FillNtupleIColumn(0, 2, pdgCode);
analysisManager->FillNtupleDColumn(0, 3, KinE/MeV);
analysisManager->FillNtupleDColumn(0, 4, position.x()/cm);
analysisManager->FillNtupleDColumn(0, 5, position.y()/cm);
analysisManager->FillNtupleDColumn(0, 6, position.z()/cm);
analysisManager->FillNtupleDColumn(0, 7, Momentum.x());
analysisManager->FillNtupleDColumn(0, 8, Momentum.y());
analysisManager->FillNtupleDColumn(0, 9, Momentum.z());
analysisManager->FillNtupleDColumn(0, 10, time/s);
analysisManager->FillNtupleDColumn(0, 11, weight);
analysisManager->FillNtupleDColumn(0, 12, pol.x());
analysisManager->FillNtupleDColumn(0, 13, pol.y());
analysisManager->FillNtupleDColumn(0, 14, pol.z());
analysisManager->FillNtupleDColumn(0, 15, VPosx/cm);
analysisManager->FillNtupleDColumn(0, 16, VPosy/cm);
analysisManager->FillNtupleDColumn(0, 17, VPosz/cm);
analysisManager->FillNtupleDColumn(0, 18, VMomentum.x());
analysisManager->FillNtupleDColumn(0, 19, VMomentum.y());
analysisManager->FillNtupleDColumn(0, 20, VMomentum.z());
analysisManager->FillNtupleDColumn(0, 21, VKinE/MeV);
analysisManager->FillNtupleSColumn(0, 22, CPName);
analysisManager->FillNtupleSColumn(0, 23, SDet1);
//analysisManager->FillNtupleSColumn(0, 24, SD2);

  analysisManager->AddNtupleRow(0);
}
if(SDet1 == "CoronaSD" && condition == true && particle == "proton"){

G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

analysisManager->FillNtupleIColumn(0, 0, TrackID);
analysisManager->FillNtupleIColumn(0, 1, ParentID);
analysisManager->FillNtupleIColumn(0, 2, pdgCode);
analysisManager->FillNtupleDColumn(0, 3, KinE/MeV);
analysisManager->FillNtupleDColumn(0, 4, position.x()/cm);
analysisManager->FillNtupleDColumn(0, 5, position.y()/cm);
analysisManager->FillNtupleDColumn(0, 6, position.z()/cm);
analysisManager->FillNtupleDColumn(0, 7, Momentum.x());
analysisManager->FillNtupleDColumn(0, 8, Momentum.y());
analysisManager->FillNtupleDColumn(0, 9, Momentum.z());
analysisManager->FillNtupleDColumn(0, 10, time/s);
analysisManager->FillNtupleDColumn(0, 11, weight);
analysisManager->FillNtupleDColumn(0, 12, pol.x());
analysisManager->FillNtupleDColumn(0, 13, pol.y());
analysisManager->FillNtupleDColumn(0, 14, pol.z());
analysisManager->FillNtupleDColumn(0, 15, VPosx/cm);
analysisManager->FillNtupleDColumn(0, 16, VPosy/cm);
analysisManager->FillNtupleDColumn(0, 17, VPosz/cm);
analysisManager->FillNtupleDColumn(0, 18, VMomentum.x());
analysisManager->FillNtupleDColumn(0, 19, VMomentum.y());
analysisManager->FillNtupleDColumn(0, 20, VMomentum.z());
analysisManager->FillNtupleDColumn(0, 21, VKinE/MeV);
analysisManager->FillNtupleSColumn(0, 22, CPName);
analysisManager->FillNtupleSColumn(0, 23, SDet1);
//analysisManager->FillNtupleSColumn(0, 24, SD2);

  analysisManager->AddNtupleRow(0);
}

if(SDet1 == "CoronaSD" && condition == true && particle == "neutron"){

G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

analysisManager->FillNtupleIColumn(0, 0, TrackID);
analysisManager->FillNtupleIColumn(0, 1, ParentID);
analysisManager->FillNtupleIColumn(0, 2, pdgCode);
analysisManager->FillNtupleDColumn(0, 3, KinE/MeV);
analysisManager->FillNtupleDColumn(0, 4, position.x()/cm);
analysisManager->FillNtupleDColumn(0, 5, position.y()/cm);
analysisManager->FillNtupleDColumn(0, 6, position.z()/cm);
analysisManager->FillNtupleDColumn(0, 7, Momentum.x());
analysisManager->FillNtupleDColumn(0, 8, Momentum.y());
analysisManager->FillNtupleDColumn(0, 9, Momentum.z());
analysisManager->FillNtupleDColumn(0, 10, time/s);
analysisManager->FillNtupleDColumn(0, 11, weight);
analysisManager->FillNtupleDColumn(0, 12, pol.x());
analysisManager->FillNtupleDColumn(0, 13, pol.y());
analysisManager->FillNtupleDColumn(0, 14, pol.z());
analysisManager->FillNtupleDColumn(0, 15, VPosx/cm);
analysisManager->FillNtupleDColumn(0, 16, VPosy/cm);
analysisManager->FillNtupleDColumn(0, 17, VPosz/cm);
analysisManager->FillNtupleDColumn(0, 18, VMomentum.x());
analysisManager->FillNtupleDColumn(0, 19, VMomentum.y());
analysisManager->FillNtupleDColumn(0, 20, VMomentum.z());
analysisManager->FillNtupleDColumn(0, 21, VKinE/MeV);
analysisManager->FillNtupleSColumn(0, 22, CPName);
analysisManager->FillNtupleSColumn(0, 23, SDet1);
//analysisManager->FillNtupleSColumn(0, 24, SD2);

  analysisManager->AddNtupleRow(0);
}

}

  return true;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Solar_SD::EndOfEvent(G4HCofThisEvent* /*hce*/)
{}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
