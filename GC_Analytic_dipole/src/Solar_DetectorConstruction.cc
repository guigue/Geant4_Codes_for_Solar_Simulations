//
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


#include "Solar_DetectorConstruction.hh"
#include "Solar_SD.hh"
#include "Solar_FieldSetup.hh"
#include "Solar_Dipole.hh"


#include "G4NistManager.hh"
#include "G4Element.hh"
#include "G4Material.hh"

#include "G4FieldManager.hh"
#include "G4UniformMagField.hh"
#include "G4TransportationManager.hh"

#include "G4PropagatorInField.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4RunManager.hh"

#include "G4SDManager.hh"
#include "G4ChordFinder.hh"

#include "G4Box.hh"
//#include "G4Orb.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4RotationMatrix.hh"

#include "G4Trd.hh"
#include "G4SystemOfUnits.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4ProductionCuts.hh"
#include "G4FieldManager.hh"
#include "G4AutoDelete.hh"

#include "G4UserLimits.hh"

#include "G4ProductionCutsTable.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Solar_DetectorConstruction::Solar_DetectorConstruction()
: G4VUserDetectorConstruction(),
  fCheckOverlaps(true),
  fMagneticLogical(0)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Solar_DetectorConstruction::~Solar_DetectorConstruction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* Solar_DetectorConstruction::Construct()
{


//Materials
G4NistManager*man = G4NistManager::Instance();
man->SetVerbose(1) ;

G4Element*H=man->FindOrBuildElement("H");
G4Element*He=man->FindOrBuildElement("He");
G4Element*C=man->FindOrBuildElement("C");
G4Element*N=man->FindOrBuildElement("N");
G4Element*O=man->FindOrBuildElement("O");

density = 4.0E-7*g/cm3;
//density = 0.2*g/cm3;
G4Material* Solar = new G4Material("Solar",density,ncomponents=5);//, kStateGas, 6000*kelvin, 1.*atmosphere);
Solar->AddElement(H, fractionmass=92.055*perCent);
Solar->AddElement(He, fractionmass=7.825*perCent);
Solar->AddElement(C, fractionmass=0.0225*perCent);
Solar->AddElement(N, fractionmass=0.0555*perCent);
Solar->AddElement(O, fractionmass=0.0421*perCent);

densityCorona = 1.0E-13*g/cm3;
//density = 0.2*g/cm3;
G4Material* SolarCorona = new G4Material("SolarCorona",densityCorona,ncomponents=5);//, kStateGas, 6000*kelvin, 1.*atmosphere);
SolarCorona->AddElement(H, fractionmass=92.055*perCent);
SolarCorona->AddElement(He, fractionmass=7.825*perCent);
SolarCorona->AddElement(C, fractionmass=0.0225*perCent);
SolarCorona->AddElement(N, fractionmass=0.0555*perCent);
SolarCorona->AddElement(O, fractionmass=0.0421*perCent);


G4Material* Vac = man->FindOrBuildMaterial("G4_Galactic");

// For an atmosphere only with hydrogen

//G4Material* Hydrogen =new G4Material("Hydrogen", 8.3748E-05*g/cm3, ncomponents=1);
//Hydrogen->AddElement(H, fractionmass=100*perCent);

//G4Material* HydroCorona = new G4Material("HydroCorona",1., 1.008, densityCorona, kStateGas);
//G4Material* HydroAtm = new G4Material("HydroAtm",1., 1.008, density, kStateGas);


G4cout << *(G4Material::GetMaterialTable()) << G4endl;


G4double Rs = 6.96e5*km;     // Solar radius (km)
G4double Rf = 0.010*Rs;      // half of feet separation (km)
//G4double lz = 0.0267632*Rs;  // loop top z-coordinate (km)
G4double lh =  0.0152632*Rs;;  // loop height (km)
//G4double dz = lz - lh;       // dipoloe depth (km)


//Dimensions
G4double solar_x = 1.2*Rf;
G4double solar_y = 1.2*Rf;
G4double solar_z = 1.2*lh/2.0;
G4double pos_x = 0.0*cm;
G4double pos_y = 0.0*cm;
G4double pos_z = 0.0*cm;


//World
G4double world_hx = 1.2*Rf;
G4double world_hy = 1.2*Rf;
G4double world_hz = 1.2*lh;
WorldBox = new G4Box("World", world_hx, world_hy, world_hz);
WorldLog =new G4LogicalVolume(WorldBox, Vac, "WorldLog");
WorldPhys = new G4PVPlacement(0,
    G4ThreeVector(pos_x,pos_y,pos_z),
    WorldLog,
    "WorldPhys",
    0,
    false,
    0);

//Atmosphere
pos_z = -1.2*lh/2.0;
AtmBox = new G4Box("Atm", solar_x, solar_y, solar_z);

AtmLog
= new G4LogicalVolume(
  AtmBox,     // its solid
  Solar,      // its Material
  "AtmLog");  // its name

AtmPhys
= new G4PVPlacement(
  0,                    // no rotation
  G4ThreeVector(pos_x, pos_y, pos_z), // its position
    AtmLog,             // its logical volume
    "AtmPhys",          // its name
    WorldLog,           // its mother volume
    false,              // no boolean operation
    0);                 // copy number


//Define Atm Region

atmRegion = new G4Region("AtmRegion");
AtmLog->SetRegion(atmRegion);
atmRegion->AddRootLogicalVolume(AtmLog);

//Set different Cuts for this region. Uncomment lines below for different cuts.

//G4ProductionCuts*atmcuts = new G4ProductionCuts();
//atmcuts->SetProductionCut(1.0*cm, "gamma");
//atmcuts->SetProductionCut(1.0*cm, "e-");
//atmcuts->SetProductionCut(1.0*cm, "e+");
//atmRegion->SetProductionCuts(atmcuts);

G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(100*keV, 1*GeV);

//fMagneticLogical = AtmLog;

//Corona
CoronaBox = new G4Box("Corona", solar_x, solar_y, solar_z);

CoronaLog
=new G4LogicalVolume(
CoronaBox,    // its solid
SolarCorona,  // its material
"CoronaLog"); // its name

pos_z = 1.2*lh/2.0;
CoronaPhys
= new G4PVPlacement(
0,                    // no rotation
G4ThreeVector(pos_x, pos_y, pos_z), // its position
    CoronaLog,        // its logical volume
    "CoronaPhys",     // its name
    WorldLog,         // its mother volume
    false,            // no boolean operation
    0);               // copy number

//Define Corona Region --for example if we want different cuts for this region

corRegion = new G4Region("CorRegion");
CoronaLog->SetRegion(corRegion);
corRegion->AddRootLogicalVolume(CoronaLog);

//Set different Cuts for this region. Uncomment lines below for different cuts.

//G4ProductionCuts*CoronaCuts = new G4ProductionCuts();
//CoronaCuts->SetProductionCut(1.0*cm, "gamma");
//CoronaCuts->SetProductionCut(1.0*cm, "e-");
//CoronaCuts->SetProductionCut(1.0*cm, "e+");
//corRegion->SetProductionCuts(CoronaCuts);


//User Limits
//....................Stop mirroring particles at the Corona..........//

SolarUserLimits = new G4UserLimits();
//SolarUserLimits->SetUserMaxTime(0.5*s);
//SolarUserLimits->SetMaxAllowedStep(30*m);
SolarUserLimits->SetUserMaxTrackLength(100000*km);
CoronaLog->SetUserLimits(SolarUserLimits);

//SolarUserLimits = new G4UserLimits ();

//SolarUserLimits->SetUserMaxTime(TimeMax);
//CoronaLog->SetUserLimits(SolarUserLimits);

//CoronaLog->SetUserLimits(new G4UserLimits (TimeMax));
//AtmLog->SetUserLimits(new G4UserLimits (1*km, TimeMax, TrackMax));
//CoronaLog->SetUserMaxTime(TimeMax);



return WorldPhys;

}



void Solar_DetectorConstruction::ConstructSDandField()
{


  auto CoronaSD
  = new Solar_SD("CoronaSD",1);
  G4SDManager::GetSDMpointer()->AddNewDetector(CoronaSD);
  SetSensitiveDetector("CoronaLog",CoronaSD);

  auto AtmSD
  = new Solar_SD("AtmSD", 2);
  G4SDManager::GetSDMpointer()->AddNewDetector(AtmSD);
  SetSensitiveDetector("AtmLog",AtmSD);


  if (!fEmFieldSetup.Get()){
    Solar_FieldSetup* field = new Solar_FieldSetup(/*0.1*tesla*/);


    G4AutoDelete::Register(field); //Kernel will delete the messenger
    fEmFieldSetup.Put(field);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
