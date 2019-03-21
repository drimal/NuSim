#include "DetectorConstruction.hh"
#include "globals.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Torus.hh"
#include "G4UnitsTable.hh"
#include "G4NistManager.hh"
#include "G4ThreeVector.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "TMath.h"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

DetectorConstruction::DetectorConstruction(){}

DetectorConstruction::~DetectorConstruction(){}

G4VPhysicalVolume* DetectorConstruction::Construct(){

  G4bool checkOverlaps = false; //true;

  G4NistManager* nist = G4NistManager::Instance();

    // define colors for visulization purpose
    G4Colour WHITE (1.0,1.0,1.0);
    G4Colour GRAY (0.5,0.5,0.5);
    G4Colour BLACK (0.0,0.0,0.0);
    G4Colour RED   (1.0,0.0,0.0);
    G4Colour GREEN (0.0,1.0,0.0);
    G4Colour BLUE  (0.0,0.0,1.0);
    G4Colour CYAN  (0.0,1.0,1.0);
    G4Colour MAGENTA (1.0,0.0,1.0);
    G4Colour YELLOW (1.0,1.0,0.0);
    G4Colour ORANGE (1.0,0.5,0.0);

// world volume

  G4NistManager* man = G4NistManager::Instance();
  G4Material* Air = man->FindOrBuildMaterial("G4_AIR");


// vacuum as a material

  G4Material* Vac = nist->FindOrBuildMaterial("G4_Galactic");


/*
   G4double worldx = 200.*m;
   G4double worldy = 200.*m;
   G4double worldz = 200.*m;
   G4Box* solid_world = new G4Box("world_box", worldx, worldy, worldz);
*/

   G4double worldRmin = 0.*m;
   G4double worldRmax = 300.0*m;
   G4double worldSPhi = 0*degree;
   G4double worldDPhi = 360*degree;
   G4double worldSTheta = 0*degree;
   G4double worldDTheta = 360*degree;

   G4Sphere* solid_world = new G4Sphere("world_sphere",
                                           worldRmin,
                                           worldRmax,
                                           worldSPhi,
                                           worldDPhi,
                                           worldSTheta,
                                           worldDTheta);


  G4LogicalVolume* logical_world = new G4LogicalVolume(solid_world, Air, "logical_vol_world", 0,0,0);
  G4VPhysicalVolume* physical_world = new G4PVPlacement(0,G4ThreeVector(0,0,0),logical_world, "physical_vol_world", 0, false, 0);

logical_world->SetVisAttributes(G4VisAttributes::GetInvisible());

    G4double a,z,density; 
    G4int ncomponents, natoms;
    G4double fractionmass;
    G4double pressure;
    G4double temperature;

// simplified Hg target

    G4Material* Hg = new G4Material("Mercury", z=80., a=200.59*g/mole, density=13.534*g/cm3);
            
    G4double targetx = 199.5*mm;
    G4double targety = 52*mm;
    G4double targetz = 250*mm;
    
  G4Box* solid_target = new G4Box("target_box", targetx, targety, targetz);
  G4LogicalVolume* logical_target = new G4LogicalVolume(solid_target,Hg,"logical_vol_target", 0,0,0);

  G4VisAttributes* targetVis = new G4VisAttributes (G4Colour(1.,0.,0.));
  logical_target->SetVisAttributes(targetVis);

    G4double targetShiftZ = -5.*cm;

  G4VPhysicalVolume* physical_target = new G4PVPlacement(0,G4ThreeVector(0,0,targetShiftZ),logical_target,"physical_vol_target",logical_world,false,0);

// define steel as a material

   G4Element* Fe  = new G4Element("Iron", "Fe" , z= 26., a=56.*g/mole);
   G4Element* C  = new G4Element("Carbon", "C" , z= 6., a=12.*g/mole);

   G4Material* Steel = new G4Material("Steel",density = 7.8*g/cm3, ncomponents=2);

   Steel->AddElement(Fe, fractionmass=98.*perCent);
   Steel->AddElement(C, fractionmass=2.*perCent);

// steel casing of Hg target - 6 plates
// steel plate 1

    G4double SteelPlate1X = 204.5*mm;
    G4double SteelPlate1Y = 52*mm;
    G4double SteelPlate1Z = 2.5*mm;
    
  G4Box* solid_SteelPlate1 = new G4Box("SteelPlate1_box", SteelPlate1X, SteelPlate1Y, SteelPlate1Z);
  G4LogicalVolume* logical_SteelPlate1 = new G4LogicalVolume(solid_SteelPlate1,Steel,"logical_SteelPlate1", 0,0,0);

  G4VisAttributes* SteelPlate1Vis = new G4VisAttributes (G4Colour(0.,1.,0.));
  logical_SteelPlate1->SetVisAttributes(SteelPlate1Vis);
  
  G4double SteelPlate1shiftZ = 202.5*mm;

  G4VPhysicalVolume* physical_SteelPlate1 = new G4PVPlacement(0,G4ThreeVector(0,0,SteelPlate1shiftZ),logical_SteelPlate1,"physical_vol_SteelPlate1",logical_world,false,0);

// Steel plate 2

    G4double SteelPlate2X = 204.5*mm;
    G4double SteelPlate2Y = 52*mm;
    G4double SteelPlate2Z = 2.5*mm;
 
  G4Box* solid_SteelPlate2 = new G4Box("SteelPlate2_box", SteelPlate2X, SteelPlate2Y, SteelPlate2Z);
  G4LogicalVolume* logical_SteelPlate2 = new G4LogicalVolume(solid_SteelPlate2,Steel,"logical_SteelPlate2", 0,0,0);

  G4VisAttributes* SteelPlate2Vis = new G4VisAttributes (G4Colour(0.,1.,0.));
  logical_SteelPlate2->SetVisAttributes(SteelPlate2Vis);

  G4double SteelPlate2shiftZ = -302.5*mm;

  G4VPhysicalVolume* physical_SteelPlate2 = new G4PVPlacement(0,G4ThreeVector(0,0,SteelPlate2shiftZ),logical_SteelPlate2,"physical_vol_SteelPlate2",logical_world,false,0);

// Steel plate 3

    G4double SteelPlate3X = 2.5*mm;
    G4double SteelPlate3Y = 52*mm;
    G4double SteelPlate3Z = 250*mm;
    
  G4Box* solid_SteelPlate3 = new G4Box("SteelPlate3_box", SteelPlate3X, SteelPlate3Y, SteelPlate3Z);
  G4LogicalVolume* logical_SteelPlate3 = new G4LogicalVolume(solid_SteelPlate3,Steel,"logical_SteelPlate3", 0,0,0);

  G4VisAttributes* SteelPlate3Vis = new G4VisAttributes (G4Colour(0.,1.,0.));
  logical_SteelPlate3->SetVisAttributes(SteelPlate3Vis);

  G4double SteelPlate3shiftX = -202*mm;
  G4double SteelPlate3shiftZ = -50*mm;

  G4VPhysicalVolume* physical_SteelPlate3 = new G4PVPlacement(0,G4ThreeVector(SteelPlate3shiftX,0,SteelPlate3shiftZ),logical_SteelPlate3,"physical_vol_SteelPlate3",logical_world,false,0);

// Steel plate 4
    
    G4double SteelPlate4X = 2.5*mm;
    G4double SteelPlate4Y = 52*mm;
    G4double SteelPlate4Z = 250*mm;
    
  G4Box* solid_SteelPlate4 = new G4Box("SteelPlate4_box", SteelPlate4X, SteelPlate4Y, SteelPlate4Z);
  G4LogicalVolume* logical_SteelPlate4 = new G4LogicalVolume(solid_SteelPlate4,Steel,"logical_SteelPlate4", 0,0,0);

  G4VisAttributes* SteelPlate4Vis = new G4VisAttributes (G4Colour(0.,1.,0.));
  logical_SteelPlate4->SetVisAttributes(SteelPlate4Vis);

  G4double SteelPlate4shiftX = 202*mm;
  G4double SteelPlate4shiftZ = -50*mm;

  G4VPhysicalVolume* physical_SteelPlate4 = new G4PVPlacement(0,G4ThreeVector(SteelPlate4shiftX,0,SteelPlate4shiftZ),logical_SteelPlate4,"physical_vol_SteelPlate4",logical_world,false,0);

// Steel plate 5
    
    G4double SteelPlate5X = 204.5*mm;
    G4double SteelPlate5Y = 2.5*mm;
    G4double SteelPlate5Z = 255*mm;
    
  G4Box* solid_SteelPlate5 = new G4Box("SteelPlate5_box", SteelPlate5X, SteelPlate5Y, SteelPlate5Z);
  G4LogicalVolume* logical_SteelPlate5 = new G4LogicalVolume(solid_SteelPlate5,Steel,"logical_SteelPlate5", 0,0,0);

  G4VisAttributes* SteelPlate5Vis = new G4VisAttributes (G4Colour(0.,1.,0.));
  logical_SteelPlate5->SetVisAttributes(SteelPlate5Vis);

  G4double SteelPlate5shiftY = 54.5*mm;
  G4double SteelPlate5shiftZ = -50*mm;

  G4VPhysicalVolume* physical_SteelPlate5 = new G4PVPlacement(0,G4ThreeVector(0,SteelPlate5shiftY,SteelPlate5shiftZ),logical_SteelPlate5,"physical_vol_SteelPlate5",logical_world,false,0);

// Steel plate 6
    
    G4double SteelPlate6X = 204.5*mm;
    G4double SteelPlate6Y = 2.5*mm;
    G4double SteelPlate6Z = 255*mm;
    
  G4Box* solid_SteelPlate6 = new G4Box("SteelPlate6_box", SteelPlate6X, SteelPlate6Y, SteelPlate6Z);
  G4LogicalVolume* logical_SteelPlate6 = new G4LogicalVolume(solid_SteelPlate6,Steel,"logical_SteelPlate6", 0,0,0);

  G4VisAttributes* SteelPlate6Vis = new G4VisAttributes (G4Colour(0.,1.,0.));
  logical_SteelPlate6->SetVisAttributes(SteelPlate6Vis);

  G4double SteelPlate6shiftY = -54.5*mm;
  G4double SteelPlate6shiftZ = -50*mm;

  G4VPhysicalVolume* physical_SteelPlate6 = new G4PVPlacement(0,G4ThreeVector(0,SteelPlate6shiftY,SteelPlate6shiftZ),logical_SteelPlate6,"physical_vol_SteelPlate6",logical_world,false,0);


// define material for Be plugs - Be(90%)+D2O(10%)

   G4Element* Be  = new G4Element("Beryllium", "Be" , z= 4., a=8.*g/mole);
   G4Element* D  = new G4Element("Deuterium", "D" , z=1., a=2.*g/mole);
   G4Element* O  = new G4Element("Oxygen", "O" , z= 8., a=16.*g/mole);

   G4Material* D2O = new G4Material("D2O",density = 1.1*g/cm3, ncomponents=2);

   D2O->AddElement(D, natoms=2);
   D2O->AddElement(O, natoms=1);
 
   G4Material* BeD2O = new G4Material("BeD2O",density = 1.85*g/cm3, ncomponents=2);
   
   BeD2O->AddElement(Be, fractionmass=90.*perCent);
   BeD2O->AddMaterial(D2O, fractionmass=10.*perCent);

// upper beryllium plug - A

   G4double UpBePlugRmin = 0.;
   G4double UpBePlugRmax = 35.*cm;
   G4double UpBePlugH = 22.5*cm;
   G4double UpBePlugPhimin = 0.*degree;
   G4double UpBePlugPhimax = 360.*degree; 
   
   G4Tubs* solid_UpBePlug = new G4Tubs ("Upper_Be_Plug", 
                                        UpBePlugRmin, 
                                        UpBePlugRmax, 
                                        UpBePlugH, 
                                        UpBePlugPhimin, 
                                        UpBePlugPhimax);
   G4LogicalVolume* logical_UpBePlug = new G4LogicalVolume(solid_UpBePlug,BeD2O,"logical_UpBePlug", 0,0,0);

  G4VisAttributes* UpBePlugVis = new G4VisAttributes (G4Colour(1.,1.,0.));
  logical_UpBePlug->SetVisAttributes(UpBePlugVis);

    G4RotationMatrix* rmUpBePlug = new G4RotationMatrix;
    rmUpBePlug -> rotateX(90*deg);

  G4double UpBePlugShiftY = 28.3*cm;

  G4VPhysicalVolume* physical_UpBePlug = new G4PVPlacement(rmUpBePlug,G4ThreeVector(0,UpBePlugShiftY,0),logical_UpBePlug,"physical_UpBePlug",logical_world,false,0);

// lower beryllium plug - A

   G4double DownBePlugRmin = 0.;
   G4double DownBePlugRmax = 35.*cm;
   G4double DownBePlugH = 22.5*cm;
   G4double DownBePlugPhimin = 0.*degree;
   G4double DownBePlugPhimax = 360.*degree;
      
   G4Tubs* solid_DownBePlug = new G4Tubs ("Down_Be_Plug", 
                                          DownBePlugRmin, 
                                          DownBePlugRmax, 
                                          DownBePlugH, 
                                          DownBePlugPhimin, 
                                          DownBePlugPhimax);
      
   G4LogicalVolume* logical_DownBePlug = new G4LogicalVolume(solid_DownBePlug,BeD2O,"logical_DownBePlug", 0,0,0);

  G4VisAttributes* DownBePlugVis = new G4VisAttributes (G4Colour(1.,1.,0.));
  logical_DownBePlug->SetVisAttributes(DownBePlugVis);

    G4RotationMatrix* rmDownBePlug = new G4RotationMatrix;
    rmDownBePlug -> rotateX(90*deg);

  G4double DownBePlugShiftY = -28.3*cm;

  G4VPhysicalVolume* physical_DownBePlug = new G4PVPlacement(rmDownBePlug,G4ThreeVector(0,DownBePlugShiftY,0),logical_DownBePlug,"physical_DownBePlug",logical_world,false,0);

// inner stainless steel(95%) + D2O(5%) cylinder

   G4Material* SteelD2O = new G4Material("SteelD2O",density = 7.8*g/cm3, ncomponents=2);
   
   SteelD2O->AddMaterial(Steel, fractionmass=95.*perCent);
   SteelD2O->AddMaterial(D2O, fractionmass=5.*perCent);

   G4double SteelD2ORmin = 35.1*cm;
   G4double SteelD2ORmax = 54.*cm;
   G4double SteelD2OH = 50.8*cm;
   G4double SteelD2OPhimin = 0.*degree;
   G4double SteelD2OPhimax = 360.*degree;
      
   G4Tubs* solid_SteelD2O = new G4Tubs ("SteelD2O", 
                                        SteelD2ORmin, 
                                        SteelD2ORmax, 
                                        SteelD2OH, 
                                        SteelD2OPhimin, 
                                        SteelD2OPhimax);

   G4LogicalVolume* logical_SteelD2O = new G4LogicalVolume(solid_SteelD2O,SteelD2O,"logical_SteelD2O", 0,0,0);

  G4VisAttributes* SteelD2OVis = new G4VisAttributes (G4Colour(1., 0., 1.));
  logical_SteelD2O->SetVisAttributes(SteelD2OVis);

    G4RotationMatrix* rmSteelD2O = new G4RotationMatrix;
    rmSteelD2O -> rotateX(90*deg);

  G4VPhysicalVolume* physical_SteelD2O = new G4PVPlacement(rmSteelD2O,G4ThreeVector(0,0,0),logical_SteelD2O,"physical_SteelD2O",logical_world,false,0);

// outer stainless steel(100%) reflector - C

   G4double SteelRmin = 54.1*cm;
   G4double SteelRmax = 510.*cm;
   G4double SteelH = 50.8*cm;
   G4double SteelPhimin = 0.*degree;
   G4double SteelPhimax = 360.*degree;

  G4Tubs* solid_Steel = new G4Tubs ("SteelOuter", 
                                    SteelRmin, 
                                    SteelRmax, 
                                    SteelH, 
                                    SteelPhimin, 
                                    SteelPhimax);

  G4LogicalVolume* logical_Steel = new G4LogicalVolume(solid_Steel,Steel,"logical_Steel", 0,0,0);

  G4VisAttributes* SteelVis = new G4VisAttributes (G4Colour(1., 0.0, 1.0));

  logical_Steel->SetVisAttributes(SteelVis);

    G4RotationMatrix* rmSteel = new G4RotationMatrix;
    rmSteel -> rotateX(90*deg);

  G4VPhysicalVolume* physical_Steel = new G4PVPlacement(rmSteel,G4ThreeVector(0,0,0),logical_Steel,"physical_Steel",logical_world,false,0);

// upper steel plate - C

   G4double UpSteelRmin = 0.*cm;
   G4double UpSteelRmax = 510.*cm;
   G4double UpSteelH = 175.*cm;
   G4double UpSteelPhimin = 0.*degree;
   G4double UpSteelPhimax = 360.*degree;

  G4Tubs* solid_UpSteel = new G4Tubs ("UpSteelOuter", 
                                      UpSteelRmin, 
                                      UpSteelRmax, 
                                      UpSteelH, 
                                      UpSteelPhimin, 
                                      UpSteelPhimax);

  G4LogicalVolume* logical_UpSteel = new G4LogicalVolume(solid_UpSteel,Steel,"logical_UpSteel", 0,0,0);

  G4VisAttributes* UpSteelVis = new G4VisAttributes (G4Colour(1.,1.,0.));
  logical_UpSteel->SetVisAttributes(UpSteelVis);

    G4RotationMatrix* rmUpSteel = new G4RotationMatrix;
    rmUpSteel -> rotateX(90*deg);
  
  G4double UpSteelShiftY = 275.8*cm;

  G4VPhysicalVolume* physical_UpSteel = new G4PVPlacement(rmUpSteel,G4ThreeVector(0,UpSteelShiftY,0),logical_UpSteel,"physical_UpSteel",logical_world,false,0);

// down steel plate - C

   G4double DownSteelRmin = 0*cm;
   G4double DownSteelRmax = 510.*cm;
   G4double DownSteelH = 175*cm;
   G4double DownSteelPhimin = 0.*degree;
   G4double DownSteelPhimax = 360.*degree;

  G4Tubs* solid_DownSteel = new G4Tubs ("DownSteelOuter", 
                                        DownSteelRmin, 
                                        DownSteelRmax, 
                                        DownSteelH, 
                                        DownSteelPhimin, 
                                        DownSteelPhimax);

  G4LogicalVolume* logical_DownSteel = new G4LogicalVolume(solid_DownSteel,Steel,"logical_DownSteel", 0,0,0);

  G4VisAttributes* DownSteelVis = new G4VisAttributes (G4Colour(1.,1.,0.));
  logical_DownSteel->SetVisAttributes(DownSteelVis);

    G4RotationMatrix* rmDownSteel = new G4RotationMatrix;
    rmDownSteel -> rotateX(90*deg);
  
  G4double DownSteelShiftY = -275.8*cm;

  G4VPhysicalVolume* physical_DownSteel = new G4PVPlacement(rmDownSteel,G4ThreeVector(0,DownSteelShiftY,0),logical_DownSteel,"physical_DownSteel",logical_world,false,0);


// circular room (made of concrete)

// concrete as a material

   G4Element* Ca  = new G4Element("Calcium","Ca", z= 20., a=40.*g/mole);
   G4Element* Si  = new G4Element("Silicon","Si", z= 14., a=28.*g/mole);
   G4Element* H = new G4Element("Hydrogen","H", z=1., a=1.*g/mole);
   G4Element* Al = new G4Element("Aluminum","Al", z= 13., a=27.*g/mole);
   G4Element* Mg = new G4Element("Magnesium","Mg", z=12., a=24.*g/mole);
   G4Element* S = new G4Element("Sulfur","S", z=16., a=32.*g/mole);
   G4Element* Mn = new G4Element("Manganese","Mn", z=25., a=55.*g/mole);
   G4Element* K = new G4Element("Potassium","K", z=19., a=39.*g/mole);
   G4Element* Na = new G4Element("Sodium","Na", z=11., a=23.*g/mole);
     
   G4Material* cement = new G4Material("Cement", density=3.9*g/cm3, ncomponents=11);
   cement->AddElement(O, fractionmass=55.714*perCent);
   cement->AddElement(Fe, fractionmass=28.026*perCent);
   cement->AddElement(H, fractionmass=10.293*perCent);
   cement->AddElement(Ca , fractionmass=3.09*perCent);
   cement->AddElement(Si , fractionmass=1.929*perCent);
   cement->AddElement(Al, fractionmass=0.55*perCent);
   cement->AddElement(Mg, fractionmass=0.236*perCent);
   cement->AddElement(S, fractionmass=0.092*perCent);
   cement->AddElement(Mn, fractionmass=0.044*perCent);
   cement->AddElement(K, fractionmass=0.014*perCent);
   cement->AddElement(Na, fractionmass=0.012*perCent);
   
 
   G4double RoomRmin = 510.1*cm;
   G4double RoomRmax = 640*cm;
   G4double RoomH = 500*cm;
   G4double RoomPhimin = 0.*degree;
   G4double RoomPhimax = 360.*degree;

  G4Tubs* solid_Room = new G4Tubs ("Room", 
                                   RoomRmin, 
                                   RoomRmax, 
                                   RoomH,
                                   RoomPhimin, 
                                   RoomPhimax);

  G4LogicalVolume* logical_Room = new G4LogicalVolume(solid_Room, cement,"logical_Room", 0,0,0);

  G4VisAttributes* RoomVis = new G4VisAttributes (G4Colour(0.,1.,1.));
  logical_Room->SetVisAttributes(RoomVis);

    G4RotationMatrix* rmRoom = new G4RotationMatrix;
    rmRoom -> rotateX(90*deg);

  G4double ConcreteRoomShiftY = 30.0*cm;

  G4VPhysicalVolume* physical_Room = new G4PVPlacement(rmRoom,G4ThreeVector(0,ConcreteRoomShiftY,0),logical_Room,"physical_Room",logical_world,false,0);

// steel shielding around the proton beamline

//plate 1

 G4double ProtonSteelPlate1X = 1.25*m; 
 G4double ProtonSteelPlate1Y = 1.25*m;
 G4double ProtonSteelPlate1Z = 10*m;


  G4Box* solid_ProtonSteelPlate1 = new G4Box("ProtonSteelPlate1_box", 
                                             ProtonSteelPlate1X, 
                                             ProtonSteelPlate1Y, 
                                             ProtonSteelPlate1Z);
  G4LogicalVolume* logical_ProtonSteelPlate1 = new G4LogicalVolume(solid_ProtonSteelPlate1,Steel,"logical_ProtonSteelPlate1", 0,0,0);

  G4VisAttributes* ProtonSteelPlate1Vis = new G4VisAttributes (G4Colour(1.,0.,0.));
  logical_ProtonSteelPlate1->SetVisAttributes(ProtonSteelPlate1Vis);

  G4double ProtonSteelPlate1shiftX = -1.75*m;
  G4double ProtonSteelPlate1shiftZ = 16.5*m;

  G4VPhysicalVolume* physical_ProtonSteelPlate1 = new
  G4PVPlacement(0,G4ThreeVector(ProtonSteelPlate1shiftX,0,ProtonSteelPlate1shiftZ),logical_ProtonSteelPlate1,"physical_vol_ProtonSteelPlate1",logical_world,false,0);

//plate 2

 G4double ProtonSteelPlate2X = 1.25*m; 
 G4double ProtonSteelPlate2Y = 1.25*m;
 G4double ProtonSteelPlate2Z = 10*m;


  G4Box* solid_ProtonSteelPlate2 = new G4Box("ProtonSteelPlate2_box", 
                                             ProtonSteelPlate2X, 
                                             ProtonSteelPlate2Y, 
                                             ProtonSteelPlate2Z);
  G4LogicalVolume* logical_ProtonSteelPlate2 = new G4LogicalVolume(solid_ProtonSteelPlate2,Steel,"logical_ProtonSteelPlate2", 0,0,0);

  G4VisAttributes* ProtonSteelPlate2Vis = new G4VisAttributes (G4Colour(1.,0.,0.));
  logical_ProtonSteelPlate2->SetVisAttributes(ProtonSteelPlate2Vis);

  G4double ProtonSteelPlate2shiftX = 1.75*m;
  G4double ProtonSteelPlate2shiftZ = 16.5*m;

  G4VPhysicalVolume* physical_ProtonSteelPlate2 = new
  G4PVPlacement(0,G4ThreeVector(ProtonSteelPlate2shiftX,0,ProtonSteelPlate2shiftZ),logical_ProtonSteelPlate2,"physical_vol_ProtonSteelPlate2",logical_world,false,0);

//plate 3

 G4double ProtonSteelPlate3X = 3.*m; 
 G4double ProtonSteelPlate3Y = 1.25*m;
 G4double ProtonSteelPlate3Z = 10*m;


  G4Box* solid_ProtonSteelPlate3 = new G4Box("ProtonSteelPlate3_box", 
                                             ProtonSteelPlate3X, 
                                             ProtonSteelPlate3Y, 
                                             ProtonSteelPlate3Z);
  G4LogicalVolume* logical_ProtonSteelPlate3 = new G4LogicalVolume(solid_ProtonSteelPlate3,Steel,"logical_ProtonSteelPlate3", 0,0,0);

  G4VisAttributes* ProtonSteelPlate3Vis = new G4VisAttributes (G4Colour(1.,0.,0.));
  logical_ProtonSteelPlate3->SetVisAttributes(ProtonSteelPlate3Vis);

  G4double ProtonSteelPlate3shiftY = 2.51*m;
  G4double ProtonSteelPlate3shiftZ = 16.5*m;

  G4VPhysicalVolume* physical_ProtonSteelPlate3 = new
  G4PVPlacement(0,G4ThreeVector(0,ProtonSteelPlate3shiftY,ProtonSteelPlate3shiftZ),logical_ProtonSteelPlate3,"physical_vol_ProtonSteelPlate3",logical_world,false,0);

//plate 4

 G4double ProtonSteelPlate4X = 3.*m; 
 G4double ProtonSteelPlate4Y = 1.25*m;
 G4double ProtonSteelPlate4Z = 10*m;


  G4Box* solid_ProtonSteelPlate4 = new G4Box("ProtonSteelPlate4_box", 
                                             ProtonSteelPlate4X, 
                                             ProtonSteelPlate4Y, 
                                             ProtonSteelPlate4Z);
  G4LogicalVolume* logical_ProtonSteelPlate4 = new G4LogicalVolume(solid_ProtonSteelPlate4,Steel,"logical_ProtonSteelPlate4", 0,0,0);

  G4VisAttributes* ProtonSteelPlate4Vis = new G4VisAttributes (G4Colour(1.,0.,0.));
  logical_ProtonSteelPlate4->SetVisAttributes(ProtonSteelPlate4Vis);

  G4double ProtonSteelPlate4shiftY = -2.51*m;
  G4double ProtonSteelPlate4shiftZ = 16.5*m;

  G4VPhysicalVolume* physical_ProtonSteelPlate4 = new
  G4PVPlacement(0,G4ThreeVector(0,ProtonSteelPlate4shiftY,ProtonSteelPlate4shiftZ),logical_ProtonSteelPlate4,"physical_vol_ProtonSteelPlate4",logical_world,false,0);


// -----------------------------------------------------here are where the geometry modifications were written in, M. McIntyre 7/19/2013--------------------------------------------------------------//

//define liquid hydrogen as a material

//G4Element* H = new G4Element("Hydrogen", "H", z=1, a= 1.01*g/mole);

G4Material* LH2 = new G4Material("LiquidHydrogen", density=67.8*mg/cm3, ncomponents=1, kStateLiquid, temperature = 20*kelvin, pressure=1.*atmosphere);

LH2 ->AddElement(H, natoms=1);




// liquid hydrogen moderator 1

G4double LiquidHydrogenModeratorRmin = 0.*cm ;
G4double LiquidHydrogenModeratorRmax = 5.*cm;
G4double LiquidHydrogenModeratorH = 18.*cm;
G4double LiquidHydrogenModeratorPhimin = 0.*degree;
G4double LiquidHydrogenModeratorPhimax = 360.*degree ;

G4Tubs* liquid_hydrogen1 = new G4Tubs ("LiquidHydrogenModerator1", 
                                       LiquidHydrogenModeratorRmin, 
                                       LiquidHydrogenModeratorRmax, 
                                       LiquidHydrogenModeratorH, 
                                       LiquidHydrogenModeratorPhimin, 
                                       LiquidHydrogenModeratorPhimax);

G4LogicalVolume* logical_hydrogen1 = new G4LogicalVolume(liquid_hydrogen1, LH2, "logical_hydrogen1", 0, 0, 0);

G4VisAttributes* LiquidHydrogenModeratorVis1 = new G4VisAttributes(G4Colour(1.,0.,0.));
logical_hydrogen1->SetVisAttributes(LiquidHydrogenModeratorVis1);

G4RotationMatrix* rmLHMod1 = new G4RotationMatrix;
    rmLHMod1 -> rotateX(90*deg);

G4double LiquidHydrogenModeratorShiftX1 = 0.*cm;
G4double LiquidHydrogenModeratorShiftY1 = -74.*cm;
G4double LiquidHydrogenModeratorShiftZ1 = 20.*cm;


G4VPhysicalVolume* physical_hydrogen1 = new G4PVPlacement(rmLHMod1, 
                                                          G4ThreeVector(LiquidHydrogenModeratorShiftX1, LiquidHydrogenModeratorShiftY1, LiquidHydrogenModeratorShiftZ1), 
                                                          logical_hydrogen1, 
                                                          "physical_vol_LiquidHydrogenModerator1",
                                                          logical_world, 
                                                          false, 
                                                          0);


// liquid hydrogen moderator 2

G4Tubs* liquid_hydrogen2 = new G4Tubs ("LiquidHydrogenModerator2", 
                                       LiquidHydrogenModeratorRmin, 
                                       LiquidHydrogenModeratorRmax, 
                                       LiquidHydrogenModeratorH, 
                                       LiquidHydrogenModeratorPhimin, 
                                       LiquidHydrogenModeratorPhimax);

G4LogicalVolume* logical_hydrogen2 = new G4LogicalVolume(liquid_hydrogen2, LH2, "logical_hydrogen2", 0, 0, 0);

G4VisAttributes* LiquidHydrogenModeratorVis2 = new G4VisAttributes(G4Colour(1.,0.,0.));
logical_hydrogen2->SetVisAttributes(LiquidHydrogenModeratorVis2);

G4RotationMatrix* rmLHMod2 = new G4RotationMatrix;
    rmLHMod2 -> rotateX(90*deg);

G4double LiquidHydrogenModeratorShiftX2 = 0.*cm;
G4double LiquidHydrogenModeratorShiftY2 = -74.*cm;
G4double LiquidHydrogenModeratorShiftZ2 = -20.*cm;


G4VPhysicalVolume* physical_hydrogen2 = new G4PVPlacement(rmLHMod2, 
                                                          G4ThreeVector(LiquidHydrogenModeratorShiftX2, LiquidHydrogenModeratorShiftY2, LiquidHydrogenModeratorShiftZ2), 
                                                          logical_hydrogen2, 
                                                          "physical_vol_LiquidHydrogenModerator2",
                                                          logical_world, 
                                                          false, 
                                                          0);


// liquid hydrogen moderator 3

G4Tubs* liquid_hydrogen3 = new G4Tubs ("LiquidHydrogenModerator3", 
                                       LiquidHydrogenModeratorRmin, 
                                       LiquidHydrogenModeratorRmax, 
                                       LiquidHydrogenModeratorH, 
                                       LiquidHydrogenModeratorPhimin, 
                                       LiquidHydrogenModeratorPhimax);

G4LogicalVolume* logical_hydrogen3 = new G4LogicalVolume(liquid_hydrogen3, LH2, "logical_hydrogen3", 0, 0, 0);

G4VisAttributes* LiquidHydrogenModeratorVis3 = new G4VisAttributes(G4Colour(1.,0.,0.));
logical_hydrogen3->SetVisAttributes(LiquidHydrogenModeratorVis3);

G4RotationMatrix* rmLHMod3 = new G4RotationMatrix;
    rmLHMod3 -> rotateX(90*deg);

G4double LiquidHydrogenModeratorShiftX3 = 0.*cm;
G4double LiquidHydrogenModeratorShiftY3 = 74.*cm;
G4double LiquidHydrogenModeratorShiftZ3 = -20.*cm;

G4VPhysicalVolume* physical_hydrogen3 = new G4PVPlacement(rmLHMod3, 
                                                          G4ThreeVector(LiquidHydrogenModeratorShiftX3, LiquidHydrogenModeratorShiftY3, LiquidHydrogenModeratorShiftZ3), 
                                                          logical_hydrogen3, 
                                                          "physical_vol_LiquidHydrogenModerator3",
                                                          logical_world, 
                                                          false, 
                                                          0);


// defining liquid water as a material

//G4double z, a, density;
//G4String name, symbol;
//G4int ncomponents, natoms;

//G4Element* H = new G4Element(name="Hydrogen", symbol="H", z=1., a = 1.01*g/mole);

//G4Element* O = new G4Element(name="Oxygen", symbol="O", z= 8., a = 16.00*g/mole);

G4Material* H2O = nist->FindOrBuildMaterial("G4_WATER");


// liquid water moderator

G4double LiquidH2ORmin = 0.*cm ;
G4double LiquidH2ORmax = 5.*cm;
G4double LiquidH2OH = 18.*cm;
G4double LiquidH2OPhimin = 0.*degree;
G4double LiquidH2OPhimax = 360.*degree ;

G4Tubs* liquid_water = new G4Tubs ("LiquidH2O",
                                   LiquidH2ORmin, 
                                   LiquidH2ORmax, 
                                   LiquidH2OH, 
                                   LiquidH2OPhimin, 
                                   LiquidH2OPhimax);

G4LogicalVolume* logical_water = new G4LogicalVolume(liquid_water, H2O, "logical_water", 0, 0, 0);

G4VisAttributes* LiquidH2OVis = new G4VisAttributes(G4Colour(1.,0.,0.));
logical_water->SetVisAttributes(LiquidH2OVis );

G4RotationMatrix* rmH2O = new G4RotationMatrix;
    rmH2O -> rotateX(90*deg);

G4double LiquidH2OShiftX = 0.*cm;
G4double LiquidH2OShiftY = 74.*cm;
G4double LiquidH2OShiftZ = 20.*cm;

G4VPhysicalVolume* physical_water = new G4PVPlacement(rmH2O, 
                                                      G4ThreeVector(LiquidH2OShiftX,LiquidH2OShiftY, LiquidH2OShiftZ), 
                                                      logical_water, 
                                                      "physical_vol_LiquidH2O",
                                                      logical_world, 
                                                      false, 
                                                      0);


// concrete ground floor beneath target

   G4double GNDFloorRmin = 0.0*m;
   G4double GNDFloorRmax = 30.0*m;
   G4double GNDFloorH = 0.4572*m;
   G4double GNDFloorPhimin = 0.*degree;
   G4double GNDFloorPhimax = 360.*degree;

	  G4RotationMatrix* rotCenter = new G4RotationMatrix;	
    rotCenter -> rotateX(90*deg);

    G4double gnd_flr_x_disp = 0*m;
    G4double gnd_flr_y_disp = -(RoomH + GNDFloorH)+ConcreteRoomShiftY;
    G4double gnd_flr_z_disp = 0*m;
      
   G4Tubs* solid_gnd_flr = new G4Tubs ("solid_gnd_flr", 
                                        GNDFloorRmin, 
                                        GNDFloorRmax, 
                                        GNDFloorH, 
                                        GNDFloorPhimin, 
                                        GNDFloorPhimax);

  G4ThreeVector gnd_flr_pos = G4ThreeVector(gnd_flr_x_disp, gnd_flr_y_disp, gnd_flr_z_disp);

    
  //G4Box* solid_gnd_flr = new G4Box("gnd_flr_box", gnd_flr_x, gnd_flr_y, gnd_flr_z);
  G4LogicalVolume* gnd_flr_logicV = new G4LogicalVolume(solid_gnd_flr,cement,"logical_vol_gnd_flr", 0,0,0);

  G4VPhysicalVolume* gnd_flr_physV = new G4PVPlacement(rotCenter,
                                                       gnd_flr_pos,
                                                       gnd_flr_logicV,
                                                      "gnd_flr_physV",
                                                       logical_world,
                                                       false,
                                                       0,
                                                       checkOverlaps);

  // basement concrete floor
    G4double base_flr_x_disp = 0*m;
    G4double base_flr_y_disp = -(RoomH + GNDFloorH + GNDFloorH + 5.6388*m);
    G4double base_flr_z_disp = 0*m;

  G4ThreeVector base_flr_pos = G4ThreeVector(base_flr_x_disp, base_flr_y_disp, base_flr_z_disp);

  G4LogicalVolume* base_flr_logicV = new G4LogicalVolume(solid_gnd_flr,cement,"logical_vol_base_flr", 0,0,0);

  G4VPhysicalVolume* base_flr_physV = new G4PVPlacement(rotCenter,
                                                       base_flr_pos,
                                                       base_flr_logicV,
                                                      "base_flr_physV",
                                                       logical_world,
                                                       false,
                                                       0,
                                                       checkOverlaps);

// steel monolith
/*
    G4double monoSteelRmax    = 5.*m;
    G4double monoSteelRmin    = 1.*m;
    G4double monoSteelHalfHt  = 6.25*m; //12.5/2 
    G4double protonBeamCutXZ  = 0.3*m; //42/2 cm is target in transverse 
    G4double transMonoZ       = 1.75*m;//6.25+1.75=8m bove; 6.25-1.75=4.5mbelow
    //minumum 1m ht below target is needed for innerSteelD2Ocyl2
    G4double monoLidTopHalfHt    = 0.5*m; // full ht=1m 
    G4double monoLidBotHalfHt    = 15.24*cm; //6.*inch; // 1 foot/2 

    G4Tubs* solid_steel_monolith = new G4Tubs ("steel_monolith",
                                               monoSteelRmin, 
                                               monoSteelRmax, 
                                               monoSteelHalfHt, 
                                               0, 
                                               twopi);

   G4LogicalVolume* steel_monolith_logicV = new G4LogicalVolume(solid_steel_monolith, Steel,"logical_vol_steel_monolith", 0,0,0);

   G4VPhysicalVolume* steel_monolith_physV = new G4PVPlacement(rotCenter,
                                                               G4ThreeVector(),
                                                               steel_monolith_logicV,
                                                               "steel_monolith_physV",
                                                               logical_world,
                                                               false,
                                                               0,
                                                               checkOverlaps);
*/
    

    // detector 1 (OscSNS far detector)

    G4double det_1_innerRadiusOfTheTube = 0.*cm;
    G4double det_1_outerRadiusOfTheTube = 4.0*m;
    G4double det_1_hightOfTheTube = 10.25*m;
    G4double det_1_startAngleOfTheTube = 0.*deg;
    G4double det_1_spanningAngleOfTheTube = 360.*deg;

    G4RotationMatrix* det1Rot = new G4RotationMatrix;
    //det1Rot.rotateX(150*deg);
    det1Rot->rotateY(-150*deg);
   //det1Rot->rotateZ(90*deg);


    G4double det_1_displace = 50*m + det_1_hightOfTheTube;
    G4double det_1_angle = -150*pi/180;

    G4double det_1_x = std::sin(det_1_angle)*det_1_displace;
    G4double det_1_y = 0.0*m;
    G4double det_1_z = -std::cos(det_1_angle)*det_1_displace;

    G4ThreeVector det_1_pos = G4ThreeVector(det_1_x, det_1_y, det_1_z);

    //G4Transform3D det_1_transform = G4Transform3D(det1Rot, det_1_pos);

    G4Tubs* det_1_solid = new G4Tubs("det_1_solid",
                                   det_1_innerRadiusOfTheTube, 
                                   det_1_outerRadiusOfTheTube,
                                   det_1_hightOfTheTube,
                                   det_1_startAngleOfTheTube, 
                                   det_1_spanningAngleOfTheTube);

    G4LogicalVolume* det_1_logicV = new G4LogicalVolume(det_1_solid, Vac, "det_1_logicV", 0,0,0);

    G4VisAttributes* det_clr = new G4VisAttributes(MAGENTA);
    det_1_logicV -> SetVisAttributes(det_clr);

    G4VPhysicalVolume* det_1_physV = new G4PVPlacement(det1Rot,
                                                       det_1_pos,                                                
                                                       det_1_logicV,
                                                      "det_1_physV",
                                                       logical_world,
                                                       false,
                                                       0,
                                                       checkOverlaps);


    // detector 2 (near detector)

    G4double det_2_innerRadiusOfTheTube = 0.*cm;
    G4double det_2_outerRadiusOfTheTube = 2.5*m;
    G4double det_2_hightOfTheTube = 2.5*m;
    G4double det_2_startAngleOfTheTube = 0.*deg;
    G4double det_2_spanningAngleOfTheTube = 360.*deg;

    G4RotationMatrix* det2Rot = new G4RotationMatrix;
    //det2Rot.rotateX(150*deg);
    det2Rot->rotateY(90*deg);
    //det2Rot->rotateZ(90*deg);


    G4double det_2_displace = 10*m + 1*m + det_2_hightOfTheTube;
    G4double det_2_angle = -90*pi/180;

    G4double det_2_x = std::sin(det_2_angle)*det_2_displace;
    G4double det_2_y = 0.0*m;
    G4double det_2_z = -std::cos(det_2_angle)*det_2_displace; 


    G4ThreeVector det_2_pos = G4ThreeVector(det_2_x, det_2_y, det_2_z);

//    G4Transform3D det_2_transform = G4Transform3D(det2Rot, det_2_pos);

    G4Tubs* det_2_solid = new G4Tubs("det_2_solid",
                                      det_2_innerRadiusOfTheTube, 
                                      det_2_outerRadiusOfTheTube,
                                      det_2_hightOfTheTube,
                                      det_2_startAngleOfTheTube, 
                                      det_2_spanningAngleOfTheTube);

    G4LogicalVolume* det_2_logicV = new G4LogicalVolume(det_2_solid, Vac, "det_2_logicV", 0,0,0);

    det_2_logicV -> SetVisAttributes(det_clr);

    G4VPhysicalVolume* det_2_physV = new G4PVPlacement(det2Rot,
                                                       det_2_pos,
                                                       det_2_logicV,
                                                      "det_2_physV",
                                                       logical_world,
                                                       false,
                                                       0,
                                                       checkOverlaps);


    // detector 3 (basement detector)

    G4double det_3_innerRadiusOfTheTube = 0.*cm;
    G4double det_3_outerRadiusOfTheTube = 1.0*m;
    G4double det_3_hightOfTheTube = 2.5*m;
    G4double det_3_startAngleOfTheTube = 0.*deg;
    G4double det_3_spanningAngleOfTheTube = 360.*deg;

    G4RotationMatrix* det3Rot = new G4RotationMatrix;
    //det3Rot.rotateX(150*deg);
    det3Rot->rotateY(110*deg);
    //det3Rot->rotateZ(90*deg);


    G4double det_3_displace = 20*m + det_3_hightOfTheTube;
    G4double det_3_angle = 110*pi/180;

    G4double det_3_x = std::sin(det_3_angle)*det_3_displace;
    G4double det_3_y = -(RoomH + GNDFloorH + 5.6388*m)+det_3_outerRadiusOfTheTube; //vertical displacement underground
    G4double det_3_z = -std::cos(det_3_angle)*det_3_displace; 


    G4ThreeVector det_3_pos = G4ThreeVector(det_3_x, det_3_y, det_3_z);

//    G4Transform3D det_3_transform = G4Transform3D(det3Rot, det_3_pos);

    G4Tubs* det_3_solid = new G4Tubs("det_3_solid",
                                      det_3_innerRadiusOfTheTube, 
                                      det_3_outerRadiusOfTheTube,
                                      det_3_hightOfTheTube,
                                      det_3_startAngleOfTheTube, 
                                      det_3_spanningAngleOfTheTube);

    G4LogicalVolume* det_3_logicV = new G4LogicalVolume(det_3_solid, Vac, "det_3_logicV", 0,0,0);

    det_3_logicV -> SetVisAttributes(det_clr);

    G4VPhysicalVolume* det_3_physV = new G4PVPlacement(det3Rot,
                                                       det_3_pos,
                                                       det_3_logicV,
                                                      "det_3_physV",
                                                       logical_world,
                                                       false,
                                                       0,
                                                       checkOverlaps);
/*
    // detector 3 (basement detector)
    G4double det_3_Rmin = 0.*m;
    G4double det_3_Rmax = 2.0*m;
    G4double det_3_SPhi = 0*degree;
    G4double det_3_DPhi = 360*degree;
    G4double det_3_STheta = 0*degree;
    G4double det_3_DTheta = 360*degree;

    G4RotationMatrix* det3Rot = new G4RotationMatrix;
    //det3Rot.rotateX(150*deg);
    det3Rot->rotateY(110*deg);
    //det3Rot->rotateZ(90*deg);


    G4double det_3_displace = 20*m;
    G4double det_3_angle = 110*pi/180;

    G4double det_3_x = std::sin(det_3_angle)*det_3_displace;
    G4double det_3_y = -(RoomH + GNDFloorH + 5.6388*m) + det_3_Rmax; //vertical displacement underground
    G4double det_3_z = -std::cos(det_3_angle)*det_3_displace;


    G4ThreeVector det_3_pos = G4ThreeVector(det_3_x, det_3_y, det_3_z);

    //    G4Transform3D det_3_transform = G4Transform3D(det3Rot, det_3_pos);

    G4Sphere* det_3_solid = new G4Sphere("det_3_sphere",
                                         det_3_Rmin,
                                         det_3_Rmax,
                                         det_3_SPhi,
                                         det_3_DPhi,
                                         det_3_STheta,
                                         det_3_DTheta);


    G4LogicalVolume* det_3_logicV = new G4LogicalVolume(det_3_solid, Vac, "det_3_logicV", 0,0,0);

    det_3_logicV -> SetVisAttributes(det_clr);

    G4VPhysicalVolume* det_3_physV = new G4PVPlacement(det3Rot,
                                                       det_3_pos,
                                                       det_3_logicV,
                                                       "det_3_physV",
                                                       logical_world,
                                                       false,
                                                       0,
                                                       checkOverlaps);



*/


/*
    // detector 4 (forward detector, 20m, 20 deg)

    G4double det_4_innerRadiusOfTheTube = 0.*cm;
    G4double det_4_outerRadiusOfTheTube = 2.5*m;
    G4double det_4_hightOfTheTube = 2.5*m;
    G4double det_4_startAngleOfTheTube = 0.*deg;
    G4double det_4_spanningAngleOfTheTube = 360.*deg;

    G4RotationMatrix* det4Rot = new G4RotationMatrix;
    //det3Rot.rotateX(150*deg);
    det4Rot->rotateY(-20*deg);
    //det3Rot->rotateZ(90*deg);


    G4double det_4_displace = 20*m + 1*m + det_4_hightOfTheTube;
    G4double det_4_angle = -20*pi/180;

    G4double det_4_x = std::sin(det_4_angle)*det_4_displace;
    G4double det_4_y = 0.0*m;
    G4double det_4_z = -std::cos(det_4_angle)*det_4_displace; 


    G4ThreeVector det_4_pos = G4ThreeVector(det_4_x, det_4_y, det_4_z);

//    G4Transform3D det_3_transform = G4Transform3D(det3Rot, det_3_pos);

    G4Tubs* det_4_solid = new G4Tubs("det_4_solid",
                                      det_4_innerRadiusOfTheTube, 
                                      det_4_outerRadiusOfTheTube,
                                      det_4_hightOfTheTube,
                                      det_4_startAngleOfTheTube, 
                                      det_4_spanningAngleOfTheTube);

    G4LogicalVolume* det_4_logicV = new G4LogicalVolume(det_4_solid, Vac, "det_4_logicV", 0,0,0);

    det_4_logicV -> SetVisAttributes(det_clr);

    G4VPhysicalVolume* det_4_physV = new G4PVPlacement(det4Rot,
                                                       det_4_pos,
                                                       det_4_logicV,
                                                      "det_4_physV",
                                                       logical_world,
                                                       false,
                                                       0,
                                                       checkOverlaps);

*/

  
  
   return physical_world;
}
