#ifndef DetectorConstruction_H 
#define DetectorConstruction_H 1

#include "G4VUserDetectorConstruction.hh"

class G4LogicalVolume;
class G4VPhysicalVolume;

class DetectorConstruction : public G4VUserDetectorConstruction{

public:

    DetectorConstruction();

    ~DetectorConstruction();

    G4VPhysicalVolume* Construct();

    //G4VPhysicalVolume* ConstructSDandField();

private:

  // logical volume for world

  G4LogicalVolume* logical_world;

  // physical volume for world

  G4VPhysicalVolume* physical_world;

};

#endif
