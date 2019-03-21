//SNSImportanceDetectorConstruction.hh

#ifndef SNSIMPORTANCEDETECTOR_HH
#define SNSIMPORTANCEDETECTOR_HH

//---------------------------------------------------------------------------//


#include "G4VUserParallelWorld.hh"
#include "globals.hh"

#include <map>
#include <vector>
//#include "G4GeometryCell.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4IStore;
class G4VIStore;
class G4VWeightWindowStore;
class SNSParallelMessenger;


//---------------------------------------------------------------------------//

class SNSParallelDetector: public G4VUserParallelWorld
{
	
public:
    enum EMGGeometryImportanceConsts {kMinImpValue = 1};
    
	SNSParallelDetector(G4String);
	virtual ~SNSParallelDetector();
	
	// overloaded from G4VUserParallelWorld:
	virtual void Construct();
	//virtual void ConstructSD();
    
    //G4VIStore*
    void CreateImportanceStore(); //caller should delete it
    //G4VWeightWindowStore *CreateWeightWindowStore();//caller should delete it
    //see B02ImportanceDetectorConstruction
    
    //void CreateImportanceStore();
    void SetImportanceValueForRegion(G4VPhysicalVolume*,
                                     G4double = kMinImpValue, G4int aRepNum = 0);
    G4VPhysicalVolume* GetWorldVolume() const { return fParallelWorld; }
    G4VPhysicalVolume& GetConstWorldVolumeAddress() const { return *fParallelWorld; }
    void SetNumVolumes( G4int );
    void SetIsVariableStep ( G4bool);
	
private:
	G4VPhysicalVolume* fParallelWorld;
	std::map<G4VPhysicalVolume* , std::pair<G4double, G4int>* > fImportanceMap;
};
//
#endif