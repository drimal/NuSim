#include "SNSParallelDetector.hh"
#include "globals.hh"
#include <sstream>
#include <map>
#include "G4ios.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4Material.hh"
#include "G4VisAttributes.hh"
#include "G4Color.hh"
#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4UserLimits.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4PVPlacement.hh"



// For Primitive Scorers
#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4SDParticleFilter.hh"
#include "G4PSNofCollision.hh"
#include "G4PSPopulation.hh"
#include "G4PSTrackCounter.hh"
#include "G4PSTrackLength.hh"
#include "G4VIStore.hh"
#include "G4IStore.hh"   // for importance biasing
#include "G4WeightWindowStore.hh"    // for weight window technique


SNSParallelDetector::SNSParallelDetector(G4String worldName)
:G4VUserParallelWorld(worldName),fParallelWorld(NULL){}

SNSParallelDetector::~SNSParallelDetector(){}

void SNSParallelDetector::Construct()
{
    
    
    fParallelWorld = GetWorld();
    G4cout << "Constructing Parallel World" << G4endl;
    
    //G4Material* dummyMat  = 0;
    
    G4cout << " SNSParallelDetector:: ParallelWorldName = "<< fParallelWorld->GetName() << G4endl;
    
    G4LogicalVolume* parallelWorldLogical = fParallelWorld->GetLogicalVolume();
    G4String shape = "cylinder";//"box"
    
    //Don't change the following keyword in comment, used in script
    G4int fNumVol                     = 10; //NUMVOL
    if(shape == "box"){
        G4double step = 5.0*m;
        G4double xHalf = 150.0*m;
        G4double yHalf = 150.0*m;
        G4double zHalf = 150.0*m;

        G4ThreeVector center(0,0,0);
        G4RotationMatrix* rot         = new G4RotationMatrix();
        G4LogicalVolume* parLog = NULL;
        G4PVPlacement *parallelPhys = NULL;
    
        for (int i=0; i < fNumVol; i++){
        
            G4cout << " outer: x:y:z " << xHalf-i*step << " " <<yHalf-i*step << " " << zHalf-i*step<< G4endl;
            std::ostringstream ss;
            ss << "slice"<< i;
            G4String sliceName = ss.str();
            G4Box* box1 = new G4Box("box1", xHalf-i*step, yHalf-i*step,zHalf-i*step);
            if(i < fNumVol){
                G4cout << " inner: x:y:z " << xHalf-(i+1)*step << " " <<yHalf-(i+1)*step << " " << zHalf-(i+1)*step<< G4endl;
            
                G4Box* box2 = new G4Box("box2", xHalf-(i+1)*step,yHalf-(i+1)*step,zHalf-(i+1)*step);
            
                G4SubtractionSolid* sub = new G4SubtractionSolid(sliceName, box1,box2,rot, center);
                parLog = new G4LogicalVolume(sub,0,sliceName+"Log");
            }
            else if(i == fNumVol){
                parLog = new G4LogicalVolume(box1,0,sliceName+"Log");
            }
        parLog->SetVisAttributes(new G4VisAttributes(G4Colour::Green()));
        parallelPhys = new G4PVPlacement(rot, center, parLog,// NOTE change pos
                                         "ParallelLog", parallelWorldLogical, false, 0);
        SetImportanceValueForRegion(parallelPhys, std::pow(2.0, i+1.0));
      }
    }
    else if(shape == "cylinder"){
          
          G4SubtractionSolid* cylSolid = NULL;
          G4Tubs* cylFull = NULL;
          G4Tubs* cylInner = NULL;
          G4double detHalfHt = 0.5*m;
          G4bool surfaceDet = false;
          detHalfHt = 5*m;
          surfaceDet = true;
          //use exampleRE06.cc to set meesenger
          G4bool cylSection          = false;
          G4double parCylHalfHt 		= 7.*m; //
          if(surfaceDet)parCylHalfHt = 11.5*m;//6 above tgt, 14 below
        G4double mainVolRad		    = 22.3*m;
          if(surfaceDet) mainVolRad  = 23.*m; //25 - detToWallHor
          G4double offBeamZ 			= 0.5*m;	  //protBoxHalfXZ= 2*m;
          G4double offTargetAxis		= 1.*m;
          G4double detToWallHor		= 1.5*m;
          G4double detToWallVert		= detHalfHt + 0.8*m;
          G4double detToTargetRad   	= mainVolRad + detToWallHor;
          
          G4double stepCylRad 		= (mainVolRad - offTargetAxis)/fNumVol-1;
          G4double stepCylHt  		= (parCylHalfHt - detToWallVert)/fNumVol-1;
          G4double posParCylZ        =  -parCylHalfHt -offBeamZ;
          if(surfaceDet) posParCylZ  =  -3.7*m;
          G4ThreeVector posParCyl(0,0, posParCylZ); //cylinder
          G4ThreeVector posParBox(18.2*m,-5.*m,-parCylHalfHt -offBeamZ);
          
          G4ThreeVector center(0,0,0);
          G4PVPlacement* parallelPhys = NULL;
          G4LogicalVolume* parLog 	= NULL;
          G4RotationMatrix* rot 	    = new G4RotationMatrix();
          
          G4double dA = 15.*degree/fNumVol; 
          G4double theta1 = 0, theta2 = 0, theta3 = 0, theta4 = 0;	
          
          theta1 = theta3 = 0;
          theta2 = theta4 = twopi;

          
          for(int i = 0; i < fNumVol; i++){
              std::ostringstream ss;
              ss << "parallelSNS"<< i;
              G4String pName = ss.str();

              if(cylSection) {
                  theta1 = 3./2.*pi+ i*dA;
                  theta2 = pi/2.-2.*i*dA;
                  theta3 = 3./2.*pi+ (i+1)*dA;
                  theta4 = pi/2.-2.*(i+1)*dA;
              }
              G4cout << " Full:  rad1 " << (i*stepCylRad+offTargetAxis)/m
              <<  " rad2 "<<(2.*detToTargetRad-offTargetAxis-i*stepCylRad)/m
              << " ht "<<(parCylHalfHt-i*stepCylHt)/m  << G4endl;
              
              
              cylFull = new G4Tubs("BasementCyl", i*stepCylRad+ offTargetAxis, 2.*detToTargetRad- offTargetAxis-i*stepCylRad,
                                 parCylHalfHt-i*stepCylHt, theta1, theta2 );
              if(i < fNumVol){
                  cylInner = new G4Tubs("CylIn",(i+1)*stepCylRad+offTargetAxis, 2.*detToTargetRad- offTargetAxis-(i+1)*stepCylRad,
                            parCylHalfHt-(i+1)*stepCylHt, theta3, theta4); //0,twopi
                  
                  G4cout << " Inner:  rad1 " << ((i+1)*stepCylRad+offTargetAxis)/m
                  <<  " rad2 "<<(2.*detToTargetRad-offTargetAxis-(i+1)*stepCylRad)/m
                  << " ht "<<(parCylHalfHt-(i+1)*stepCylHt)/m << G4endl;
                  
                  cylSolid =  new G4SubtractionSolid("SubCyl",cylFull, cylInner, rot, center);
                  parLog = new G4LogicalVolume(cylSolid,0,pName+"Log");
              }
              else if(i == fNumVol){
                parLog = new G4LogicalVolume(cylFull,0,pName+"Log");
              }
              parLog->SetVisAttributes(new G4VisAttributes(G4Colour::Cyan()));
              parallelPhys = new G4PVPlacement(rot, posParCyl, parLog,// NOTE change pos
                                               pName, parallelWorldLogical, false, 0);
              SetImportanceValueForRegion(parallelPhys, std::pow(2.0, i+1.0));
          }
    }
    
}
//B4Box* main = new G4Box("main",10*m,10*m,10*m);
//G4LogicalVolume* log = new G4LogicalVolume(main,mat, "mainlog");
//new G4VPlacement(rot,pos,log,"MainVolume",mainLogical,false,0);


//G4VIStore*
void SNSParallelDetector::CreateImportanceStore()
{
    G4cout << "SNSParallelDetector:: Creating Importance Store "<< G4endl;
    
    G4IStore *istore = G4IStore::GetInstance(GetName());
    
    G4GeometryCell parallelVolCell(GetConstWorldVolumeAddress(), 0);
    istore->AddImportanceGeometryCell(1, parallelVolCell);
    
    for(std::map<G4VPhysicalVolume*, std::pair<G4double, G4int>* >::iterator i =
        fImportanceMap.begin();i != fImportanceMap.end(); i++)
      {
        istore->AddImportanceGeometryCell((i->second)->first,
                                          *(i->first), (i->second)->second);
        
        G4cout <<" Volume added: "<< (i->first)->GetName() << " (Replica "
        <<(i->second)->second<< ") of imp value: "<< (i->second)->first <<G4endl;
      }
    // return istore;
}

//---------------------------------------------------------------------------//

void SNSParallelDetector::SetImportanceValueForRegion(G4VPhysicalVolume* aVol,  G4double anImpValue,   G4int aRepNum)
{
    // sets importance value (default kMinImpValue)
    // anImpValue: 0 -> particle will be killed
    //           : > 0 -> allowed (higher means greater importance)
    //           : < 0 -> not allowed!
    
    if (anImpValue< 0) G4cout << "Error in importance value..."<< G4endl;
    else
      {
        std::map<G4VPhysicalVolume*, std::pair<G4double, G4int>* > :: const_iterator anIt;
        anIt = fImportanceMap.find(aVol);
        if(anIt == fImportanceMap.end())
          {
            if(anImpValue < kMinImpValue && anImpValue > 0.) anImpValue = kMinImpValue;
            std::pair<G4double, G4int>* aPair =
            new std::pair<G4double, G4int>(anImpValue, aRepNum);
            fImportanceMap.insert(
                                  std::pair<G4VPhysicalVolume*, std::pair<G4double, G4int>* >(aVol,aPair));
            G4cout << "Set " << aVol->GetName() << " (: " << aRepNum
            << ") of importance: " << anImpValue << G4endl;
          }
        else
          {
            G4cout << "WARNING: Importance numberfor: " << aVol->GetName()
            << " had been set" << G4endl;
          }
      }
}
