#include "G4RunManager.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "SteppingAction.hh"
#include "G4UIterminal.hh"
#include "G4UImanager.hh"
#include "G4UItcsh.hh"
#include "G4VisExecutive.hh"
#include "SNSRootManager.hh"
#include "time.h"
#include "Randomize.hh"
/*
#include "G4ImportanceBiasing.hh"
#include "G4ParallelWorldPhysics.hh"
#include "G4GeometrySampler.hh"
#include "G4GeometryManager.hh"
#include "G4IStore.hh"
#include "G4VWeightWindowStore.hh"
#include "G4WeightWindowAlgorithm.hh"
*/

#include "QGSP_BERT.hh"

//#include "QGSP_BERT.hh"

using namespace std;

//Function prototypes
G4bool IsValidInput (int,                   //argc
                     char**,                //argv
                     int &,                 //run Number
                     int &,                  //Events Processed
                     G4String &);            // macro file for batch mode

int main(int argc, char **argv){
    long seed = time(0);
    G4Random::setTheSeed(seed);
    
    //long seed = time(0);
    //CLHEP::HepRandom::setTheSeed(seed);
    int runNumber = 0;
    int NEvents = 0;
    G4String macroName = "";
    //G4bool useParallel = false;
    
    //check if the input is valid
    if (!IsValidInput (argc, argv, runNumber, NEvents,macroName)) {
        return 0;
    }
    
    //cout <<useParallel << endl;
    
    //If the input is valid, set output file name
    G4String outFile = Form("sns_out_%04d.root", runNumber);
  
    G4RunManager* rm = new G4RunManager();
    
    DetectorConstruction* detector = new DetectorConstruction();
    rm->SetUserInitialization(detector);
    rm->SetUserInitialization(new QGSP_BERT);
    rm->Initialize();
    G4VUserPrimaryGeneratorAction* pga = new PrimaryGeneratorAction();
    
    rm->SetUserAction(pga);
    
    SNSRootManager* rootman = new SNSRootManager(outFile);
    RunAction* ra = new RunAction(rootman, runNumber, NEvents);
    EventAction* ea = new EventAction(ra,rootman);
    SteppingAction* sa = new SteppingAction(ra,ea,rootman);
    
    rm->SetUserAction(ra);
    rm->SetUserAction(ea);
    rm->SetUserAction(sa);

    
    G4VisManager* visManager = new G4VisExecutive();
    visManager->Initialize();
    
    G4UIsession* session = new G4UIterminal(new G4UItcsh);
    G4UImanager* UI = G4UImanager::GetUIpointer();
    
    if (argc == 4) {
        // batch mode
        G4String command = "/control/execute ";
        UI->ApplyCommand(command+macroName);
    }
    else {
        // interactive mode : define UI session
        UI->ApplyCommand("/control/execute vis.mac");
        session->SessionStart();
	delete session;
    }
      
    delete visManager;
    delete rm;

    return 0;
}

G4bool IsValidInput (int argc, char** argv, int &runNumber, int & NEvents, G4String &macroName){
    G4bool result = true;
    //useParallel = false;
    if (argc < 3 || argc > 4) {
        result = false;
        G4cout << "Error: Invalid Inputs specified " << endl;
        G4cout << "Usage: SNS [runNumber] [EventsProcessed][run macro] " << endl;
    }
    else {
        result = true;
        runNumber = atoi (argv[1]);
        NEvents  =  atoi (argv[2]);
        //useParallel =  atoi(argv[3]);
        if (argc == 4){
            macroName = argv[3];
        }
        
    }
    return result;
}
