#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "RunAction.hh"
#include "globals.hh"
#include "SNSRootManager.hh"
#include <vector>
#include <set>
#include <iomanip>
#include <map>

typedef  struct part_t_def {
  part_t_def(): name(""),sec_name(""),proc(""),parent(""), ltime(0), gtime(0), KinEnergy(0), trackID(0), volFlag(0), partFlag(0), stepFlag(0), procFlag(0) {}
  G4String name, sec_name, proc, parent; 
  G4double ltime, gtime;    
  G4double KinEnergy; 
  G4int id, track_id;          
  G4int trackID;
  

  G4int volFlag, partFlag, stepFlag, procFlag;
  // add here what you want else
  
} part_t;
   
typedef std::vector<part_t> v_part_t;

typedef std::vector<part_t> vec_part;

typedef  struct part_t_def_1 {
  part_t_def_1(): name1(""), parentTrackID1(0), KinEnergy1(0), ltime1(0), gtime1(0), volFlag1(0), primaryFlag1(0), stepFlag1(0), procFlag1(0), DIF1stepEnergy(0), DIFkillstepEnergy(0), didItDIF(0) {}
  
  G4String name1;
  G4double KinEnergy1;
  G4double ltime1, gtime1;     

  G4int volFlag1, primaryFlag1, procFlag1;

  // determing 1st and last step energies for primaries that DIF only
  G4int parentTrackID1;
  G4int stepFlag1;
  G4double DIF1stepEnergy, DIFkillstepEnergy;
  G4bool didItDIF=false;
  
} part_t_1;

// vector containing particles and information per each particle  
typedef std::vector<part_t_1> v_part_t_1;

class RunAction;
class SNSRootManager;

class EventAction : public G4UserEventAction
{
public:

  EventAction(RunAction*, SNSRootManager*);
  ~EventAction();
  
  void   BeginOfEventAction(const G4Event*);
  void   EndOfEventAction(const G4Event*);
  
  map<int,int> last_secondary_id;
  // Here in this container the key is track ID, the value is (current) number 
  // of secondary particles, which is updated for each step
  
  v_part_t ParticlesBorn;
  void AddParticle(part_t p){ParticlesBorn.push_back(p);}
      
  v_part_t_1 ParticlesBorn1;
  void AddParticle1(part_t_1 p1) {ParticlesBorn1.push_back(p1);}


  v_part_t_1 primary_information_vector;

  void AddPrimary(part_t_1 p_p) {primary_information_vector.push_back(p_p);}

  

  void AddPrimaryTracking(part_t_1 primary, float energy){

  // just fill histos here.
    if (primary.didItDIF == true){

       rootman->FillEnergy_1STEP(primary.primaryFlag1, 
                                 energy);

       rootman->FillEnergy_KILLSTEP(primary.primaryFlag1, 
                                    primary.DIFkillstepEnergy);

    //cout << "E1:" << setw(10) << energy << setw(10) << "Ekill:" << setw(10) << primary.DIFkillstepEnergy << setw(10) << "step:" << setw(10) << primary.stepFlag1 << endl;

    }

  } 

  

  // multidimensional arrays of sets 
  G4bool tracking_logic (G4int id, G4int det, G4int neutrino) {return tracking[det][neutrino].insert(id).second;}

  G4bool energy_tracking_logic (G4int track_id, G4int step, G4int primary) {return tracking_dif_energy[step][primary].insert(track_id).second;}


  G4bool InFlight(G4int id) {return inflight.insert(id).second;}


private:

  RunAction*  ra;
  SNSRootManager* rootman;

  // multidimensional arrays of sets 
  std::set<G4int> tracking [NDet][NNeutrinos];

  std::set<G4int> tracking_dif_energy [NPrimaries][2];

  std::set<G4int> inflight;

  
};

#endif

    
