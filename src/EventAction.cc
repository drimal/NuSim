#include "EventAction.hh"
#include "RunAction.hh"
#include "SNSRootManager.hh"
#include "G4Event.hh"
#include "Randomize.hh"

EventAction::EventAction(RunAction* run, SNSRootManager* Rman): ra(run), rootman(Rman){}

EventAction::~EventAction(){}


void EventAction::BeginOfEventAction(const G4Event* evt)
{

  // outputing the current event ID number
  G4int evtID = evt->GetEventID();
    if(evtID % 500000 == 0) {
        cout << "processed:"<< evtID << " events" << endl;
    }
  //G4cout << setw(5) << "Event" << setw(10) << evtID << G4endl;

  // clearing the tracking information for neutrinos
  for (int det=0; det<NDet; det++){
    for (int nu=0; nu<NNeutrinos; nu++){
      tracking [det][nu].clear();
    }
  }

  // clearing the tracking information for primaries
  for (int step=0; step<2; step++){
    for (int p=0; p<NPrimaries; p++){
      tracking_dif_energy [step][p].clear();
    }
  }

  inflight.clear();
  
  ParticlesBorn.clear();
  ParticlesBorn1.clear();

  primary_information_vector.clear();
  
  last_secondary_id.clear();
}

void EventAction::EndOfEventAction(const G4Event* evt)
{
  


  for (unsigned i=0; i<primary_information_vector.size(); i++) {
    part_t_1  par = primary_information_vector[i];

    // confirming that it was a DIF
    if (par.didItDIF == true){


    }

  }









							   
}// end of event action
