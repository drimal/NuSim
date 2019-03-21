#ifndef SteppingAction_h
#define SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "SNSRootManager.hh"

#include "G4Proton.hh"
#include "G4AntiProton.hh"
#include "G4PionZero.hh"
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"
#include "G4NeutrinoE.hh"
#include "G4AntiNeutrinoE.hh"
#include "G4NeutrinoMu.hh"
#include "G4AntiNeutrinoMu.hh"
#include "G4NeutrinoTau.hh"
#include "G4AntiNeutrinoTau.hh"
#include "G4Neutron.hh"
#include "G4AntiNeutron.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4Alpha.hh"
#include "G4KaonPlus.hh"
#include "G4KaonMinus.hh"
#include "G4KaonZero.hh"
#include "G4KaonZeroShort.hh"
#include "G4KaonZeroLong.hh"
#include "G4He3.hh"
#include "G4Quarks.hh"
#include "G4Alpha.hh"
#include "G4Lambda.hh"
#include "G4AntiLambda.hh"
#include "G4TauPlus.hh"
#include "G4TauMinus.hh"
#include "G4Deuteron.hh"
#include "G4Upsilon.hh"
#include "G4Triton.hh"
#include "G4Eta.hh"




#include "G4ios.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4TrackStatus.hh"
#include "CLHEP/Units/SystemOfUnits.h"

class RunAction;
class EventAction;
class SNSRootManager;

class SteppingAction : public G4UserSteppingAction
{
public:
  SteppingAction(RunAction*, EventAction*, SNSRootManager* rm);
  ~SteppingAction();
  
  void UserSteppingAction(const G4Step*);
  
  G4String GetTrackStatus(G4TrackStatus stat);
  
private:
  RunAction* ra;
  EventAction* ea;   
  SNSRootManager* rm;
};

#endif
