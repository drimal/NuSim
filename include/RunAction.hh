#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"
#include "SNSRootManager.hh"
#include "G4String.hh"
#include <string>


// Presetting the array dimensions
const int NPart=13;
const int NProc=57;
const int NParents=33;
const int NTracks=100000;
const int NmuMinusCaptures=70000;
const int NmuMinusCaptureSec=6;
const int NDecayTypes=2;
const int NDestructProc=26;
const int NLogicalNewSec=2;
const int NDet=4;
const int NNeutrinos=4;
const int NPrimaries=10;

class G4Run;

class RunAction : public G4UserRunAction
{
public:

  RunAction(SNSRootManager*, int, int);
  ~RunAction();
  
  void BeginOfRunAction(const G4Run*);
  void EndOfRunAction(const G4Run*);
  
  

  G4int ProtonsOnTarget;
  G4int EventsProcessed;
  G4int runNumber;

  // a string array containing all the process names, used for the output file
  G4String processNames[NProc] = {
    "Transportation",
    "hBrems",
    "annihil",
    "conv",
    "photonNuclear",
    "hadElastic",
    "msc",
    "hPairProd",
    "CoulombScat",
    "muIoni",
    "electronNuclear",
    "hIoni",
    "eIoni",
    "phot",
    "muBrems",
    "positronNuclear",
    "ionIoni",
    "eBrem",
    "compt",
    "muPairProd",
    "pi+Inelastic",
    "kaon0LInelastic",
    "sigma-Inelastic",
    "xi-Inelastic",
    "omega-Inelastic",
    "anti_deuteronInelastic",
    "hFritiofCaptureAtRest",
    "tInelastic",
    "nKiller",
    "neutronInelastic",
    "pi-Inelastic",
    "kaon0SInelastic",
    "anti_sigma-Inelastic",
    "anti_xi-Inelastic",
    "anti_omega-Inelastic",
    "anti_tritonInelastic",
    "hBertiniCaptureAtRest",
    "He3Inelastic",
    "nCapture",
    "kaon+Inelastic",
    "lambdaInelastic",
    "sigma+Inelastic",
    "xi0Inelastic",
    "anti_protonInelastic",
    "anti_He3Inelastic",
    "muMinusCaptureAtRest",
    "alphaInelastic",
    "protonInelastic",
    "kaon-Inelastic",
    "anti-lambdaInelastic",
    "anti_sigma+Inelastic",
    "anti_xi0Inelastic",
    "anti_neutronInelastic",
    "anti_alphaInelastic",
    "dInelastic",
    "ionInelastic",
    "Decay",
  };


  // a string array containing all the particle names that we are interested in, used for the output file
  G4String partNames[NPart] = {
    "pi+",
    "pi-",
    "mu+",
    "mu-",
    "nu_mu",
    "anti_nu_mu",
    "nu_e",
    "anti_nu_e",
    "neutron",
    "proton",
    "e-",
    "gamma",
    "e+"
  };

  // a string array containing all the parent particle names, used for the output file
  G4String parentNames[NParents] = {
    "Proton",
    "Anti-Proton",
    "Neutron",
    "Anti-Neutron",
    "Gamma",
    "Pion+",
    "Pion-",
    "Pion0",
    "Muon+",
    "Muon-",
    "Kaon+",
    "Kaon-",
    "Kaon0",
    "Kaon0Short", 
    "Kaon0Long",
    "Lambda",
    "Anti-Lambda", 
    "e-",
    "e+",
    "Tau+",
    "Tau-",
    "Alpha",
    "He3",
    "Deuteron",
    "Upsilon",
    "Eta",
    "Triton",
    "NuMu",
    "Anti-NuMu",
    "NuE",
    "Anti-NuE",
    "Nuclear",
    "Unknown"
  };

  // a string array containing the type of capture process mu- engaged in
  G4String muCapProcessNames[4] = {
    "Ordinary muon capture",
    "Radiative muon capture",
    "K-shell capture",
    "Normal decay mode"
  };

  // a string array containing the possible product names of mu- capture
  G4String muCapProductNames[8] = {
    "nu_mu",
    "nu_e",
    "anti_nu_e",
    "gamma",
    "e-",
    "neutron",
    "Z"
  };

  // a string array containing the detector names
  G4String volNames[NDet] = {
    "detector_oscsns",
    "detector_near",
    "detector_forward_15",
    "detector_other_space"
};

  // a string array containing the decay type names
  G4String decayNames[NDecayTypes+1] = {
    "DAR",
    "DIF",
    "all"
};

  // all counter
  G4double All [NPart]={0};
  void CountAll(int partFlag);

  // counting all the processes a particle engaged in
  G4double Processes [NParents][NProc]={{0}};
  void CountProcesses(int parentFlag, int procFlag);

  // counting all the processes a particle engaged in, with added discrimination for number of new secondaries
  G4double ProcessesTracking [NParents][NProc][NLogicalNewSec]={{{0}}};
  void CountProcessesTracking(int parentFlag, int procFlag, int newsecFlag);

  // counting all the parents of specific particles
  G4double ParentParticles [NPart][NParents] = {{0}};
  void CountParentParticles(int partFlag, int parentFlag);

  // counting all the processes that produced specific particles
  G4double ProcessProducedParticle [NPart][NProc] = {{0}};
  void CountProcessProducedParticle(int partFlag, int procFlag);

  // counting the specific decays that produced a particle
  G4double DecayProducedParticle [NDecayTypes][NPart][NParents] = {{{0}}};
  void CountDecayProducedParticle(int stepstatFlag, int partFlag, int parentFlag);

  // counting the absolute particle deaths
  G4double DestroyedParticle [NDestructProc] = {0};
  void CountDestroyedParticle(int deathFlag);

  // counting the number of pi+'s produced per each pi+Inelastic process
  G4int track=0;

  G4double TrackData [NTracks] = {0};
  void CountTrackData(int pionPlusInelastic_pi_plus);

  // counting the number of times pi+Inelastic was called, inside the secondary loop
  G4double secLoopProcesses [NPart][NProc] = {{0}};
  void CountsecLoopProcesses(int parentFlag, int procFlag);

  // counting mu minus capture information
  G4int muCapInstance=0;
  G4double mucap_secondaries [NmuMinusCaptures][NmuMinusCaptureSec] = {{0}};
  G4String muCapZ [NmuMinusCaptures] = {""};
  void CountmuProcessProducts(double mucap_sec[], string mucap_z);

  // counting the number of unsuccessful interaction attempts
  G4double FailedProcesses [NParents][NProc] = {{0}};
  void CountFailedProcesses(int parentFlag, int procFlag);

  // counting the number of particles strictly from p+Hg interaction, primary protons ONLY
  G4double ProtonPlusHgProduction [NPart];
  void CountProtonPlusHgProduction(int partFlag);

private:
  SNSRootManager* rm;
};

#endif


