#include "SteppingAction.hh"
#include "G4Step.hh"
#include "G4PionMinus.hh"
#include "G4PionPlus.hh"
#include "G4PionZero.hh"
#include "G4StepPoint.hh"
#include "G4SteppingManager.hh"
#include "G4TrackStatus.hh"
#include "G4StepStatus.hh"

#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>

using namespace CLHEP;

float energy;
G4String parent_part[10000];
G4int parent_stat[10000];
G4double parent_posz[10000];
G4double parent_energy[10000];

SteppingAction::SteppingAction (RunAction* ract, EventAction* eact, SNSRootManager* rman):ra(ract), ea(eact), rm(rman) { }

SteppingAction::~SteppingAction(){ }

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{

  // get volume of the current step
  G4VPhysicalVolume* volume = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
  G4String vol_name=volume->GetName();

  // get the curent step number
  G4int cur_step = aStep->GetTrack()->GetCurrentStepNumber();

  // get the last step number 
  //G4int last_step = aStep->GetTrack()->Get

  // get track ID of the particle 
  G4ParticleDefinition* part = aStep->GetTrack()->GetDefinition();
  G4int trackID= aStep->GetTrack()->GetTrackID();

  // get the particle name and MCNP ID
  G4String part_name = part -> GetParticleName();
  G4int part_mcnpid = part->GetPDGEncoding();
   
  if (!ea->last_secondary_id[trackID]) // once for each track
    ea->last_secondary_id.insert(make_pair(trackID,0));
  
  // current process
  const G4VProcess *process = aStep->GetPostStepPoint()->GetProcessDefinedStep();
  G4String proc_name = (process) ? process->GetProcessName() : "no_process";

    // list of secondary particles
  const G4TrackVector* secondary = fpSteppingManager->GetSecondary();
  
  // number of new secondaries:
  G4int new_sec = 0;
  new_sec = (*secondary).size()-ea->last_secondary_id[trackID];
 
  // pre step point and position
   G4StepPoint* pre_point = aStep->GetPreStepPoint();
   G4ThreeVector pre_pos = pre_point->GetPosition();

    
  // current point and position
  G4StepPoint* cur_point = aStep->GetPostStepPoint();
  G4ThreeVector cur_pos = cur_point->GetPosition();
    
  // DAR and DIF flags
  G4StepStatus stepStatus = fpSteppingManager->GetStep()->GetPostStepPoint()->GetStepStatus();
  G4String stepstat = (stepStatus==fAtRestDoItProc) ? "AtRest" : "NotAtRest";

  // getting status of the track
  G4TrackStatus trackStat = aStep->GetTrack()->GetTrackStatus();

  // Special counters for muon captures
  G4double mucap_sec [6] = {0};      
  G4String mucap_z="";

  // Special counters tracking pi+Inelastic
  G4int pionPlusInelastic_pi_plus=0;

  // Array flags
  G4int partFlag=1001, procFlag=1001, parentFlag=1001, stepstatFlag=1001, muCapFlag=1001, deathFlag=0, logic=0, newsecFlag=1001, volFlag=1001, detInterest=1001, nuFlag=1001, primaryFlag=10001;

  G4bool endTrackLogic=false;

  // setting a bool to determine the end of life step for a track
  if (trackStat == fStopAndKill){endTrackLogic=true;}

  // setting the volume flag
  if (vol_name == "det_1_physV"){volFlag=0;}
  else if (vol_name == "det_2_physV"){volFlag=1;}
  else if (vol_name == "det_3_physV"){volFlag=2;}
  else {volFlag=3;}


  // setting the new secondaries flag, it's an index to show if any secondaries were produced by a process
  if (new_sec == 0){newsecFlag=0;}
  else if (new_sec == 1){newsecFlag=1;}
  else if (new_sec > 1){newsecFlag=2;}
  else {
    // print warning if newsec less than 0 or greater than 1
    cout << "WARNING: newsecFlag =" << setw(10) << newsecFlag << endl; 
    newsecFlag=10001;
  }

   
  // setting the step status / decay type flag
  if (stepstat == "AtRest"){stepstatFlag=0;}
  else if (stepstat == "NotAtRest"){stepstatFlag=1;}
  else {
    // print warning if stepstat less than 0 or greater than 1
    cout << "WARNING: stepstatFlag =" << setw(10) << stepstat << endl; 
    stepstatFlag=10001;
  }

  // setting the primary flag, i.e., primaries we're interested in tracking
  // this flag is used only for filling the ROOT histos
  if (part == G4Proton::Proton()){primaryFlag=0;} 
  else if (part == G4Neutron::Neutron()){primaryFlag=1;} 
  else if (part == G4PionPlus::PionPlus()){primaryFlag=2;} 
  else if (part == G4PionMinus::PionMinus()){primaryFlag=3;} 
  else if (part == G4MuonPlus::MuonPlus()){primaryFlag=4;} 
  else if (part == G4MuonMinus::MuonMinus()){primaryFlag=5;}
  else if (part == G4NeutrinoMu::NeutrinoMu()){primaryFlag=6;}
  else if (part == G4AntiNeutrinoMu::AntiNeutrinoMu()){primaryFlag=7;}
  else if (part == G4NeutrinoE::NeutrinoE()){primaryFlag=8;}
  else if (part == G4AntiNeutrinoE::AntiNeutrinoE()){primaryFlag=9;}
  else {primaryFlag = 10001;}

  // setting the parent flag
  if (part == G4Proton::Proton()){parentFlag=0;} 
  else if (part == G4AntiProton::AntiProton()){parentFlag=1;} 
  else if (part == G4Neutron::Neutron()){parentFlag=2;} 
  else if (part == G4AntiNeutron::AntiNeutron()){parentFlag=3;} 
  else if (part == G4Gamma::Gamma()){parentFlag=4;} 
  else if (part == G4PionPlus::PionPlus()){parentFlag=5;} 
  else if (part == G4PionMinus::PionMinus()){parentFlag=6;} 
  else if (part == G4PionZero::PionZero()){parentFlag=7;} 
  else if (part == G4MuonPlus::MuonPlus()){parentFlag=8;} 
  else if (part == G4MuonMinus::MuonMinus()){parentFlag=9;} 
  else if (part == G4KaonPlus::KaonPlus()){parentFlag=10;} 
  else if (part == G4KaonMinus::KaonMinus()){parentFlag=11;} 
  else if (part == G4KaonZero::KaonZero()){parentFlag=12;} 
  else if (part == G4KaonZeroShort::KaonZeroShort()){parentFlag=13;} 
  else if (part == G4KaonZeroLong::KaonZeroLong()){parentFlag=14;} 
  else if (part == G4Lambda::Lambda()){parentFlag=15;} 
  else if (part == G4AntiLambda::AntiLambda()){parentFlag=16;} 
  else if (part == G4Electron::Electron()){parentFlag=17;} 
  else if (part == G4Positron::Positron()){parentFlag=18;} 
  else if (part == G4TauPlus::TauPlus()){parentFlag=19;} 
  else if (part == G4TauMinus::TauMinus()){parentFlag=20;} 
  else if (part == G4Alpha::Alpha()){parentFlag=21;} 
  else if (part == G4He3::He3()){parentFlag=22;} 
  else if (part == G4Deuteron::Deuteron()){parentFlag=23;} 
  else if (part == G4Upsilon::Upsilon()){parentFlag=24;} 
  else if (part == G4Eta::Eta()){parentFlag=25;}  
  else if (part == G4Triton::Triton()){parentFlag=26;}  
  else if (part == G4NeutrinoMu::NeutrinoMu()){parentFlag=27; nuFlag=0;}
  else if (part == G4AntiNeutrinoMu::AntiNeutrinoMu()){parentFlag=28; nuFlag=1;}
  else if (part == G4NeutrinoE::NeutrinoE()){parentFlag=29; nuFlag=2;}
  else if (part == G4AntiNeutrinoE::AntiNeutrinoE()){parentFlag=30; nuFlag=3;}
  else if (part_mcnpid > 1000000000){parentFlag=31;}
  else {
    if (parentFlag != 1001){
      // print warning if parent flag is less than 0 or greater than 25
      if (parentFlag < 0 || parentFlag > NParents-1){cout << "PARENT FLAG WARNING:" << setw(10) << "parentFlag =" << setw(10) << parentFlag << "MCNP ID =" << setw(16) << part_mcnpid << endl;}
    }
    parentFlag=32;
    //output cout if parent flag == 1001
  }


  // Setting the process flag
  // this list comprises ALL possible physics processes as of Geant4 10.0
  if (proc_name == "Transportation"){procFlag=0;}
  else if (proc_name == "hBrems"){procFlag=1;}
  else if (proc_name == "annihil"){procFlag=2;}
  else if (proc_name == "conv"){procFlag=3;}
  else if (proc_name == "photonNuclear"){procFlag=4;}
  else if (proc_name == "hadElastic"){procFlag=5;}
  else if (proc_name == "msc"){procFlag=6;}
  else if (proc_name == "hPairProd"){procFlag=7;}
  else if (proc_name == "CoulombScat"){procFlag=8;}
  else if (proc_name == "muIoni"){procFlag=9;}
  else if (proc_name == "electronNuclear"){procFlag=10;}
  else if (proc_name == "hIoni"){procFlag=11;}
  else if (proc_name == "eIoni"){procFlag=12;}
  else if (proc_name == "phot"){procFlag=13;}
  else if (proc_name == "muBrems"){procFlag=14;}
  else if (proc_name == "positronNuclear"){procFlag=15;}
  else if (proc_name == "ionIoni"){procFlag=16;}
  else if (proc_name == "eBrem"){procFlag=17;}
  else if (proc_name == "compt"){procFlag=18;}
  else if (proc_name == "muPairProd"){procFlag=19;}
  else if (proc_name == "pi+Inelastic"){procFlag=20;}
  else if (proc_name == "kaon0LInelastic"){procFlag=21;}
  else if (proc_name == "sigma-Inelastic"){procFlag=22;}
  else if (proc_name == "xi-Inelastic"){procFlag=23;}
  else if (proc_name == "omega-Inelastic"){procFlag=24;}
  else if (proc_name == "anti_deuteronInelastic"){procFlag=25;}
  else if (proc_name == "hFritiofCaptureAtRest"){procFlag=26;}
  else if (proc_name == "tInelastic"){procFlag=27;}
  else if (proc_name == "nKiller"){procFlag=28;}
  else if (proc_name == "neutronInelastic"){procFlag=29;}
  else if (proc_name == "pi-Inelastic"){procFlag=30;}
  else if (proc_name == "kaon0SInelastic"){procFlag=31;}
  else if (proc_name == "anti_sigma-Inelastic"){procFlag=32;}
  else if (proc_name == "anti_xi-Inelastic"){procFlag=33;}
  else if (proc_name == "anti_omega-Inelastic"){procFlag=34;}
  else if (proc_name == "anti_tritonInelastic"){procFlag=35;}
  else if (proc_name == "hBertiniCaptureAtRest"){procFlag=36;}
  else if (proc_name == "He3Inelastic"){procFlag=37;}
  else if (proc_name == "nCapture"){procFlag=38;}
  else if (proc_name == "kaon+Inelastic"){procFlag=39;}
  else if (proc_name == "lambdaInelastic"){procFlag=40;}
  else if (proc_name == "sigma+Inelastic"){procFlag=41;}
  else if (proc_name == "xi0Inelastic"){procFlag=42;}
  else if (proc_name == "anti_protonInelastic"){procFlag=43;}
  else if (proc_name == "anti_He3Inelastic"){procFlag=44;}
  else if (proc_name == "muMinusCaptureAtRest"){procFlag=45;}
  else if (proc_name == "alphaInelastic"){procFlag=46;}
  else if (proc_name == "protonInelastic"){procFlag=47;}
  else if (proc_name == "kaon-Inelastic"){procFlag=48;}
  else if (proc_name == "anti-lambdaInelastic"){procFlag=49;}
  else if (proc_name == "anti_sigma+Inelastic"){procFlag=50;}
  else if (proc_name == "anti_xi0Inelastic"){procFlag=51;}
  else if (proc_name == "anti_neutronInelastic"){procFlag=52;}
  else if (proc_name == "anti_alphaInelastic"){procFlag=53;}
  else if (proc_name == "dInelastic"){procFlag=54;}
  else if (proc_name == "ionInelastic"){procFlag=55;}
  else if (proc_name == "Decay"){procFlag=56;}
  else {
    if (procFlag != 1001){
      // output statement if process flag actually is greater than 56 or less than 0
      if (procFlag < 0 || procFlag > NProc-1){cout << "PROCESS FLAG WARNING:" << setw(10) << "procFlag =" << setw(10) << procFlag << "Process Name =" << setw(16) << proc_name << endl;}
    }
    procFlag=10001;
      
  }

    parent_part[trackID] = part_name;
    int parID = aStep->GetTrack()->GetParentID();

    if (procFlag == 30){
        cout << aStep->GetTrack()->GetCreatorProcess()->GetProcessName() << " " << parent_part[parID] << endl;
    }



/*
    parent_part[trackID] = part_name;
    parent_stat[trackID] = stepstatFlag;
    parent_posz[trackID] = cur_pos.z()/cm;
    parent_energy[trackID] = aStep->GetTrack()->GetKineticEnergy()/MeV;

    int parID = aStep->GetTrack()->GetParentID();

   // if(pre_pos.z()/cm > 20 && primaryFlag > 1 && primaryFlag < 6 && proc_name == "Decay"  && stepstat == "NotAtRest"){
        if(procFlag == 56 && endTrackLogic == true && primaryFlag < NPrimaries){
            //rm->Fill1D(pre_pos.z()/cm);

            if(cur_pos.z()/cm > 20){
                G4cout << vol_name << setw(10) << part_name << setw(10) << cur_step << setw(10) << pre_pos.z()/cm << setw(10) << stepstat << setw(25) << aStep->GetTrack()->GetCreatorProcess()->GetProcessName() << setw(10) << parent_part[parID] <<" "<< parID << " " << parent_posz[parID] <<" " << parent_energy[parID] << endl;
            }

         // rm->Fill3D(pre_pos.x()/cm, pre_pos.y()/cm, pre_pos.z()/cm);
        //rm->Fill2D(pre_pos.x()/cm, pre_pos.y()/cm);
        //rm->Fill1D(pre_pos.z()/cm);
        //if(pre_pos.z()/cm > 20){
        //}
        
    }

    //int parID = aStep->GetTrack()->GetParentID();

    if(nuFlag != 1001 && cur_step == 1){

        if(pre_pos.z()/cm > 20){
            //G4cout << cur_step <<" "<< part_name << " " << parent_part[parID] << " " << parent_stat[parID] <<" " <<aStep->GetTrack()->GetCreatorProcess()->GetProcessName() <<" "<< pre_pos.z()/cm << endl;
            rm->Fill1D(parent_stat[parID]);
            rm->Fill2D(cur_pos.x()/cm, cur_pos.y()/cm);
        }

    }
*/
  // before scanning through the secondaries


  // currently only concerned with tracking neutrinos in the detector volumes
  // converting to minimized primary flag system SOON...
  if (nuFlag<4){

    part_t_1 primary;

    primary.name1 = aStep->GetTrack()->GetDefinition()->GetParticleName();
    primary.parentTrackID1 = aStep->GetTrack()->GetParentID();
    primary.KinEnergy1 = aStep->GetTrack()->GetKineticEnergy()/MeV;
    primary.volFlag1 = volFlag;
    primary.primaryFlag1 = parentFlag;
    primary.procFlag1 = procFlag;
    primary.stepFlag1 = cur_step;
    primary.ltime1 = aStep->GetTrack()->GetLocalTime()/ns;
    primary.gtime1 = aStep->GetTrack()->GetGlobalTime()/ns;

    ea->AddParticle1(primary);

    if (ea->tracking_logic(trackID, volFlag, nuFlag)){

      // 1D histo for tracking energy of neutrinos that made it to a detector
      rm->FillEnergy(nuFlag, 
                     volFlag, 
                     aStep->GetTrack()->GetKineticEnergy()/MeV, 
                     0);

      // 1D histo for tracking the time it took for neutrinos to make it to a detector (not relative to p+Hg currently)
      rm->FillTime(nuFlag, 
                   volFlag, 
                   aStep->GetTrack()->GetLocalTime()/ns, 
                   0);

    }
  }

// staying within some bounds so I don't record every single particle that rolls through here

if (primaryFlag < NPrimaries){

   part_t_1 primary_;


   primary_.primaryFlag1 = primaryFlag;
   primary_.parentTrackID1 = trackID;
   primary_.stepFlag1 = cur_step;



   if (cur_step == 1){energy = aStep->GetTrack()->GetKineticEnergy()/MeV;}

   if (endTrackLogic == true){primary_.DIFkillstepEnergy = aStep->GetTrack()->GetKineticEnergy()/MeV;}

   if ((procFlag == 56) && (stepstatFlag == 1)){primary_.didItDIF = true; 
                                                primary_.didItDIF = true;}


   ea->AddPrimaryTracking(primary_, energy);

}

  if (parentFlag < NParents){
  
    // counting all the processes a particle engaged in, also counts the processes that kill a primary
    ra->CountProcesses(parentFlag, procFlag);
  
    // added discrimination for determining number of secondaries produced by a parent interaction
    ra->CountProcessesTracking(parentFlag, procFlag, newsecFlag);
  }


  // counting the number of unsuccessful interactions, used for failed inelastic collisions
  if (parentFlag < NParents && parentFlag >= 0){
    if (procFlag < NProc){
      if (new_sec == 0){
	      ra->CountFailedProcesses(parentFlag, procFlag);
      }
    }
  }


  // filling the death counter histos, using the destructive process flags
  // staying within parimay flag bounds, and proccess flag bounds
  if (primaryFlag < NPrimaries && primaryFlag >= 0){
    if (procFlag >= 20 && procFlag < NProc){
      rm->FillDestroyed(primaryFlag, procFlag);
    }
  }

/*
  // trying to detect an unsucessful pi+Inelastic event
  if (procFlag == 20){
    if ( new_sec <= 1){
      cout << "*****************CRITICAL WARNING*****************: pi+Inelastic interaction was unsucessfull" << setw(30) << "# of secondaries =" << setw(20) << new_sec << endl;
    }
  }

  // trying to determine times for p+Hg interaction

  // tracking the SOURCE protons from the particle beam
  if (trackID == 1){
  cout << setw(10) << "particle:" << setw(10) << part_name << setw(10) << "TrackID:" << setw(10) << trackID << setw(10) << "Step:" << setw(6) << cur_step << setw(12) << "Process:" << setw(20) << proc_name << setw(20)
  << setw(15) << "volume:" << setw(30) << vol_name << setw(10) << "LT:" << setw(15) << aStep->GetTrack()->GetLocalTime()/ns << setw(10) << "GT:" << setw(15) << aStep->GetTrack()->GetGlobalTime()/ns << setw(10) 
  << "E:" << setw(10) << aStep->GetTrack()->GetKineticEnergy()/MeV << endl;
  }


  // looking at the output for mu- capture time in Hg target
  if (procFlag == 45){
  cout << setw(10) << "particle:" << setw(10) << part_name << setw(10) << "TrackID:" << setw(10) << trackID << setw(10) << "Step:" << setw(6) << cur_step << setw(12) << "Process:" << setw(22) << proc_name << setw(20)
  << setw(15) << "volume:" << setw(30) << vol_name << setw(10) << "LT:" << setw(15) << aStep->GetTrack()->GetLocalTime()/ns << setw(10) << "GT:" << setw(15) << aStep->GetTrack()->GetGlobalTime()/ns << setw(10) 
  << "PT:" << setw(10) << aStep->GetTrack()->GetProperTime()/ns << setw(10) << "E:" << setw(10) << aStep->GetTrack()->GetKineticEnergy()/MeV << endl;
  }
  */


  // now begin scanning through the secondary list
  for (size_t lp = ea->last_secondary_id[trackID]; lp<(*secondary).size(); lp++) {
    G4Track *secondary_track = (*secondary)[lp];
    G4ParticleDefinition* sec_part = secondary_track->GetDefinition();
    G4String sec_name=sec_part->GetParticleName();
    G4int sec_trackID= secondary_track->GetTrackID();

    // setting the secondary particle type flag
    if (sec_name == "pi+"){partFlag=0;}
    else if (sec_name == "pi-"){partFlag=1;}
    else if (sec_name == "mu+"){partFlag=2;}
    else if (sec_name == "mu-"){partFlag=3;}
    else if (sec_name == "nu_mu"){partFlag=4;}
    else if (sec_name == "anti_nu_mu"){partFlag=5;}
    else if (sec_name == "nu_e"){partFlag=6;}
    else if (sec_name == "anti_nu_e"){partFlag=7;}
    else if (sec_name == "neutron"){partFlag=8;} 
    else if (sec_name == "proton"){partFlag=9;}
    else if (sec_name == "e-"){partFlag=10;}
    else if (sec_name == "gamma"){partFlag=11;}
    else if (sec_name == "e+"){partFlag=12;}
    else {
      if (partFlag != 1001 && partFlag != 10001){
	// output statement if particle flag actually is greater than 12 or less than 0
	if (partFlag < 0 || partFlag > NPart-1){cout << "SECONDARY PARTICLE FLAG WARNING:" << setw(10) << "partFlag =" << setw(10) << partFlag << "Secondary Name =" << setw(16) << sec_name << endl;}
      }      partFlag=10001;
    }

    // setting the properties of the secondary particle using Bari's structure
    // right now discrimating for mu- and target volume
    if (partFlag == 3 && vol_name == "physical_vol_target"){

      part_t particle;

      particle.name=aStep->GetTrack()->GetDefinition()->GetParticleName();
      particle.sec_name=secondary_track->GetDefinition()->GetParticleName();
      particle.proc=process->GetProcessName();
      particle.id=sec_part->GetPDGEncoding();
      particle.ltime=aStep->GetTrack()->GetLocalTime()/ns;
      particle.gtime=aStep->GetTrack()->GetGlobalTime()/ns;
      particle.KinEnergy=aStep->GetTrack()->GetKineticEnergy()/MeV;
      particle.trackID=aStep->GetTrack()->GetTrackID();
      particle.volFlag=volFlag;
      particle.partFlag=partFlag;
      particle.stepFlag = cur_step;
      particle.procFlag = procFlag;
	  
      ea->AddParticle(particle);
    }
     
    // making sure I don't go out of array bounds for filling histo arrays
    if ((partFlag < NPart) && (procFlag < NProc)){

      // processes that generated a secondary particle
      rm->FillProcess(partFlag, procFlag);

      // MCNP ID of each secondary generated
      rm->FillMCNPID(partFlag, part->GetPDGEncoding());
   
      // initial/emission energy of secondaries generated
      rm->FillEnergy(partFlag, secondary_track->GetKineticEnergy()/MeV);

    }


    // making sure I don't go out of my array bounds
    if (partFlag < NPart){

      // counting all the processes that produced specific secondaries
      ra->CountProcessProducedParticle(partFlag, procFlag);

      // generically counting all secondary particle types created
      ra->CountAll(partFlag);

      // counting secondary particles from only the initial p+Hg interaction
      if (parentFlag == 0 && trackID == 1 && vol_name == "physical_vol_target"){
	      ra->CountProtonPlusHgProduction(partFlag);
      }

      // recording the emission/creation energy and time relative to p+Hg in a 2D histo
      rm->Fill2DEnergyTime(partFlag, 
		                       secondary_track->GetKineticEnergy()/MeV, 
		                       secondary_track->GetGlobalTime()/ns, 
		                       1);

      // recording the time of emission relative to p+Hg (initTime=0) in a 1D histo
        //if(partFlag > 3 && partFlag < 8 && secondary_track->GetGlobalTime()/ns < 1){
            //cout << sec_name << " " << secondary_track->GetGlobalTime()/ns << endl;
        //}
      rm->FillTime(partFlag, 
                   secondary_track->GetGlobalTime()/ns,
                   1);

      // recording the initial energy of emitted secondary in a 1D histo
      rm->FillEnergy(partFlag,
                     secondary_track->GetKineticEnergy()/MeV);

      if (parentFlag < NParents){

	// counting the specific decays that produced a particle
	if ((procFlag == 56) && (stepstatFlag <= 1)){
	  ra->CountDecayProducedParticle(stepstatFlag, partFlag,  parentFlag);
	}

	// counting all the parents of specific particles
	ra->CountParentParticles(partFlag, parentFlag);

      }
    }

/*
    // some specific checking of neutrons and pi+ 	
    if (sec_name == "neutron" ||  partFlag == 0){

      // checking to see if there are any parent particles I have not included
      if (parentFlag > NParents-1){
	cout << "SECONDARY PARENT WARNING:" << setw(20) << "Parent= " << setw(10) << part_name << setw(20) << "Secondary=" << setw(10) << sec_name << setw(20) << "MCNP Parent ID =" << setw(15) 
	     << part->GetPDGEncoding() << setw(30) << "Process =" << setw(20) << proc_name << endl;
      } 

      // checking to see if there are any processes that may be hiding...
      if (procFlag > NProc){cout << "SECONDARY PRODUCTION PROCESS WARNING:" << setw(20) << "Secondary Name =" << setw(10) << sec_name << setw(20) << "Produced by Process ="  << setw(15) << proc_name << endl;}
      
    }
*/

    // counting secondaries from mu- capture
    if (proc_name == "muMinusCaptureAtRest"){ 
      if (partFlag == 10){mucap_sec[0]++;}
      if (partFlag == 7){mucap_sec[1]++;}
      if (partFlag == 6){mucap_sec[2]++;}
      if (partFlag == 4){mucap_sec[3]++;}
      if (partFlag == 11){mucap_sec[4]++;}
      if (partFlag == 8){mucap_sec[5]++;}
      mucap_z = (mucap_z + sec_name + " ");
    }


    // looking specifically into absolute deaths, i.e., the parent did not produce a secondary of it's same type (for now just looking into pi+Inelastic)
    // also looking into the number of pi+'s produced per each pi+Inelastic collision
    // set a logic variable, 1 means (pi+ -> pi+) and 0 means (pi+ -> not pi+)
    if (procFlag == 20){
      if (partFlag == 0){
	      logic=1;
	      pionPlusInelastic_pi_plus++;
      } 
    }


    //if (procFlag == 38 && partFlag == 8){cout << "nCapture Event Produced neutron" << setw(20) << sec_name << endl; }


  }// end of secondary list cycle





  // passing the number of pi+ secondaries per each pi+Inelastic from inside the loop into an array
  if ((procFlag == 20) && (logic == 1)){ra->CountTrackData(pionPlusInelastic_pi_plus);}


  // counting the number of absolute deaths, i.e., the parent did not produce a secondary of it's same type
  // here only looking into pi+Inelastic
  if (procFlag == 20){
    if (logic == 0){
      if(new_sec > 0){     
	      ra->CountDestroyedParticle(deathFlag); 
	      //cout << ra-> DestroyedParticle[0] << setw(20) << proc_name << setw(20) << "logic =" << setw(10) << logic << endl;
      }
    }
  }


  // making sure it's the end of track, last step
  // making sure I don't go out of array bounds for filling histo arrays
  if ((endTrackLogic==true) && (primaryFlag < NPrimaries) && (procFlag < NProc)){

    // recording the lifetime of a primary that DAR or DIF, in a 1D histo
    // making sure I'm within DAR/DIF flag bounds
    // making sure I'm only counting decays
    if ((stepstatFlag < NDecayTypes) && (procFlag==56)){
        rm->FillTime(primaryFlag,
                     stepstatFlag,
                     aStep->GetTrack()->GetLocalTime()/ns);


    // recording the time when a primary DAR/DIF relative to p+Hg
    rm->FillTime(primaryFlag, 
                 stepstatFlag, 
                 0, 
                 aStep->GetTrack()->GetGlobalTime()/ns);


    // recording the energy of primaries at DAR/DIF
    rm->FillEnergy(primaryFlag, 
                   stepstatFlag, 
                   0, 
                   aStep->GetTrack()->GetKineticEnergy()/MeV);

    }

    // recording the lifetime of a primary in a 1D histo
    rm->FillTime(primaryFlag, 
                 aStep->GetTrack()->GetLocalTime()/ns,
                 0);

  // making sure sure the process occurs "Not At Rest", i.e., In Flight
  // making sure the process is actually Decay
    if((stepstatFlag==1) && (procFlag==56)){

      rm->Fill2DEnergyTime(primaryFlag, 
		                       aStep->GetTrack()->GetKineticEnergy()/MeV, 
		                       aStep->GetTrack()->GetLocalTime()/ns, 
		                       0);
    }

  }


  // passing the mu- capture array information from inside the secondaries loop into a storage array
  if (proc_name == "muMinusCaptureAtRest"){
    ra->CountmuProcessProducts(mucap_sec, mucap_z);
  }

	
  // update the number of secondary particles for the current track ID
  ea->last_secondary_id[trackID] = (*secondary).size();
}//end of stepping action



