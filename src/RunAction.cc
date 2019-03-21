#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include <iomanip>
#include <iostream>
#include <fstream>

RunAction::RunAction(SNSRootManager* rman, int run, int NEvents) : rm(rman) {
    runNumber = run;
    EventsProcessed = NEvents;
    ProtonsOnTarget = NEvents;
}

RunAction::~RunAction(){}


// all counter
void RunAction::CountAll(int partFlag){ All [partFlag] ++;}

// counting all the processes a particle engaged in
void RunAction::CountProcesses(int parentFlag, int procFlag){ Processes [parentFlag][procFlag] ++;}

// counting all the processes a particle engaged in, with added discrimination for number of new secondaries
void RunAction::CountProcessesTracking(int parentFlag, int procFlag, int newsecFlag){ ProcessesTracking [parentFlag][procFlag][newsecFlag] ++;}

// counting the specific decays that produced a particle
void RunAction::CountDecayProducedParticle(int stepstatFlag, int partFlag, int parentFlag){ DecayProducedParticle [stepstatFlag][partFlag][parentFlag] ++;}

// counting all the parents of specific particles
void RunAction::CountParentParticles(int partFlag, int parentFlag){ParentParticles [partFlag][parentFlag] ++;}

// counting all the processes that produced specific particles
void RunAction::CountProcessProducedParticle(int partFlag, int procFlag){ProcessProducedParticle [partFlag][procFlag] ++;}

// counting the pi+Inelastic collisions that resulted in a death and no secondary pi+'s
void RunAction::CountDestroyedParticle(int deathFlag){DestroyedParticle [deathFlag] ++;}

// counting the number of pi+'s produced per each pi+Inelastic process
void RunAction::CountTrackData(int pionPlusInelastic_pi_plus){
  TrackData [track] = pionPlusInelastic_pi_plus;
  track++;

}

// counting the number of times pi+Inelastic was called, inside the secondary loop
void RunAction::CountsecLoopProcesses(int parentFlag, int procFlag){ secLoopProcesses [parentFlag] [procFlag] ++;}

// counting mu minus capture information
void RunAction::CountmuProcessProducts(double mucap_sec[], string mucap_z){
  for (int d=0; d<6; d++){
    mucap_secondaries[muCapInstance][d] = mucap_sec[d];
  }
  muCapZ [muCapInstance] = mucap_z;
  muCapInstance++;
}

// counting the number of unsuccessful interaction attempts
void RunAction::CountFailedProcesses(int parentFlag, int procFlag){FailedProcesses [parentFlag] [procFlag] ++;}

// counting the number of particles strictly from p+Hg interaction, primary protons ONLY
void RunAction::CountProtonPlusHgProduction(int partFlag){ProtonPlusHgProduction [partFlag] ++;}




void RunAction::BeginOfRunAction(const G4Run* aRun)
{

  //G4cout<<"Enter the number of incident protons:"<<G4endl;
  //G4cin>>ProtonsOnTarget;

  G4cout<< " " <<G4endl;

  //G4cout<<"Enter the number of events being processed:"<<G4endl;
  //G4cout<<"(# you entered after the /run/beamOn command)"<<G4endl;
  //G4cin>>EventsProcessed;

  // setting the number of protons on target equal to the number I'm actually processing
  //ProtonsOnTarget = EventsProcessed;

  rm->Init();

}



void RunAction::EndOfRunAction(const G4Run* aRun)
{

    cout << "run Number = " << runNumber << endl;
   // saves the histogram information to a ROOT file
    rm->Save();

  // saves the counter information to a formatted txt file
    ofstream myfile;
    G4String fileName = Form("counters_%03d.txt", runNumber);


    ofstream myfile_;
    G4String fileName_ = Form("data_%03d.txt", runNumber);


    ofstream myfile_1;
    G4String fileName_1 = Form("dpiInelasticprocess_%03d.txt", runNumber);
    
    myfile.open(fileName);
    myfile_.open(fileName_);
    myfile_1.open(fileName_1);

  myfile <<"_______________________________________________________________________________________________________________________________________________"<< G4endl;
  myfile <<" "<< G4endl;
  myfile <<"       /////////              //////////      ///              //          //    //////////////       //  //      "<< G4endl;
  myfile <<"     //        //            //             //  //            ////        //          //            //   //       "<< G4endl;
  myfile <<"    //          //          //             //    //          // //       //          //           //    //        "<< G4endl;
  myfile <<"   //                      //             //     //         //   //     //          //          //     //         "<< G4endl;
  myfile <<"   //                     //////////     //////////        //    //    //          //          ////////////       "<< G4endl;
  myfile <<"   //        ////////    //             //       //       //      //  //          //              //              "<< G4endl;
  myfile <<"    //         //       //             //        //      //       // //          //              //               "<< G4endl;
  myfile <<"     //       //       //             //         //     //        ////          //              //                "<< G4endl;
  myfile <<"       ///////        //////////     //          //    //          //          //              //                 "<< G4endl;
  myfile <<"_______________________________________________________________________________________________________________________________________________"<< G4endl;
  myfile <<""<< G4endl;
  myfile <<"                                           SIMULATION OUTPUT FILE                                                 "<< G4endl;
  myfile <<"_______________________________________________________________________________________________________________________________________________"<< G4endl;
  myfile <<""<< G4endl;
  myfile <<""<< G4endl;

  myfile << "Simulation parameters :" << G4endl;
  myfile <<""<< G4endl;
  myfile << setw(45) << "Total number of incident protons = "<< setw(8) << ProtonsOnTarget <<G4endl;
  myfile << setw(45) << "Total number of events = "<< setw(8) << EventsProcessed <<G4endl;
  myfile <<""<< G4endl;
  myfile <<""<< G4endl;


  myfile <<"_______________________________________________________________________________________________________________________________________________ "<< G4endl;
  myfile <<""<< G4endl;
  myfile <<"                                         OUTSIDE THE SECONDARY LOOP                                 "<< G4endl;
  myfile <<"_______________________________________________________________________________________________________________________________________________ "<< G4endl;


  myfile <<""<< G4endl;
  myfile <<""<< G4endl;
  myfile <<"                                  Processes a particle engaged in                                                 "<< G4endl;
  myfile <<""<< G4endl;
  myfile <<""<< G4endl;


  myfile << setw(40) << "" << setw(20) << parentNames[5] << setw(16) << parentNames[6] << setw(16) << parentNames[2] << setw(16) 
	 << parentNames[8] << setw(16) << parentNames[9] << setw(16) << parentNames[10] << G4endl;

  for (int j=0; j<NProc; j++){
    myfile << setw(40) << processNames[j] << setw(20) << Processes[5][j] << setw(16) << Processes[6][j] << setw(16) << Processes[2][j] << setw(16)
	   << Processes[8][j] << setw(16) << Processes[9][j] << setw(16) << Processes[10][j] << G4endl;
  }

  myfile <<""<< G4endl;
  myfile <<"                                  Processes engaged sum                               "<< G4endl;
  myfile <<""<< G4endl;

  // determing the sum of processes a particle engaged in
  G4double process_engaged_sum [NParents] = {0};

  for (int p=0; p<=NPart; p++){
    for (int proc=0; proc<=NProc; proc++){
      process_engaged_sum[p] = process_engaged_sum[p] + Processes[p][proc];
    }
  }


  myfile << setw(40) << "" << setw(20) << process_engaged_sum[5] << setw(16) << process_engaged_sum[6]  << setw(16) << process_engaged_sum[2] << setw(16) << process_engaged_sum[8] << setw(16) 
	 << process_engaged_sum[9] << setw(16) << process_engaged_sum[10] << setw(16) << G4endl;



  myfile <<""<< G4endl;
  myfile <<""<< G4endl;

  myfile <<"_______________________________________________________________________________________________________________________________________________"<< G4endl;
  myfile <<""<< G4endl;
  myfile <<"                                       INSIDE THE SECONDARY LOOP                                   "<< G4endl;
  myfile <<"_______________________________________________________________________________________________________________________________________________"<< G4endl;

  myfile <<""<< G4endl;
  myfile <<""<< G4endl;

  myfile <<""<< G4endl;
  myfile <<""<< G4endl;
  myfile <<"                                  All counter results for particular particles                             "<< G4endl;
  myfile <<""<< G4endl;
  myfile <<""<< G4endl;

  myfile << setw(40) << "" << setw(20) << partNames[0] << setw(16) << partNames[1] << setw(16) << partNames[2] << setw(16) 
	 << partNames[3]  << setw(16) << partNames[4] << setw(16) << partNames[5] << setw(16)
	 << partNames[6]  << setw(16) << partNames[7] << setw(16) << partNames[10] << setw(16)
	 << partNames[12] << setw(16) << partNames[8] << G4endl;

  myfile << setw(40) << "" << setw(20) << All[0] << setw(16) << All[1] << setw(16) << All[2] << setw(16) 
	 << All[3]  << setw(16) << All[4] << setw(16) << All[5] << setw(16)
	 << All[6]  << setw(16) << All[7] << setw(16) << All[10] << setw(16)
	 << All[12] << setw(16) << All[8] << G4endl;



  myfile <<""<< G4endl;
  myfile <<""<< G4endl; 
  myfile <<"                                  Processes that produced a particle                               "<< G4endl;
  myfile <<""<< G4endl;
  myfile <<""<< G4endl;

  myfile << setw(40) << "" << setw(20) << partNames[0] << setw(16) << partNames[1] << setw(16) << partNames[2] << setw(16) 
	 << partNames[3]  << setw(16) << partNames[4] << setw(16) << partNames[5] << setw(16)
	 << partNames[6]  << setw(16) << partNames[7] << setw(16) << partNames[10] << setw(16)
	 << partNames[12] << setw(16) << partNames[8] << G4endl;

  for (int j=0; j<NProc; j++){
    myfile << setw(40) << processNames[j] << setw(20) << ProcessProducedParticle[0][j] << setw(16) << ProcessProducedParticle[1][j] << setw(16) << ProcessProducedParticle[2][j] << setw(16)
	   << ProcessProducedParticle[3][j]  << setw(16) << ProcessProducedParticle[4][j] << setw(16) << ProcessProducedParticle[5][j] << setw(16)
	   << ProcessProducedParticle[6][j]  << setw(16) << ProcessProducedParticle[7][j] << setw(16) << ProcessProducedParticle[10][j] << setw(16)
	   << ProcessProducedParticle[12][j] << setw(16) << ProcessProducedParticle[8][j] << G4endl;
  }


  myfile <<""<< G4endl;
  myfile <<"                                  Processes produced sum                               "<< G4endl;
  myfile <<""<< G4endl;

  // determining the sum of all processes that produced a particular particle
  G4double process_sum [NProc] = {0};

  for (int d=0; d<=NPart; d++){
    for (int proc=0; proc<=NProc; proc++){
      process_sum[d] = process_sum[d] + ProcessProducedParticle[d][proc];
    }
  }


  myfile << setw(40) << "" << setw(20) << process_sum[0] << setw(16) << process_sum[1] << setw(16) << process_sum[2] << setw(16) << process_sum[3] << setw(16) 
	 << process_sum[4]  << setw(16) << process_sum[5] << setw(16) << process_sum[6] << setw(16) << process_sum[7] << setw(16)
	 << process_sum[10] << setw(16) << process_sum[12] << setw(16) << process_sum[8] << G4endl;


  myfile <<""<< G4endl;
  myfile <<""<< G4endl;
  myfile <<"                                  Parents that produced a particle                               "<< G4endl;
  myfile <<""<< G4endl;
  myfile <<""<< G4endl;

  myfile << setw(40) << "" << setw(20) << partNames[0] << setw(16) << partNames[1] << setw(16) << partNames[2] << setw(16) 
	 << partNames[3]  << setw(16) << partNames[4] << setw(16) << partNames[5] << setw(16)
	 << partNames[6]  << setw(16) << partNames[7] << setw(16) << partNames[10] << setw(16)
	 << partNames[12] << setw(16) << partNames[8] << G4endl;
  for (int j=0; j<NParents; j++){
    myfile << setw(40) << parentNames[j] << setw(20) << ParentParticles[0][j] << setw(16) << ParentParticles[1][j] << setw(16) << ParentParticles[2][j] << setw(16)
	   << ParentParticles[3][j]  << setw(16) << ParentParticles[4][j] << setw(16) << ParentParticles[5][j] << setw(16)
	   << ParentParticles[6][j]  << setw(16) << ParentParticles[7][j] << setw(16) << ParentParticles[10][j] << setw(16)
	   << ParentParticles[12][j] << setw(16) << ParentParticles[8][j] << G4endl;
  }

  myfile <<""<< G4endl;
  myfile <<"                                  Parents sum                               "<< G4endl;
  myfile <<""<< G4endl;

  // determining the sum of all parents of particular particles
  G4double parents_sum [NPart] = {0};

  for (int d=0; d<=NPart; d++){
    for (int p=0; p<NParents; p++){
      parents_sum[d] = parents_sum[d] + ParentParticles[d][p];
    }
  }

  myfile << setw(40) << "" << setw(20) << parents_sum[0]  << setw(16) << parents_sum[1]  << setw(16) << parents_sum[2] << setw(16) << parents_sum[3]  << setw(16) 
	 << parents_sum[4]  << setw(16) << parents_sum[5]  << setw(16) << parents_sum[6] << setw(16) << parents_sum[7]  << setw(16)
	 << parents_sum[10] << setw(16) << parents_sum[12] << setw(16) << parents_sum[8] << G4endl;


  myfile <<""<< G4endl;
  myfile <<""<< G4endl;
  myfile <<"                                  Decays that produced a particle                               "<< G4endl;
  myfile <<""<< G4endl;
  myfile <<""<< G4endl;


  myfile <<""<< G4endl;
  myfile <<"                                  Decay At Rest                               "<< G4endl;
  myfile <<""<< G4endl;

  myfile << setw(40) << "" << setw(20) << partNames[0] << setw(16) << partNames[1] << setw(16) << partNames[2]  << setw(16) 
	 << partNames[3] << setw(16) << partNames[4] << setw(16) << partNames[5]  << setw(16)
	 << partNames[6] << setw(16) << partNames[7] << setw(16) << partNames[10] << setw(16)
	 << partNames[12] << G4endl;

  for (int j=0; j<NParents; j++){
    myfile << setw(40) << parentNames[j] << setw(20) << DecayProducedParticle[0][0][j] << setw(16) << DecayProducedParticle[0][1][j] << setw(16) << DecayProducedParticle[0][2][j]  << setw(16)
	   << DecayProducedParticle[0][3][j] << setw(16) << DecayProducedParticle[0][4][j] << setw(16) << DecayProducedParticle[0][5][j]  << setw(16)
	   << DecayProducedParticle[0][6][j] << setw(16) << DecayProducedParticle[0][7][j] << setw(16) << DecayProducedParticle[0][10][j] << setw(16)
	   << DecayProducedParticle[0][12][j] << G4endl;
  }



  myfile <<""<< G4endl;
  myfile <<""<< G4endl;

  myfile <<""<< G4endl;
  myfile <<"                                  Decay In Flight                               "<< G4endl;
  myfile <<""<< G4endl;

  myfile << setw(40) << "" << setw(20) << partNames[0] << setw(16) << partNames[1] << setw(16) << partNames[2] << setw(16) 
	 << partNames[3] << setw(16) << partNames[4] << setw(16) << partNames[5] << setw(16)
	 << partNames[6] << setw(16) << partNames[7] << setw(16) << partNames[10] << setw(16)
	 << partNames[12] << G4endl;

  for (int j=0; j<NParents; j++){
    myfile << setw(40) << parentNames[j] << setw(20) << DecayProducedParticle[1][0][j] << setw(16) << DecayProducedParticle[1][1][j] << setw(16) << DecayProducedParticle[1][2][j] << setw(16)
	   << DecayProducedParticle[1][3][j] << setw(16) << DecayProducedParticle[1][4][j] << setw(16) << DecayProducedParticle[1][5][j] << setw(16)
	   << DecayProducedParticle[1][6][j] << setw(16) << DecayProducedParticle[1][7][j] << setw(16) << DecayProducedParticle[1][10][j] << setw(16)
	   << DecayProducedParticle[1][12][j] << G4endl;
  }

  // determining the sum of all decays that produced a particular particle
  G4double decay_sum [NPart] = {0};

  for (int d=0; d<=NPart; d++){
    for (int p=0; p<NParents; p++){
      decay_sum[d] = decay_sum[d] + DecayProducedParticle[1][d][p] + DecayProducedParticle[0][d][p];
    }
  }


  myfile <<""<< G4endl;
  myfile <<"                                  DIF + DAR Sums                               "<< G4endl;
  myfile <<""<< G4endl;


  myfile << setw(40) << "" << setw(20)  << decay_sum[0] << setw(16)  << decay_sum[1] << setw(16) << decay_sum[2]  << setw(16) 
	 << decay_sum[3] << setw(16)  << decay_sum[4] << setw(16) << decay_sum[5]  << setw(16) 
	 << decay_sum[6] << setw(16)  << decay_sum[7] << setw(16) << decay_sum[10] << setw(16) 
	 << decay_sum[12] << G4endl;


  myfile <<""<< G4endl;
  myfile <<""<< G4endl;


  myfile <<""<< G4endl;
  myfile <<"                                  mu- Capture instances and products                              "<< G4endl;
  myfile <<""<< G4endl;

  myfile << setw(40) << " " << setw(10) << "#" << setw(20) << muCapProductNames [0] << setw(16) << muCapProductNames [1] << setw(16) << muCapProductNames [2] << setw(16)
	 << muCapProductNames [3] << setw(16) << muCapProductNames [4] << setw(16) << muCapProductNames [5] << setw(16)
	 << muCapProductNames [6] << G4endl;


  for (int inst=0; inst < muCapInstance; inst++){
    myfile << setw(40) << "Instance" << setw(10) << inst << setw(20) << mucap_secondaries[inst][3] << setw(16) << mucap_secondaries[inst][2] << setw(16) 
	   << mucap_secondaries[inst][1] << setw(16) << mucap_secondaries[inst][4] << setw(16) 
	   << mucap_secondaries[inst][0] << setw(16) << mucap_secondaries[inst][5] << setw(45)
	   << G4endl;
  }


  myfile <<""<< G4endl;
  myfile <<""<< G4endl;
  myfile <<""<< G4endl;
  myfile <<"                                  mu- Capture products sum                              "<< G4endl;
  myfile <<""<< G4endl;

  // determining the sum of mu- capture products
  G4double mcapture_sums [6] = {0};

  for (int d=0; d<6; d++){
    for (int inst=0; inst<muCapInstance; inst++){
      mcapture_sums[d] = mcapture_sums[d] + mucap_secondaries[inst][d];
    }
  }


  myfile << setw(40) << "" << setw(10) << "" << setw(20) << mcapture_sums[3]  << setw(16) << mcapture_sums[2] << setw(16) << mcapture_sums[1] << setw(16)
	 << mcapture_sums[4] << setw(16) << mcapture_sums[0]  << setw(16) << mcapture_sums[5]   << setw(16)
	 <<  "sum_capture_z" << G4endl;


  myfile <<""<< G4endl;
  myfile <<""<< G4endl;

  myfile <<"                                  mu- Capture instances and products (strings)                             "<< G4endl;

  myfile << setw(40) << " " << setw(10) << "#" << setw(20) << G4endl;

  for (int j=0; j <= muCapInstance; j++){
    myfile << setw(40) << "Instance" << setw(10) << j << setw(10) << "     " << muCapZ [j] << G4endl;
  }

  myfile <<""<< G4endl;
  myfile <<""<< G4endl;
  myfile <<""<< G4endl;



  myfile <<""<< G4endl;
  myfile <<"                                  This section provides a digestable breakdown of the above raw simulation numbers  "<< G4endl;
  myfile <<""<< G4endl;


  myfile <<""<< G4endl;
  myfile <<"                                  Percent breakdown of neutrino and pion parents                               "<< G4endl;
  myfile <<""<< G4endl;


  myfile << setw(40) << "" << setw(20) << partNames[0] << setw(16) << partNames[1] << setw(16) << partNames[4] << setw(16) 
	 << partNames[5] << setw(16) << partNames[6] << setw(16) << partNames[7] << G4endl;



  for (int j=0; j<NParents; j++){
    myfile << setw(40) << parentNames[j] << setw(20) << ParentParticles[0][j]/parents_sum[0]*100 << setw(1) << "%" << setw(16) 
	   << ParentParticles[1][j]/parents_sum[1]*100 << setw(1) << "%" << setw(16) 
	   << ParentParticles[4][j]/parents_sum[4]*100 << setw(1) << "%" << setw(16)
	   << ParentParticles[5][j]/parents_sum[5]*100 << setw(1) << "%" << setw(16) 
	   << ParentParticles[6][j]/parents_sum[6]*100 << setw(1) << "%" << setw(16) 
	   << ParentParticles[7][j]/parents_sum[7]*100 << setw(1) << "%" << G4endl;
  }


  myfile <<""<< G4endl;
  myfile <<""<< G4endl;


  myfile <<""<< G4endl;
  myfile <<"                                  Percent breakdown of neutrinos from pion & muon DAR/DIF                              "<< G4endl;
  myfile <<""<< G4endl;



  myfile <<""<< G4endl;
  myfile <<"                                  Decay At Rest                               "<< G4endl;
  myfile <<""<< G4endl;

  myfile << setw(40) << "" << setw(20) << partNames[4] << setw(16) 
	 << partNames[5] << setw(16) << partNames[6] << setw(16) << partNames[7] << G4endl;


  for (int j=5; j<=9; j++){
    myfile << setw(40) << parentNames[j] << setw(20) << DecayProducedParticle[0][4][j]/decay_sum[4]*100  << setw(1) << "%" << setw(16)
	   << DecayProducedParticle[0][5][j]/decay_sum[5]*100  << setw(1) << "%" << setw(16)
	   << DecayProducedParticle[0][6][j]/decay_sum[6]*100  << setw(1) << "%" << setw(16) 
	   << DecayProducedParticle[0][7][j]/decay_sum[7]*100  << setw(1) << "%" << G4endl;
  }


  myfile <<""<< G4endl;
  myfile <<""<< G4endl;

  myfile <<""<< G4endl;
  myfile <<"                                  Decay In Flight                               "<< G4endl;
  myfile <<""<< G4endl;

  myfile << setw(40) << "" << setw(20) << partNames[4] << setw(16) 
	 << partNames[5] << setw(16) << partNames[6] << setw(16) << partNames[7] << G4endl;

  for (int j=5; j<=9; j++){
    myfile << setw(40) << parentNames[j] << setw(20) << DecayProducedParticle[1][4][j]/decay_sum[4]*100       << setw(1) << "%" << setw(16)
	   << DecayProducedParticle[1][5][j]/decay_sum[5]*100  << setw(1) << "%" << setw(16) 
	   << DecayProducedParticle[1][6][j]/decay_sum[6]*100        << setw(1) << "%" << setw(16) 
	   << DecayProducedParticle[1][7][j]/decay_sum[7]*100   << setw(1) << "%" << G4endl;
  }

  myfile <<""<< G4endl;
  myfile <<""<< G4endl;
  myfile <<""<< G4endl;


  myfile <<""<< G4endl;
  myfile <<"                                  Destructive Processes a Particle Engaged In                              "<< G4endl;
  myfile <<""<< G4endl;
  myfile <<" NOTE: Some of these processes may have also resulted in the creation of a particle, but they do all result in the loss of the primary/parent particle,                       "<< G4endl;
  myfile <<"       so this list is just the counted number of particles (according to each particle type) that did die as a result of these interactions.                  "<< G4endl;
  myfile <<""<< G4endl;

  myfile << setw(40) << "" << setw(20) << parentNames[5] << setw(16) << parentNames[6] << setw(16) << parentNames[8] << setw(16) << parentNames[9] << setw(16) 
	 << parentNames[2] << G4endl;

  for (int proc=20; proc<=56; proc++){

    if (proc == 28){myfile << setw(40) << processNames[proc] << setw(20) << Processes[5][proc] << setw(16) << Processes[6][proc]  << setw(16) 
                                                      << Processes[8][proc] << setw(16) << Processes[9][proc]  << setw(16)
	                                                    << Processes[2][proc] << G4endl;}

    else if (proc != 28){myfile << setw(40) << processNames[proc] << setw(20) << Processes[5][proc] - FailedProcesses[5][proc] << setw(16) << Processes[6][proc] - FailedProcesses[6][proc]  << setw(16) 
                                                         << Processes[8][proc] - FailedProcesses[8][proc] << setw(16) << Processes[9][proc] - FailedProcesses[9][proc]  << setw(16)
	                                                       << Processes[2][proc] - FailedProcesses[2][proc] << G4endl;}

  }


  // determining the sum of all destructive processes that a particle engaged in
  G4double destructive_sum [NPart] = {0};
  
  for (int p=0; p<=NParents; p++){
    for (int proc=20; proc<NProc; proc++){
      destructive_sum[p] = destructive_sum[p] + Processes[p][proc];
      if (proc != 28){destructive_sum[p] = destructive_sum[p] - FailedProcesses[p][proc];}
    }
  }


  myfile <<""<< G4endl;
  myfile <<"                                  Destructive Processes Sum                              "<< G4endl;
  myfile <<""<< G4endl;

  myfile << setw(40) << "" << setw(20) << destructive_sum[5] << setw(16) << destructive_sum[6]  << setw(16) << destructive_sum[8] << setw(16)
	 << destructive_sum[9] << setw(16) << destructive_sum[2]  << G4endl;

  myfile <<""<< G4endl;
  myfile <<""<< G4endl; 
  myfile <<"                                 (Particle Produced By Parent Sums) - (Destructive Processes Sums For That Particle)                           "<< G4endl;
  myfile <<""<< G4endl;

  myfile << setw(40) << "" << setw(20) << (parents_sum[0] - destructive_sum[5])  << setw(16) << (parents_sum[1] - destructive_sum[6])  << setw(16) << (parents_sum[2] - destructive_sum[8]) << setw(16)
	 << (parents_sum[3] - destructive_sum[9])  << setw(16) << (parents_sum[8] - destructive_sum[2]) << G4endl;

  myfile <<""<< G4endl;
  myfile <<""<< G4endl; 
  myfile <<"                                 (Particle Produced Process) - (Destructive Processes Sums For That Particle)                           "<< G4endl;
  myfile <<""<< G4endl;

  myfile << setw(40) << "" << setw(20) << (process_sum[0] - destructive_sum[5])  << setw(16) << (process_sum[1] - destructive_sum[6])  << setw(16) << (process_sum[2] - destructive_sum[8]) << setw(16)
	 << (process_sum[3] - destructive_sum[9])  << setw(16) << (process_sum[8] - destructive_sum[2]) << G4endl;


  myfile <<""<< G4endl; 
  myfile <<""<< G4endl; 

/*
  myfile <<""<< G4endl;
  myfile <<""<< G4endl;
  myfile <<"                                  Processes a particle engaged in with Discrimination for Secondaries Produced                                                "<< G4endl;
  myfile <<""<< G4endl;
  myfile <<""<< G4endl;

  myfile <<"                                  Processes that did not Produce Secondaries                                                "<< G4endl;
  myfile <<""<< G4endl;
  myfile <<""<< G4endl;

  myfile << setw(40) << "" << setw(20) << parentNames[5] << setw(16) << parentNames[6] << setw(16) << parentNames[2] << setw(16) 
	 << parentNames[8] << setw(16) << parentNames[9] << setw(16) << parentNames[10] << G4endl;

  for (int j=0; j<NProc; j++){
    myfile << setw(40) << processNames[j] << setw(20) << ProcessesTracking[5][j][0] << setw(16) << ProcessesTracking[6][j][0] << setw(16) << ProcessesTracking[2][j][0] << setw(16)
	   << ProcessesTracking[8][j][0] << setw(16) << ProcessesTracking[9][j][0] << setw(16) << ProcessesTracking[10][j][0] << G4endl;
  }

  myfile <<""<< G4endl;
  myfile <<""<< G4endl;
  myfile <<"                                  Proccesses that Produced EXACTLY 1 Secondaries                                                "<< G4endl;
  myfile <<""<< G4endl;
  myfile <<""<< G4endl;

  for (int j=0; j<NProc; j++){
    myfile << setw(40) << processNames[j] << setw(20) << ProcessesTracking[5][j][1] << setw(16) << ProcessesTracking[6][j][1] << setw(16) << ProcessesTracking[2][j][1] << setw(16)
	   << ProcessesTracking[8][j][1] << setw(16) << ProcessesTracking[9][j][1] << setw(16) << ProcessesTracking[10][j][1] << G4endl;
  }

  myfile <<""<< G4endl;
  myfile <<""<< G4endl;
  myfile <<"                                  Proccesses that Produced more than 1 Secondaries                                                "<< G4endl;
  myfile <<""<< G4endl;
  myfile <<""<< G4endl;

  for (int j=0; j<NProc; j++){
    myfile << setw(40) << processNames[j] << setw(20) << ProcessesTracking[5][j][2] << setw(16) << ProcessesTracking[6][j][2] << setw(16) << ProcessesTracking[2][j][2] << setw(16)
	   << ProcessesTracking[8][j][2] << setw(16) << ProcessesTracking[9][j][2] << setw(16) << ProcessesTracking[10][j][2] << G4endl;
  }



  myfile <<""<< G4endl;
  myfile <<"                                  Processes Tracking Net                               "<< G4endl;
  myfile <<""<< G4endl;

  // determing the sum of processes a particle engaged in
  G4double process_engaged_tracking_sum [NParents] = {0};

  for (int p=0; p<=NPart; p++){
    for (int proc=0; proc<=NProc; proc++){
      process_engaged_tracking_sum[p] = process_engaged_tracking_sum[p] - ProcessesTracking[p][proc][0] + ProcessesTracking[p][proc][1];
    }
  }


  myfile << setw(40) << "" << setw(20) << process_engaged_tracking_sum[5] << setw(16) << process_engaged_tracking_sum[6]  << setw(16) << process_engaged_tracking_sum[2] << setw(16) 
                                       << process_engaged_tracking_sum[8] << setw(16) << process_engaged_tracking_sum[9] << setw(16)  << process_engaged_tracking_sum[10] << setw(16) << G4endl;




*/

  myfile <<""<< G4endl; 
  myfile <<""<< G4endl; 
  myfile <<""<< G4endl;
  myfile <<""<< G4endl; 
  myfile <<"                                  Decay Reactions                             "<< G4endl;
  myfile <<""<< G4endl;

  myfile <<""<< G4endl;
  myfile <<"                                  pi+   ->   mu+   nu_mu   ->   e+   nu_e   anti_nu_mu                            "<< G4endl;
  myfile <<""<< G4endl;
  myfile <<"                                  *Includes counters for the alternate decays of mu+, see Physical Review D                "<< G4endl;
  myfile <<"                                  **If e+, nu_e, anti_nu_mu > pi+, mu+, nu_mu, check to see if K+ -> mu+, or {any} -> mu+               "<< G4endl;
  myfile <<""<< G4endl;
  myfile << setw(40) << "" << setw(20) << "pi+"            << setw(16) <<  "mu+" << setw(16) << "nu_mu" <<  setw(16) << "e-" <<  setw(16) << "gamma" <<  setw(16) 
	 << "e+" << setw(16) << "nu_e" << setw(16) << "anti_nu_mu" << G4endl;
  myfile << setw(40) << "" << setw(20) << Processes[5][56] << setw(16) << DecayProducedParticle[1][2][5]  + DecayProducedParticle[0][2][5]  << setw(16) 
	 << DecayProducedParticle[1][4][5]  + DecayProducedParticle[0][4][5]  << setw(16)
	 << DecayProducedParticle[1][10][8] + DecayProducedParticle[0][10][8] << setw(16)
	 << DecayProducedParticle[1][11][8] + DecayProducedParticle[0][11][8] << setw(16)
	 << DecayProducedParticle[1][12][8] + DecayProducedParticle[0][12][8] << setw(16)
	 << DecayProducedParticle[1][6][8]  + DecayProducedParticle[0][6][8]  << setw(16)
	 << DecayProducedParticle[1][5][8]  + DecayProducedParticle[0][5][8]  << G4endl;
  myfile <<""<< G4endl; 
  myfile <<""<< G4endl; 
  myfile <<"                                  pi-   ->   mu-   anti_nu_mu   ->   e-   anti_nu_e   nu_mu                            "<< G4endl;
  myfile <<""<< G4endl;
  myfile << setw(40) << "" << setw(20) << "pi-"            << setw(16) <<  "mu-" << setw(16) << "anti_nu_mu" <<  setw(16) << "e-" << setw(16) << "anti_nu_e" << setw(16) << "nu_mu" << G4endl;
  myfile << setw(40) << "" << setw(20) << Processes[6][56] << setw(16) << DecayProducedParticle[1][3][6]  + DecayProducedParticle[0][3][6]  << setw(16) 
	 << DecayProducedParticle[1][3][6]  + DecayProducedParticle[0][3][6]  << setw(16)
	 << DecayProducedParticle[1][10][9] + DecayProducedParticle[0][10][9] << setw(16)
	 << DecayProducedParticle[1][7][9]  + DecayProducedParticle[0][7][9]  << setw(16)
	 << DecayProducedParticle[1][4][9]  + DecayProducedParticle[0][4][9]  << G4endl;
  myfile <<""<< G4endl; 
  myfile <<"                                 *Checking number of mu- captured               "<< G4endl;
  myfile <<""<< G4endl; 
  myfile << setw(40) << "" << setw(20) << "muMinusCaptureAtRest" << setw(16) << Processes[9][45] << setw(16) << G4endl;
  myfile <<""<< G4endl; 
  myfile <<""<< G4endl; 



  myfile <<""<< G4endl; 
  myfile <<"                                 # X secondaries from p+Hg / # protons on target              "<< G4endl;
  myfile <<""<< G4endl; 

  myfile << setw(40) << "" << setw(20) << partNames[0] << setw(16) << partNames[1] << setw(16) << partNames[8] << setw(16) 
	 << partNames[9] << G4endl;

  myfile << setw(40) << "" << setw(20) << ProtonPlusHgProduction[0]/ProtonsOnTarget << setw(16) << ProtonPlusHgProduction[1]/ProtonsOnTarget 
                           << setw(16) << ProtonPlusHgProduction[8]/ProtonsOnTarget << setw(16) << ProtonPlusHgProduction[9]/ProtonsOnTarget << G4endl;

  myfile <<""<< G4endl; 
  myfile <<""<< G4endl; 

  myfile <<""<< G4endl; 
  myfile <<""<< G4endl; 












  myfile_1 <<""<< G4endl; 
  myfile_1 <<""<< G4endl; 
  myfile_1 <<""<< G4endl; 
  myfile_1 <<"                                  pi+Inelastic Analysis                      "<< G4endl;
  myfile_1 <<""<< G4endl;
  myfile_1 <<""<< G4endl;
  myfile_1 <<""<< G4endl;

  //determing the sum of all pi+'s from pi+Inelastic
  G4double sum_pi_plus_inelastic_produced=0;

  for (int j=0; j<=NTracks; j++){
    sum_pi_plus_inelastic_produced = sum_pi_plus_inelastic_produced + TrackData[j];
  }


  // determing the number of instances when pi+Inelastic produced single or multiple pi+'s
  G4double multiplicity_pi [10] = {0};

  for (int j=0; j<=NTracks; j++){
    if(TrackData[j] == 1){multiplicity_pi [0] ++;}
    if(TrackData[j] == 2){multiplicity_pi [1] ++;}
    if(TrackData[j] == 3){multiplicity_pi [2] ++;}
    if(TrackData[j] == 4){multiplicity_pi [3] ++;}
    else{multiplicity_pi [4] ++;}
  }


  // determining the number of non-zero entries of the track data array
  G4double nonZeroTracks=0;

  for (int j=0; j<=NTracks; j++){
    if (TrackData[j] != 0){nonZeroTracks++;}
  }


  myfile_1 << setw(80) << "# of times pi+Inelastic called for pi+ =" << setw(20) << Processes[5][20] << G4endl;
  myfile_1 << setw(80) << "(outside secondary loop, OL)" << setw(20) << G4endl;
  myfile_1 <<""<< G4endl; 
  myfile_1 <<""<< G4endl; 
  myfile_1 << setw(80) << "# of pi+Inelastic Absolute Deaths" << setw(20) << G4endl;
  myfile_1 << setw(80) << "(pi+Inelastic -> no pi+, but something else) =" << setw(20) << DestroyedParticle [0] << G4endl;
  myfile_1 <<""<< G4endl; 
  myfile_1 <<""<< G4endl; 
  myfile_1 << setw(80) << "# of pi+ Produced by pi+Inelastic (orig. code) =" << setw(20) << ProcessProducedParticle[0][20] << G4endl;
  myfile_1 <<""<< G4endl; 
  myfile_1 << setw(80) << "sum of array elements for (pi+Inelastic -> N pi+) =" << setw(20) <<  sum_pi_plus_inelastic_produced << G4endl;
  myfile_1 <<""<< G4endl; 
  myfile_1 <<""<< G4endl; 
  myfile_1 << setw(80) << "# Instances where (pi+Inelastic ->  multiple pi+'s) =" << setw(20) << (nonZeroTracks - multiplicity_pi[0]) << G4endl;
  myfile_1 <<""<< G4endl; 
  myfile_1 << setw(80) << "# Instances where (pi+Inelastic ->  1 pi+'s) =" << setw(20) << multiplicity_pi[0] << G4endl;
  myfile_1 <<""<< G4endl; 
  myfile_1 << setw(80) << "# Instances where (pi+Inelastic ->  2 pi+'s) =" << setw(20) << multiplicity_pi[1] << G4endl;
  myfile_1 <<""<< G4endl; 
  myfile_1 << setw(80) << "# Instances where (pi+Inelastic ->  3 pi+'s) =" << setw(20) << multiplicity_pi[2] << G4endl;
  myfile_1 <<""<< G4endl; 
  myfile_1 <<""<< G4endl; 
  myfile_1 << setw(80) << "# of Non-zero entries in TrackData array" << setw(20) << G4endl;
  myfile_1 << setw(80) << "i.e., total # of times (pi+ -> any # of pi+'s) =" << setw(20) << nonZeroTracks << G4endl;
  myfile_1 <<""<< G4endl; 
  myfile_1 <<""<< G4endl; 
  myfile_1 << setw(80) << "# of times pi+Inelastic successful =" << setw(20) << nonZeroTracks + DestroyedParticle[0] << G4endl;
  myfile_1 << setw(80) << "(inside secondary loop, IL)" << setw(20) << G4endl;
  myfile_1 <<""<< G4endl; 
  myfile_1 <<""<< G4endl; 

  // this calculation should give the number of unsuccessful pi+Inelastic events
  myfile_1 << setw(80) << "OL - IL =" << setw(20) << Processes[5][20] - nonZeroTracks - DestroyedParticle[0] << G4endl;




}// end of run action




