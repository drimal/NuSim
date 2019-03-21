#include "SNSRootManager.hh"
#include <iostream>
#include <stdexcept>
#include <sstream>

#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>

/*
IGNORE FOR NOW, LEAVING HERE ONLY FOR REFERENCE 

// emission energy of parent decays that produced a particle
// inside secondary loop
TH1F *Energy_DAR_DIF_[(NDecayTypes)+1][(NParents)][(NPart)];

// need to adjust for only emission energy of daughters
TH1F *Energy_Emission[(NPart)];

*/


// declaring some things globally

// dimensioning integers for histogram arrays
const int NPart=13;
const int NProc=57;
const int NParents=33;
const int NTracks=100000;
const int NmuMinusCaptures=70000;
const int NDecayTypes=2;
const int NDestructProc=26;
const int NDet=4;
const int NNeutrinoTypes=4;
const int NPrimaries=10;

// 1D histogram counting all the parents of specific particles
TH1F *ParentParticles_[(NPart)];

// 1D histogram counting all the processes that produced specific particles
TH1F *ProcessProducedParticle_[(NPart)];

// 1D histogram recording MCNP ID's
TH1F *MCNPID_[(NPart)];
TH1F *MCNPID_L[(NPart)];

// 1D histogram of destructive process histograms
TH1F *Kill_[(NPrimaries)];

// 2D histogram of energy of primary that DIF, and end of track lifetime of that primary
TH2F *Energy_DIF[(NPrimaries)];

// 2D histogram of the initial time and energy of all emitted secondaries
TH2F *Energy_Time_Creation[(NPart)];

// 1D histogram of initial energy of emitted secondaries
TH1F *Energy_initial[(NPart)];

// 1D histogram of track lifetimes of particles
TH1F *Time_ofLife[(NPrimaries)];

// 1D histogram of the time of emission of secondaries relative to the p+Hg interaction (initial time = 0)
TH1F *Time_atEmission[(NPart)];

// 1D histogram of the track lifetime of particles that die from decay only, DAR/DIF seperated
TH1F *Time_ofLife_tillDecayed[(NDecayTypes)][(NPrimaries)];

// 1D histogram of the track lifetime of particles that die from decay only, any decay type, DIF + DAR
TH1F *Time_ofLife_tillAnyDecay[(NPrimaries)];

// 1D histogram of the time when particles DAR/DIF seperated relative to p+Hg interaction, GlobalTime
TH1F *Time_elapsedSince_ProtonHg_Interaction[(NDecayTypes)][(NPrimaries)];

// 1D histogram of the time when particles DAR + DIF, relative to p+Hg interaction, GlobalTime
TH1F *Time_elapsedSince_ProtonHg_Interaction_AnyDecay[(NPrimaries)];

// 1D histogram of the arrival times of neutrinos in the detector volumes
TH1F *Time_Nu[(NDet)][(NNeutrinoTypes)];

// 1D histogram of the energy of neutrinos that make it to the detector volumes
TH1F *Energy_Nu[(NDet)][(NNeutrinoTypes)];

// 1D histograms of the energy of a primary that DIF at the 1st step and last step
TH1F *Energy_Primary_DIF_1Step[(NPrimaries)];
TH1F *Energy_Primary_DIF_KillStep[(NPrimaries)];

// 1D histograms of the energy of all particles that DAR/DIF seperated
TH1F *Energy_Primaries_DAR_DIF[(NDecayTypes)][(NPrimaries)];

// 1D histograms of the energy of all particles that decay, any decay type, DIF + DAR
TH1F *Energy_Primaries_AnyDecay[(NPrimaries)];
TH2F *h2d;
TH3F *h3d;
TH1F *h1d;

SNSRootManager* SNSRootManager::fgInstance = 0;

SNSRootManager::SNSRootManager(string fout)
  : fOutFileName(fout)
{
  if (fgInstance) {
    std::ostringstream message;
    message << "SNSRootManager::SNSRootManager: singleton instance already exists.";
    throw std::runtime_error(message.str());
    return;
  } 
  fHlist = new TList();
  fgInstance = this;
}

SNSRootManager::~SNSRootManager() 
{
  fgInstance = 0;
}

void SNSRootManager::Init()
{

  // initialize the histo arrays in seperate loops, so that they appear in ordered groups
    h3d = new TH3F("h2dxyz", "x_y_z; x [cm]; y [cm]; z[cm]", 200, -25, 25, 200, -20, 20, 200, -60, 60);
    h3d->SetDirectory(0);fHlist->Add(h3d);
    
    h2d = new TH2F("h2dxy", "x_y; x [cm]; y [cm]", 500, -60, 60, 500, -20, 20);
    h2d->SetDirectory(0);fHlist->Add(h2d);
    
    h1d = new TH1F("h1", "distance; z [cm]; counts", 500, -50, 50);
    h1d->SetDirectory(0);fHlist->Add(h1d);

    

  for (int d=0; d<NPart; d++){
    // processes that produced a particle histos
    ProcessProducedParticle_ [d] = new TH1F(Form("processes_created_%s", partNames[d]),"created; process; counts",57,-1,57);
    ProcessProducedParticle_ [d] -> SetDirectory(0); fHlist->Add(ProcessProducedParticle_[d]);
  }


  for (int d=0; d<NPart; d++){
    // MCNP ID histos
    MCNPID_ [d] = new TH1F(Form("mcnp_id_for_parents_of_%s", partNames[d]),"parent ID's; ID; Counts",5000, -2500,2500);
    MCNPID_ [d] -> SetDirectory(0); fHlist->Add(MCNPID_ [d]);
  }


  for (int d=0; d<NPart; d++){
    // Long MCNP ID histos
    MCNPID_L [d] = new TH1F(Form("long_mcnp_id_for_parents_of_%s", partNames[d]),"parent ID's; ID; Counts",60000,1000000000,1000060000);
    MCNPID_L [d] -> SetDirectory(0); fHlist->Add(MCNPID_L [d]);
  }


  for (int p=0; p<NPrimaries; p++){
    // Destructive process histo
    Kill_ [p] = new TH1F(Form("processes_destroyed_%s", primaryNames[p]),"Destructive processes; Process; Counts",58,-1,57);
    Kill_ [p] -> SetDirectory(0); fHlist->Add(Kill_ [p]);
  }


  for (int p=0; p<NPrimaries; p++){
    //2D histogram of energy of primaries that DIF, and end of track lifetime of that primary
    // need for pi+/-, mu+/- absolutely
    Energy_DIF[p] = new TH2F(Form("DIF_energy_lifetime_of_%s", primaryNames[p]),"Energy vs. Lifetime; ns; MeV",10000,0,10000,700,0,700);
    Energy_DIF[p] -> SetDirectory(0); fHlist->Add(Energy_DIF [p]);
  }


  for (int d=0; d<NPart; d++){
    // 2D histogram of the initial time and energy of all emitted secondaries
    // n, pi+/-, mu+/-, all nu flavors
    Energy_Time_Creation[d] = new TH2F(Form("initial_energy_time_of_%s", partNames[d]),"Initial Energy vs. Creation time; ns; MeV",1000,0,1000, 700,0, 700);
    Energy_Time_Creation[d] -> SetDirectory(0); fHlist->Add(Energy_Time_Creation [d]);
  }

  for (int d=0; d<NPart; d++){
    // 1D histogram of the initial energy of emitted secondaries
    Energy_initial[d] = new TH1F(Form("initial_energy_of_%s", partNames[d]),"Initial Energy; MeV; Counts",700,0,700);
    Energy_initial[d] -> SetDirectory(0); fHlist->Add(Energy_initial [d]);
  }

  for (int p=0; p<NPrimaries; p++){
    // 1D histogram of total particle lifetimes
    Time_ofLife[p] = new TH1F(Form("lifetime_of_%s", primaryNames[p]),"Lifetimes; ns; Counts",10000,0,10000);
    Time_ofLife[p] -> SetDirectory(0); fHlist->Add(Time_ofLife [p]);
  }

  for (int d=0; d<NPart; d++){
    // 1D histogram of the time of emission of secondaries relative to the p+Hg interaction (initial time = 0)
    Time_atEmission[d] = new TH1F(Form("creation_time_of_%s", partNames[d]),"Time created rel. to p+Hg; ns; Counts", 50000, 0, 20000);
    Time_atEmission[d] -> SetDirectory(0); fHlist->Add(Time_atEmission [d]);
  }


  for (int decay=0; decay<NDecayTypes; decay++){

    for (int p=0; p<NPrimaries; p++){ 

    // 1D histogram of lifetimes for particles that DAR/DIF (Time since created until decay, LocalTime)
    Time_ofLife_tillDecayed[decay][p] = new TH1F(Form("lifetime_%s_of_%s", decayNames[decay], primaryNames[p]), Form("%s decay lifetime; ns; Counts", primaryNames[p]),1000,0,10000);
    Time_ofLife_tillDecayed[decay][p] -> SetDirectory(0); fHlist->Add(Time_ofLife_tillDecayed[decay][p]);
    }
  }

    for (int p=0; p<NPrimaries; p++){ 
    // 1D histogram of lifetimes for all particles that Decay (OF ANY TYPE, DIF + DAR) (Time since created until decay, LocalTime)
    Time_ofLife_tillAnyDecay[p] = new TH1F(Form("lifetime_any_decay_of_%s", primaryNames[p]), Form("DAR + DIF decay lifetime; ns; Counts", primaryNames[p]),1000,0,10000);
    Time_ofLife_tillAnyDecay[p] -> SetDirectory(0); fHlist->Add(Time_ofLife_tillAnyDecay[p]);
    }

  for (int decay=0; decay<NDecayTypes; decay++){

    for (int p=0; p<NPrimaries; p++){ 
    // 1D histogram of the time when particles DAR/DIF seperated relative to p+Hg interaction, GlobalTime
    Time_elapsedSince_ProtonHg_Interaction[decay][p] = new TH1F(Form("global_existence_time_%s_of_%s", decayNames[decay], primaryNames[p]), Form("%s decay lifetime; ns; Counts", primaryNames[p]),1000,0,10000);
    Time_elapsedSince_ProtonHg_Interaction[decay][p] -> SetDirectory(0); fHlist->Add(Time_elapsedSince_ProtonHg_Interaction[decay][p]);
    }
  }

    for (int p=0; p<NPrimaries; p++){ 
    // 1D histogram of lifetimes for all particles that Decay (OF ANY TYPE, DIF + DAR) relative to the p+Hg interaction, GlobalTime
    Time_elapsedSince_ProtonHg_Interaction_AnyDecay[p] = new TH1F(Form("global_exsistence_time_any_decay_of_%s", primaryNames[p]), Form("DAR + DIF decay lifetime; ns; Counts", primaryNames[p]),1000,0,10000);
    Time_elapsedSince_ProtonHg_Interaction_AnyDecay[p] -> SetDirectory(0); fHlist->Add(Time_elapsedSince_ProtonHg_Interaction_AnyDecay[p]);
    }





// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  for (int decay=0; decay<NDecayTypes; decay++){

    for (int p=0; p<NPrimaries; p++){ 
    // 1D histogram of the energy when particles DAR/DIF seperated
    Energy_Primaries_DAR_DIF[decay][p] = new TH1F(Form("primary_energy_%s_of_%s", decayNames[decay], primaryNames[p]), Form("%s Energy; MeV; Counts", primaryNames[p]),1000,0,1000);
    Energy_Primaries_DAR_DIF[decay][p] -> SetDirectory(0); fHlist->Add(Energy_Primaries_DAR_DIF[decay][p]);
    }
  }

    for (int p=0; p<NPrimaries; p++){ 
    // 1D histogram of energy for all particles that Decay (OF ANY TYPE, DIF + DAR)
    Energy_Primaries_AnyDecay[p] = new TH1F(Form("primary_energy_any_decay_of_%s", primaryNames[p]), Form("energy DAR + DIF; MeV; Counts", primaryNames[p]),1000,0,1000);
    Energy_Primaries_AnyDecay[p] -> SetDirectory(0); fHlist->Add(Energy_Primaries_AnyDecay[p]);
    }
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




  for (int det=0; det<NDet; det++){

    for (int nu=0; nu<NNeutrinoTypes; nu++){
      // 1D histogram of the energy of neutrinos that make it to a detector volume
      Energy_Nu [det][nu] = new TH1F(Form("%s_det_all_%s_energy", volNames[det], nuNames[nu]),"Energies; MeV; Counts",1000,0,300);
      Energy_Nu [det][nu] -> SetDirectory(0); fHlist->Add(Energy_Nu [det][nu]);
    }

    for (int nu=0; nu<NNeutrinoTypes; nu++){
      // 1D histogram of the arrival times of the neutrinos that make it to a detector volume
      Time_Nu [det][nu] = new TH1F(Form("%s_det_all_%s_time", volNames[det], nuNames[nu]),"Emitted times; ns; Counts",1000,0,1000);
      Time_Nu [det][nu] -> SetDirectory(0); fHlist->Add(Time_Nu [det][nu]);
    }

  }


//TEMPORARY:
  for (int p=0; p<NPrimaries; p++){
    // 1D histograms of the energy of a primary that DIF at the 1st step
    Energy_Primary_DIF_1Step[p] = new TH1F(Form("DIF_1step_energy_of_%s", primaryNames[p]),"Energy; MeV; Counts",700,0, 700);
    Energy_Primary_DIF_1Step[p] -> SetDirectory(0); fHlist->Add(Energy_Primary_DIF_1Step [p]);
  }


  for (int p=0; p<NPrimaries; p++){
    // 1D histograms of the energy of a primary that DIF at the last step
    Energy_Primary_DIF_KillStep[p] = new TH1F(Form("DIF_killstep_energy_of_%s", primaryNames[p]),"Energy; MeV; Counts",700,0, 700);
    Energy_Primary_DIF_KillStep[p] -> SetDirectory(0); fHlist->Add(Energy_Primary_DIF_KillStep [p]);
  }


}//end of initialization


void SNSRootManager::Fill3D(Double_t x, Double_t y, Double_t z){
    h3d->Fill(x, y, z);
}

void SNSRootManager::Fill2D(Double_t x, Double_t y){
    h2d->Fill(x, y);
}

void SNSRootManager::Fill1D(Double_t x){
    h1d->Fill(x);
}

// methods for filling ROOT histo arrays
// need to add file error checkers, and some other debugging stuff

//TEMPORARY:
void SNSRootManager::FillEnergy_1STEP(Int_t primaryFlag, Double_t energy){ 
  Energy_Primary_DIF_1Step[primaryFlag] -> Fill(energy); 
}
//TEMPORARY:
void SNSRootManager::FillEnergy_KILLSTEP(Int_t primaryFlag, Double_t energy){ 
  Energy_Primary_DIF_KillStep[primaryFlag] -> Fill(energy); 
}


void SNSRootManager::FillProcess(Int_t partFlag, Int_t procFlag)
{
  int d = partFlag;
  int proc = procFlag;

  ProcessProducedParticle_[d] -> Fill(proc);

}


void SNSRootManager::FillMCNPID(Int_t partFlag, Int_t mcnpID)
{
  int d = partFlag;

  if (mcnpID > 5000){MCNPID_L[d] -> Fill(mcnpID);}
  else {MCNPID_[d] -> Fill(mcnpID);}

}


void SNSRootManager::FillDestroyed(Int_t parentFlag, Int_t procFlag){
  int p = parentFlag;
  int proc = procFlag;

  Kill_ [p] -> Fill(proc);

}



void SNSRootManager::Fill2DEnergyTime(Int_t particleFlag, Double_t energy, Double_t time, Int_t pType){

if (pType == 0){
   int p = particleFlag;
   Energy_DIF[p]->Fill(time, energy);
}

if (pType == 1){
  int d = particleFlag;
  Energy_Time_Creation[d]->Fill(time, energy);
}

}


void SNSRootManager::FillEnergy(Int_t partFlag, Double_t energy){ 
  Energy_initial [partFlag] -> Fill(energy); 
}


void SNSRootManager::FillTime(Int_t particleFlag, Double_t time, Int_t timeFlag){ 

// track lifetime of primaries, time recorded at last step
if (timeFlag == 0){
  Time_ofLife [particleFlag] -> Fill(time); 
}

// emission time of secondaries relative to p+Hg (initTime=0)
if (timeFlag == 1){
  Time_atEmission [particleFlag] -> Fill(time); 
}

}

void SNSRootManager::FillTime(Int_t primaryFlag, Int_t decayFlag, Double_t time){ 

// track lifetime of primaries, DAR/DIF seperated
Time_ofLife_tillDecayed [decayFlag][primaryFlag] -> Fill(time); 

// track lifetime of primaries, of any type, DAR + DIF 
Time_ofLife_tillAnyDecay[primaryFlag] -> Fill(time); 

}

//NEW ====================================================================================================================
void SNSRootManager::FillTime(Int_t primaryFlag, Int_t decayFlag, Int_t gFlag, Double_t time){ 

// time when particles DAR/DIF separated, relative to the p+Hg Interaction, GlobalTime
Time_elapsedSince_ProtonHg_Interaction [decayFlag][primaryFlag] -> Fill(time); 

// time when particles deacy, DAR + DIF, relative to the p+Hg Interaction, GlobalTime
Time_elapsedSince_ProtonHg_Interaction_AnyDecay [primaryFlag] -> Fill(time); 

}

//NEW ====================================================================================================================
void SNSRootManager::FillEnergy(Int_t primaryFlag, Int_t decayFlag, Int_t gFlag, Double_t energy){ 

// 1D histogram of the energy when particles DAR/DIF seperated
Energy_Primaries_DAR_DIF [decayFlag][primaryFlag] -> Fill(energy); 

// 1D histogram of energy for all particles that Decay (OF ANY TYPE, DIF + DAR)
Energy_Primaries_AnyDecay [primaryFlag] -> Fill(energy); 

}

void SNSRootManager::FillEnergy(Int_t nuFlag, Int_t volFlag, Double_t energy, Int_t eFlag){
  // tracking energy of neutrinos that made it to a detector
  Energy_Nu[volFlag][nuFlag] -> Fill(energy); 
}

void SNSRootManager::FillTime(Int_t nuFlag, Int_t volFlag, Double_t time, Int_t tFlag){
  // tracking the time it took for neutrinos to make it to a detector (not relative to p+Hg currently)
  Time_Nu[volFlag][nuFlag] -> Fill(time); 
}





// Bari's methods for filling and outputing ROOT data


// saving/writing the ROOT output file
// MODIFY THIS SECTION TO CHANGE NAME NUMBER IF OUTPUT ALREADY EXISTS
void SNSRootManager:: Save()
{
    
  TFile *f = new TFile(fOutFileName.c_str(),"recreate","SNS simulation output");
  if (!f || f->IsZombie()) {
    std::ostringstream msg_err;
    msg_err << "Can't (re)create file='" << fOutFileName << "'\n";
    throw std::runtime_error(msg_err.str().c_str());
  }

  if (fHlist) fHlist->Write();

  f->Write(); f->Close();
  delete f; delete fHlist;
  cout << "\nData has been saved in file='" << fOutFileName << "'\n\n";
}  



// fill histo function
void SNSRootManager::FillHisto(string hname, Double_t x) {
  TH1F *h = (TH1F*) fHlist->FindObject(hname.c_str());
  if (h)
    h->Fill(x);
  else
    cout << "SNSRootManager::FillHisto: Can't find TH1F histo='" << hname << "'\n";
}


// setting histogram bin content
void SNSRootManager::SetHistoBin(Int_t nbin, string hname, Double_t x) {
  TH1F *h = (TH1F*) fHlist->FindObject(hname.c_str());
  if (h)
    h->SetBinContent(nbin,x);
  else
    cout << "SNSRootManager::SetHistoBin: Can't find TH1F histo='" << hname << "'\n";
}


// getting histogram bin content
void SNSRootManager::GetHistoBin(Int_t nbin, string hname,Double_t &x) {
  TH1F *h = (TH1F*) fHlist->FindObject(hname.c_str());
  if (h)
    x = h->GetBinContent(nbin);
  else
    cout << "SNSRootManager::GetHistoBin: Can't find TH1F histo='" << hname << "'\n";
}


// fill histo functions...

void SNSRootManager::FillHisto(string hname, Double_t x, Double_t y) {
  TH2F *h = (TH2F*) fHlist->FindObject(hname.c_str());
  if (h)
    h->Fill(x,y);
  else
    cout << "SNSRootManager::FillHisto: Can't find TH2F histo='" << hname << "'\n";
}

void SNSRootManager::FillHisto(string hname, Double_t x, Double_t y, Double_t z) {
  TH3F *h = (TH3F*) fHlist->FindObject(hname.c_str());
  if (h)
    h->Fill(x,y,z);
  else
    cout << "SNSRootManager::FillHisto: Can't find TH3F histo='" << hname << "'\n";
}

void SNSRootManager::FillHisto(string hname, Int_t x) {
  TH1I *h = (TH1I*) fHlist->FindObject(hname.c_str());
  if (h)
    h->Fill(x);
  else
    cout << "SNSRootManager::FillHisto: Can't find TH1I histo='" << hname << "'\n";
}
  
void SNSRootManager::FillHisto(string hname, Int_t x, Int_t y) {
  TH2I *h = (TH2I*) fHlist->FindObject(hname.c_str());
  if (h)
    h->Fill(x,y);
  else
    cout << "SNSRootManager::FillHisto: Can't find TH2I histo='" << hname << "'\n";
}

void SNSRootManager::FillHisto(string hname, Int_t x, Int_t y, Int_t z) {
  TH3I *h = (TH3I*) fHlist->FindObject(hname.c_str());
  if (h)
    h->Fill(x,y,z);
  else
    cout << "SNSRootManager::FillHisto: Can't find TH3I histo='" << hname << "'\n";
}

void SNSRootManager::FillHisto(string hname, Float_t x) {
  SNSRootManager::FillHisto(hname, (Double_t) x);
}
void SNSRootManager::FillHisto(string hname, Float_t x, Float_t y) {  
  SNSRootManager::FillHisto(hname, (Double_t) x, (Double_t) y);
}
void SNSRootManager::FillHisto(string hname, Float_t x, Float_t y, Float_t z) {  
  SNSRootManager::FillHisto(hname, (Double_t) x, (Double_t) y, (Double_t) z);
}


