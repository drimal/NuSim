#include "PrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "globals.hh"
#include "G4GeneralParticleSource.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

//#include "math.h"

PrimaryGeneratorAction::PrimaryGeneratorAction(){

   gun = new G4GeneralParticleSource();
}

PrimaryGeneratorAction::~PrimaryGeneratorAction(){

    delete gun;
}

#include "G4Event.hh"
#include "Randomize.hh"
#include "CLHEP/Random/RandomEngine.h"

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent){


  G4bool verbose = false;
 //G4bool verbose = true;

  G4double PulseDuration  = 695*ns;
  G4double PulseFrequency = 60*hertz;
  G4long    NPulses       = 10;                   // Number of pulses to simulate

  G4double TimeInPulse=G4UniformRand()*PulseDuration;
  G4int PulseID = CLHEP::RandFlat::shootInt(NPulses);

  G4double IniTime = TimeInPulse + PulseID/PulseFrequency;

  if (verbose) 
	 G4cout << "pga: initial time =" << IniTime/second << " sec = " <<
		"time in pulse=" << TimeInPulse/ns << " ns + pulse shift=" << PulseID/PulseFrequency/second <<
		" sec; " << PulseID <<  "-th Pulse\n";

  
  //gun->SetParticleTime(IniTime);
  //gun->SetParticleTime(0.0*s);

  gun->GeneratePrimaryVertex(anEvent);
}
