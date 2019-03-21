#ifndef SNSRootManager_h
#define SNSRootManager_h 1

#include <TTree.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TList.h>

using namespace std;

class SNSRootManager
{
	public:
		SNSRootManager(string fout);
		virtual ~SNSRootManager();     

		// methods
		void Init();
		void Save();

		// methods for ROOT objects
		void FillHisto(string hname, Double_t x);
		void FillHisto(string hname, Double_t x, Double_t y);
		void FillHisto(string hname, Double_t x, Double_t y, Double_t z);
		void FillHisto(string hname, Int_t x);
		void FillHisto(string hname, Int_t x, Int_t y);
		void FillHisto(string hname, Int_t x, Int_t y, Int_t z);
		void FillHisto(string hname, Float_t x);
		void FillHisto(string hname, Float_t x, Float_t y);
		void FillHisto(string hname, Float_t x, Float_t y, Float_t z);
		void SetHistoBin(Int_t nbin, string hname,Double_t x);
		void GetHistoBin(Int_t nbin, string hname,Double_t &x);


		// methods for filling histogram arrays
		void FillProcess(Int_t partFlag, Int_t procFlag);
		void FillMCNPID(Int_t partFlag, Int_t mcnpID);
		void FillDestroyed(Int_t parentFlag, Int_t procFlag);
    void Fill2DEnergyTime(Int_t particleFlag, Double_t energy, Double_t time, Int_t pType);


		void FillEnergy(Int_t partFlag, Double_t energy);
		void FillEnergy(Int_t parentFlag, Int_t volFlag, Double_t energy, Int_t eFlag);
    void FillEnergy(Int_t primaryFlag, Int_t decayFlag, Int_t gFlag, Double_t energy);

    void FillTime(Int_t particleFlag, Double_t time, Int_t timeFlag);
    void FillTime(Int_t primaryFlag, Int_t decayFlag, Double_t time);
		void FillTime(Int_t parentFlag, Int_t volFlag, Double_t time, Int_t tFlag);
    void FillTime(Int_t primaryFlag, Int_t decayFlag, Int_t gFlag, Double_t time);

    void FillEnergy_1STEP(Int_t primaryFlag, Double_t energy);
    void FillEnergy_KILLSTEP(Int_t primaryFlag, Double_t energy);

    void Fill3D(Double_t x, Double_t y, Double_t z);
    void Fill2D(Double_t x, Double_t y);
    void Fill1D(Double_t x);

		// a string array containing all the secondary particle names that we are interested in, used for the output file
		char *partNames[13] = {
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


		char *parentNames[33] = {
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

		// a string array containing all the primary particle names that we are interested in tracking
		char *primaryNames[10] = {
			"Proton",
			"Neutron",
			"Pion+",
			"Pion-",
			"Muon+",
			"Muon-",
      "NuMu",
      "Anti-NuMu",
      "NuE",
      "Anti-NuE",
		};

    char *volNames[4] = {
    "oscsns",
    "near",
    "forward",
    "all_space"
};

    char *nuNames[4] = {
    "nu_mu",
    "anti_nu_mu",
    "nu_e",
    "anti_nu_e"
};

  // a string array containing the decay type names
    char *decayNames[3] = {
    "DAR",
    "DIF",
    "all"
};

	private:
		// data members
		static  SNSRootManager* fgInstance; //Singleton instance

		// data members
		TH1F                 *fHist;
		TList                *fHlist;
		string               fOutFileName;

	public:
		static SNSRootManager* Instance() {return fgInstance;}

};

#endif




