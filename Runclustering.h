#ifndef Runclustering_H
#define Runclustering_H

#include <cmath>
#include <iostream>
#include <vector>
#include <array>
#include <stdexcept>


#include "TBNtupleAnalyzer.h"
#include "TDirectory.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TLorentzVector.h"
#include "TProfile.h"
#include "TPRegexp.h"
#include "TObjString.h"
#include "TSystem.h"

#include "CLUEAlgo.h"
#include "CLUE3DAlgo.h"

class Runclustering : public TBNtupleAnalyzer {
 public:
  /**
   * \param listOfFilePaths list of paths to files holding trees to load
   * \param datatype can be either data or simulation, used for DWC cuts
   * \param shiftRechits use ce_clean_x_shifted and impactX_unshifted columns (suitable for data) if true, use ce_clean_x_unshifted and impactX_unshifted columns (suitable only for Monte Carlo) if false
   * \param filterDwc if true, filter events base on Delay Wire Chamber cuts
  */
  Runclustering(std::vector<std::string> listOfFilePaths,
                  const char *outFileName,
                  ClueAlgoParameters clueParams, Clue3DAlgoParameters clue3DParams,
                  std::string datatype, bool shiftRechits, bool filterDwc
                 ); 
                                              
  ~Runclustering();
  Long64_t LoadTree(Long64_t entry);
  void EventLoop(unsigned filterMinLayerClusterSize);

  std::vector<bool> *noise_flag;
  int event_count[7] = {};
  int count = 0, count_afterCuts = 0;

  std::string datatype;
  bool filterDwc;
  
  TFile* output_file_;
  TH1F *h_beamenergy;
  TH1F *h_nrechits;
  ClueAlgoParameters clueParams_;
  Clue3DAlgoParameters clue3DParams_;

  UShort_t currentNtupleNumber; ///< The ntuple nb currently loaded by the chain, inferred from the filename (-1 if could not infer)
};
#endif

#ifdef Runclustering_cxx



Runclustering::Runclustering(
    std::vector<std::string> listOfFilePaths, const char *outFileName,
    ClueAlgoParameters clueParams, Clue3DAlgoParameters clue3DParams, std::string datatype, bool shiftRechits, bool filterDwc) 
    : clueParams_(clueParams), clue3DParams_(clue3DParams), datatype(datatype), filterDwc(filterDwc), currentNtupleNumber(-1) { 
  
  TChain *tree = new TChain("relevant_branches");
  for (std::string path : listOfFilePaths) {
    tree->Add(path.c_str());
  }
  std::cout << "No. of Entries in chain  : " << tree->GetEntries() << std::endl
    << "Initiating analysis of dataset " << std::endl;

  TBNtupleAnalyzer::Init(tree, shiftRechits); //Init branch pointers for reading TTree 

  gSystem->Exec(TString("mkdir -p $(dirname \"") + outFileName +  "\")"); //Create directories as needed
  output_file_ = new TFile(outFileName, "RECREATE");
  if (output_file_->IsZombie())
    throw std::runtime_error("Could not open output file");
  //BookHistogram(outFileName, config, energy);
}


Long64_t Runclustering::LoadTree(Long64_t entry) {
  // Set the environment to read one entry
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (!fChain->InheritsFrom(TChain::Class())) return centry;
  TChain *chain = (TChain *)fChain;

  //Match at end of string underscore + numbers + .root 
  TPRegexp regexp_match_ntuple_nb("_([0-9]+)\\.root$");
  if (chain->GetTreeNumber() != fCurrent) {
    fCurrent = chain->GetTreeNumber();
    //    Notify();
    TString filename = fChain->GetCurrentFile()->GetName();

    //Get the ntuple number from file
    TObjArray *subStrL = regexp_match_ntuple_nb.MatchS(filename);
    if (subStrL->GetEntries() >= 1) {
      currentNtupleNumber = ((TObjString*) subStrL->Last())->GetString().Atoi();
    } else {
      currentNtupleNumber = -1;
    }
    delete subStrL;
  }


  return centry;
}

Runclustering::~Runclustering() {
  if (!fChain) return;
  delete fChain->GetCurrentFile();
  delete output_file_;
}

#endif

