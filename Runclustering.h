#ifndef Runclustering_H
#define Runclustering_H

#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include <array>


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

#include "CLUEAlgo.h"
#include "CLUE3DAlgo.h"

class Runclustering : public TBNtupleAnalyzer {
 public:
  /**
   * \param inputFileList path to a file listing all paths to the input root files
   * \param energy not used (set in EventLoop)
  */
  Runclustering(const TString &inputFileList,
                  const char *outFileName,
                  ClueAlgoParameters clueParams, Clue3DAlgoParameters clue3DParams
                 ); 
                                              
  ~Runclustering();
  Bool_t FillChain(TChain *chain, const TString &inputFileList); ///< Add all trees from the file list to the given TChain
  Long64_t LoadTree(Long64_t entry);
  void EventLoop();  //, const char *,const char *);

  std::vector<bool> *noise_flag;
  int event_count[7] = {};
  int count = 0, count_afterCuts = 0;
  
  TFile* output_file_;
  TH1F *h_beamenergy;
  TH1F *h_nrechits;
  ClueAlgoParameters clueParams_;
  Clue3DAlgoParameters clue3DParams_;

};
#endif

#ifdef Runclustering_cxx



Runclustering::Runclustering(
    const TString &inputFileList, const char *outFileName,
    ClueAlgoParameters clueParams, Clue3DAlgoParameters clue3DParams) 
    : clueParams_(clueParams_), clue3DParams_(clue3DParams) { 
  
  TChain *tree = new TChain("relevant_branches");
  if (!FillChain(tree, inputFileList)) {
    std::cerr << "Cannot get the tree " << std::endl;
  } else {
    std::cout << "Initiating analysis of dataset " << std::endl;
  }


  TBNtupleAnalyzer::Init(tree); //Init branch pointers for reading TTree 

  output_file_ = new TFile(outFileName, "RECREATE");
  //BookHistogram(outFileName, config, energy);
}


Bool_t Runclustering::FillChain(
    TChain *chain,
    const TString &inputFileList) {
                                   

  ifstream infile(inputFileList, ifstream::in);
  std::string buffer;

  if (!infile.is_open()) {
    std::cerr << "** ERROR: Can't open '" << inputFileList << "' for input"
              << std::endl;
    return kFALSE;
  }

  std::cout << "TreeUtilities : FillChain " << std::endl;
  while (1) {
    infile >> buffer;
    if (!infile.good()) break;

    chain->Add(buffer.c_str());

  }
  std::cout << "No. of Entries in chain  : " << chain->GetEntries()
            << std::endl;


  return kTRUE;
}

Long64_t Runclustering::LoadTree(Long64_t entry) {
  // Set the environment to read one entry
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (!fChain->InheritsFrom(TChain::Class())) return centry;
  TChain *chain = (TChain *)fChain;
  if (chain->GetTreeNumber() != fCurrent) {
    fCurrent = chain->GetTreeNumber();
    //    Notify();
  }


  return centry;
}

Runclustering::~Runclustering() {
  if (!fChain) return;
  delete fChain->GetCurrentFile();
  delete output_file_;
}

#endif

