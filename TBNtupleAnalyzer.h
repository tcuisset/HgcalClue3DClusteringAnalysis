///////////////////////////o///////////////////////////////
// This class has been automatically generated on
// Tue Feb 14 12:16:39 2023 by ROOT version 6.18/04
// from TTree relevant_branches/relevant_branches
// found on file: ntuple_selection_data_em_415.root
//////////////////////////////////////////////////////////

#ifndef TBNtupleAnalyzer_h
#define TBNtupleAnalyzer_h


#include <TChain.h>
#include <TFile.h>
#include <TROOT.h>
// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
using namespace std;

class TBNtupleAnalyzer {
public :
   TBNtupleAnalyzer(TTree *  = 0) : fChain(0) {}
  ~TBNtupleAnalyzer(){}
   // Int_t    Cut(Long64_t entry);
   Int_t GetEntry(Long64_t entry, Int_t getall = 0) {
     return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0;
   }
   //Long64_t LoadTree(Long64_t entry);
   void     Init(TTree *tree, bool shiftRechits);
   //void     Loop();
   Bool_t   Notify();
   //void     Show(Long64_t entry = -1);

   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          event;
   UInt_t          run;
   UInt_t          NRechits;
   vector<unsigned int> *ce_clean_detid;
   vector<float>   *ce_clean_x;
   vector<float>   *ce_clean_y;
   vector<float>   *ce_clean_z;
   vector<unsigned int> *ce_clean_layer;
   vector<float>   *ce_clean_energy_MIP; ///< Rechit energies in MIP
   vector<float>   *ce_clean_energy_MeV; ///< Rechit energies in MeV
   Float_t         beamEnergy;
   Float_t         trueBeamEnergy; // only filled if present in input dataset, ie for simulation
   vector<float>   *impactX;
   vector<float>   *impactY;

   Double_t DWC_b_x; ///< from delay wire chambers :  track offsets = impact onto EE
   Double_t DWC_b_y;
   Float_t DWC_trackChi2_X; ///< from delay wire chambers : chi2 of extrapolated tracks, straight line model
   Float_t DWC_trackChi2_Y;

   //   TBNtupleAnalyzer(TTree * tree = 0);
//   Tree *fChain;  //! pointer to the analyzed TTree or TChain
   // TTree          *fChain2;   //!pointer to the analyzed TTree or TChain
   // TTree          *fChain3;   //!pointer to the analyzed TTree or TChain
//   Int_t fCurrent;  //! current Tree number in a TChain  
   

};

#endif

#ifdef TBNtupleAnalyzer_cxx


/**
 * \param shiftRechits : if true, use ce_clean_x_shifted and impactX_unshifted columns (suitable for data).
 *          if false, use ce_clean_x_unshifted and impactX_unshifted columns (suitable only for Monte Carlo)
*/
void TBNtupleAnalyzer::Init(TTree *tree, bool shiftRechits)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   ce_clean_detid = 0;
   ce_clean_x = 0;
   ce_clean_y = 0;
   ce_clean_z = 0;
   ce_clean_layer = 0;
   ce_clean_energy_MeV = 0;
   impactX = 0;
   impactY = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("event", &event);
   fChain->SetBranchAddress("run", &run);
   fChain->SetBranchAddress("NRechits", &NRechits);
   fChain->SetBranchAddress("ce_clean_detid", &ce_clean_detid);
   TString rechitsShiftString(shiftRechits ? "_shifted" : "_unshifted");
   fChain->SetBranchAddress(TString("ce_clean_x")+rechitsShiftString, &ce_clean_x);
   fChain->SetBranchAddress(TString("ce_clean_y")+rechitsShiftString, &ce_clean_y);
   fChain->SetBranchAddress("ce_clean_z", &ce_clean_z);
   fChain->SetBranchAddress("ce_clean_layer", &ce_clean_layer);
   fChain->SetBranchAddress("ce_clean_energy", &ce_clean_energy_MIP);
   fChain->SetBranchAddress("ce_clean_energy_MeV", &ce_clean_energy_MeV);
   fChain->SetBranchAddress("beamEnergy", &beamEnergy);

   if (fChain->GetBranch("trueBeamEnergy")) {
      fChain->SetBranchAddress("trueBeamEnergy", &trueBeamEnergy);
   }

   //Delay wire chambers
   fChain->SetBranchAddress("impactX_unshifted", &impactX);
   fChain->SetBranchAddress("impactY_unshifted", &impactY);
   fChain->SetBranchAddress("DWC_b_x", &DWC_b_x);
   fChain->SetBranchAddress("DWC_b_y", &DWC_b_y);
   fChain->SetBranchAddress("DWC_trackChi2_X", &DWC_trackChi2_X);
   fChain->SetBranchAddress("DWC_trackChi2_Y", &DWC_trackChi2_Y);

   Notify();
}

Bool_t TBNtupleAnalyzer::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef TBNtupleAnalyzer_cxx
