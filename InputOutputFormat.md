# Format of input ROOT tree for ClusteringAnalysis
see details : <https://indico.cern.ch/event/770743/contributions/3204606/attachments/1747143/2829163/Reco_ntuples_06Nov2018.pdf> and <https://gitlab.cern.ch/cms-hgcal-tb/TestBeam/-/wikis/samples/ntuples-description>
TTree name : relevant_branches
Columns :
* event / i
* run /i 
* NRecHits /i : branch entirely identical for whole event with nb of reconstructed hits
* ce_clean_detid : vector<uint>
* ce_clean_x : vector<float> -> Point position
* ce_clean_y : vector<float>
* ce_clean_z : vector<float>
* ce_clean_layer : vector<uint> -> Layer number of hit
* ce_clean_energy : vector<float> -> Energy in MIP
* ce_clean_energy_MeV : vector<float> -> Energy in MeV
* beamEnergy /F : nominal beam energy (20, 50, ... 300 GeV)

For simulation only :
* trueBeamEnergy : true particle momentum,


## impact_shifted
For each event : vector<float> of length 40 (nb of layers)
* impactX_shifted  -> to be compared against ce_clean_x (shifts have been applied to this purpose)
* impactY_shifted       

For data only : 
* ce_clean_x_shifted -> shifted hits positions to take into account misalignment of different layers of HGCAL prototype  
* impactX_unshifted -> impact values, without any shift, to compare against ce_clean_x_shifted

## DWC information (see <https://gitlab.cern.ch/cms-hgcal-tb/TestBeam/-/wikis/samples/ntuples-description>)
 * DWC_b_x and DWC_b_y : track offsets = impact onto EE to make cuts on events where beam is too far of center of DWC and detector (these are taken without any correction from the original ntuples, so they might be mirrored in data)
 * DWC_trackChi2_X and DWC_trackChi2_Y : chi2 of extrapolated tracks, straight line model

# Format of output
## 2nd output root file : CLUE_clusters.root
one TTree name clusters
beamEnergy : float
NRechits : uint
impactX, impactY : vector<float> of length 40 (nb of layers). Cuurently it is always unshifted
DWC_passes_cuts : bool : whether the event passes DWC cuts (always true if --filter-dwc was passed to runclustering, since events that fail cuts were discarded from dataset)

rechits_* -> PointsCloud (vector branch) :
  std::vector<float> x;  // Whether it is shifted or not depends on --shift-rechits command line parameter
  std::vector<float> y;
  std::vector<float> z;
  std::vector<unsigned int> layer; 
  std::vector<float> weight; ///< Weights of all points (ie energy) -> Called energy in tree
  
  std::vector<float> rechits_energy_MIP; ///< Rechits energy in MIPs (straight copy from ce_clean_energy)

  // Output variables
  std::vector<float> rho; ///< Local energy density
  std::vector<float> delta; ///< Distance to nearest highest
  std::vector<char> pointType; /// 0 is Follower, 1 is Seed, 2 is Outlier

clus2D_* -> ClustersSoA (vector branch) : list of all 2D clusters info
  std::vector<float> x;
  std::vector<float> y;
  std::vector<float> z;
  std::vector<float> energy;
  std::vector<int> layer;
  std::vector<int> size; ///< Number of hits per cluster
  std::vector<std::vector<int>> hitidxs; ///< List of hits IDs per cluster
  also rho, delta, pointType
  
clus3D_* -> ClustersSoA (vector branch) : same but 3D clusters
