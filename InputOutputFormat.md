# Format of input ROOT tree for ClusteringAnalysis
TTree name : relevant_branches
Columns :
* event / i
* run /i 
* NRecHits /i
* ce_clean_detid : vector<uint>
* ce_clean_x : vector<float> -> Point position
* ce_clean_y : vector<float>
* ce_clean_z : vector<float>
* ce_clean_layer : vector<uint> -> Layer number of hit

branch: ce_clean_energy_MeV   44126381

branch: beamEnergy                 600

## impact_shifted
For each event : vector<float> of length 40 (nb of layers)
branch: impactX_shifted       
branch: impactY_shifted       



* Instance
* event.eve

* NRechits :branch entirely identical for whole event with nb of reconstructed hits

# Format of output
## 2nd output root file : CLUE_clusters.root
one TTree name clusters
beamEnergy : float
NRechits : uint
impactX, impactY : vector<float> of length 40 (nb of layers)

rechits_* -> PointsCloud (vector branch) :
  std::vector<float> x;
  std::vector<float> y;
  std::vector<float> z;
  std::vector<unsigned int> layer; 
  std::vector<float> weight; ///< Weights of all points (ie energy) -> Called energy in tree

  // Output variables
  std::vector<float> rho; ///< Local energy density
  std::vector<float> delta; ///< Distance to nearest highest

clus2D_* -> ClustersSoA (vector branch) : list of all 2D clusters info
  std::vector<float> x;
  std::vector<float> y;
  std::vector<float> z;
  std::vector<float> energy;
  std::vector<int> layer;
  std::vector<int> size; ///< Number of hits per cluster
  std::vector<std::vector<int>> hitidxs; ///< List of hits IDs per cluster
  also rho, delta
  
clus3D_* -> ClustersSoA (vector branch) : same but 3D clusters
