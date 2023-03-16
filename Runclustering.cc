#define Runclustering_cxx
#include "Runclustering.h"

#include <TF1.h>
#include <math.h>

#include <cstring>
#include <string>
#include <iostream>
#include <vector>

#include "CLUEAlgo.h"
#include "CLUE3DAlgo.h"
#include "Math/SMatrix.h"
#include "Math/SVector.h"
using namespace std;


void dumpSoA(const PointsCloud & points) {
  std::cout << "SoA DUMP" << std::endl;
  for (int c = 0 ; c < points.x.size(); ++c) {
    std::cout
      << "id " << c
      << " x = " << points.x[c]
      << " y = " << points.y[c]
      << " z = " << points.z[c]
      << " layer = " << points.layer[c]
      << " energy = " << points.weight[c]
      << " rho = " << points.rho[c]
      << " delta = " << points.delta[c]
      << " isSeed = " << points.isSeed[c]
      << " clusterIdx = " << points.clusterIndex[c]
      << " nearestHigher = " << points.nearestHigher[c]
      << std::endl;
  }
}

int main(int argc, char *argv[]) 
{
  if (argc < 3) {
    cerr << "Please give at least 2 arguments : " 
         << "runList "
         << " "
         << "outputFileName" << endl
         << "and optionally 4 additional parameters for CLUE3D : "
         << "dc rhic outlierDeltaFactor densitySiblingLayers" << endl
         << "(for now only em section params + parameters for CLUE2D are kept as default)";
    return -1;
  }
  const char *inputFileList = argv[1]; //List of input files produced by TestBeamReconstruction
  const char *outFileName = argv[2];

  constexpr float MIP2GeV[3] = {0.0105, 0.0812, 0.12508};

  ClueAlgoParameters clueParameters;
  clueParameters.dc = {1.3f, 3.f * sqrt(2.f) + 0.1};
  clueParameters.rhoc = {4.f * MIP2GeV[0], 4.f * MIP2GeV[2]};
  clueParameters.outlierDeltaFactor = 2.f;

  Clue3DAlgoParameters clue3DParameters;
  if (argc >= 7) {
    clue3DParameters.dc = {std::stof(argv[3]), -1.};
    clue3DParameters.rhoc = {std::stof(argv[4]), -1.};
    clue3DParameters.outlierDeltaFactor = std::stof(argv[5]);
    clue3DParameters.densitySiblingLayers = std::stoi(argv[6]);
    cout << "Using custom CLUE3D parameters : " << clue3DParameters << endl;
  } else {
    clue3DParameters.dc = {1.3f, 3.f * sqrt(2.f) + 0.1};
    clue3DParameters.rhoc = {4.f * MIP2GeV[0], 4.f * MIP2GeV[2]};
    clue3DParameters.outlierDeltaFactor = 2.f;
    clue3DParameters.densitySiblingLayers = 2;
    cout << "Using default CLUE3D parameters : " << clue3DParameters << endl;
  }
  
  Runclustering tbCLUS(inputFileList, outFileName, clueParameters, clue3DParameters);
  tbCLUS.EventLoop();
  return 0;
}

void Runclustering::EventLoop() {

  bool NTUPLEOUT = false;
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t hgc_jentry = 0;

  cout << "nentries " << nentries << endl;
  cout << "Analyzing dataset " << endl;

  Long64_t nb = 0;
  int decade = 0;
  bool DEBUG = false;

  // counter
  int nHgc = 0, nAhc = 0, nrechits = 0;

  if (DEBUG) cout << "DEBUG: Entering event Loop" << endl;
  
  Long64_t jentry;
  
  TFile *f = output_file_;
  TTree clusters_tree("clusters", "clusters");

  // Create a unique PointsCloud object and (re)-use it to fill the output
  // ntuple.
  PointsCloud pcloud; ///< PointsCloud of all hits per event (3D)
  PointsCloud pcloud2d; ///< PointsCloud of all 2D clusters

  // Create a SoA for the output clusters
  ClustersSoA clusters_soa;
  ClustersSoA clusters3d_soa;

  float Esum_allRecHits_inGeV;
  //h_beamenergy->Fill(beamEnergy);
  //h_nrechits->Fill(NRechits);
  // Create the branches in the output ntuple.
  clusters_tree.Branch("beamEnergy", &beamEnergy);
  clusters_tree.Branch("ntupleNumber", &currentNtupleNumber);
  clusters_tree.Branch("NRechits", &NRechits);
  clusters_tree.Branch("impactX", &impactX_shifted); //These are vector<float> of size 40 (nb of layers)
  clusters_tree.Branch("impactY", &impactY_shifted);

  clusters_tree.Branch("rechits_x", &pcloud.x);
  clusters_tree.Branch("rechits_y", &pcloud.y);
  clusters_tree.Branch("rechits_z", &pcloud.z);
  clusters_tree.Branch("rechits_energy", &pcloud.weight);
  clusters_tree.Branch("rechits_layer", &pcloud.layer);
  clusters_tree.Branch("rechits_rho", &pcloud.rho);
  clusters_tree.Branch("rechits_delta", &pcloud.delta);
  clusters_tree.Branch("rechits_isSeed", &pcloud.isSeed);
  clusters_tree.Branch("rechits_isOutlier", &pcloud.isOutlier);

  clusters_tree.Branch("clus2D_x", &clusters_soa.x);
  clusters_tree.Branch("clus2D_y", &clusters_soa.y);
  clusters_tree.Branch("clus2D_z", &clusters_soa.z);
  clusters_tree.Branch("clus2D_energy", &clusters_soa.energy);
  clusters_tree.Branch("clus2D_layer", &clusters_soa.layer);
  clusters_tree.Branch("clus2D_size", &clusters_soa.size);
  clusters_tree.Branch("clus2D_idxs", &clusters_soa.hitidxs);
  clusters_tree.Branch("clus2D_rho", &pcloud2d.rho);
  clusters_tree.Branch("clus2D_delta", &pcloud2d.delta);
  clusters_tree.Branch("clus2D_isSeed", &pcloud2d.isSeed);
  clusters_tree.Branch("clus2D_isOutlier", &pcloud2d.isOutlier);

  clusters_tree.Branch("clus3D_x", &clusters3d_soa.x);
  clusters_tree.Branch("clus3D_y", &clusters3d_soa.y);
  clusters_tree.Branch("clus3D_z", &clusters3d_soa.z);
  clusters_tree.Branch("clus3D_energy", &clusters3d_soa.energy);
  clusters_tree.Branch("clus3D_layer", &clusters3d_soa.layer);
  clusters_tree.Branch("clus3D_size", &clusters3d_soa.size);
  clusters_tree.Branch("clus3D_idxs", &clusters3d_soa.hitidxs);

  for (jentry = 0; jentry < nentries; jentry++, hgc_jentry++) {

    // Reset PointsCloud
    pcloud.reset();
    pcloud2d.reset();
    // ==============print number of events done == == == == == == == =
    double progress = 10.0 * jentry / (1.0 * nentries);
    int k = int(progress);
    if (k > decade) cout << 10 * k << " %" << endl;
    decade = k;

    // ===============read this entry == == == == == == == == == == ==

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) {
      break;
      cout << "Breaking" << endl;
    }
    nb = fChain->GetEntry(jentry);

    // rescaling energy from MeV to GeV
    vector<float> scaled_energy;
    scaled_energy = *ce_clean_energy_MeV;
    double en_scale{0.001};
    std::transform(scaled_energy.begin(), scaled_energy.end(), scaled_energy.begin(), [&en_scale](auto& c){return c*en_scale;});



    // Compute clusters using CLUE
    //NLAYERS defined in LayersTilesConstants
    std::array<LayerTiles, NLAYERS> tiles; ///< Array of LayerTiles (2D layer of tiles) for each layer
    //make pointcloud
    pcloud.x = *ce_clean_x;
    pcloud.y = *ce_clean_y;
    pcloud.z = *ce_clean_z;
    pcloud.weight = scaled_energy;
    //    updateLayersAndEnergies(pcloud, comb_z_boundaries, MIP2GeV);
    pcloud.layer = *ce_clean_layer;
    pcloud.resizeOutputContainers(ce_clean_x->size());
    // Fill "TILES"
    compute_histogram(tiles, pcloud);
    // Calculate density quantities for points
    calculate_density(tiles, pcloud, clueParams_.dc);
    // Calculate nearest higher density  point
    calculate_distanceToHigher(tiles, pcloud, clueParams_.outlierDeltaFactor, clueParams_.dc);
    // get seeds and followers
    auto total_clusters = findAndAssign_clusters(pcloud, clueParams_.outlierDeltaFactor, clueParams_.dc, clueParams_.rhoc);
    std::vector<Cluster> clusters = getClusters(total_clusters, pcloud);
    // Fill in the clusters_SoA
    clusters_soa.load(clusters);

    // Compute clusters using CLUE3D [similar sequence as the 2D part, starting pointcloud here composed of the 2D clusters just made above]
    std::array<LayerTiles, NLAYERS> tiles2d;
    pcloud2d.x = clusters_soa.x ;
    pcloud2d.y = clusters_soa.y ;
    pcloud2d.z = clusters_soa.z ;
    vector<int> clusters_soao ; ///< Layer of each 2D cluster, indexed by cluster ID
    clusters_soao = clusters_soa.layer;
    std::vector<unsigned int> clusters_soau(std::begin(clusters_soao), std::end(clusters_soao));
    pcloud2d.layer = clusters_soau ;
    pcloud2d.weight = clusters_soa.energy ;
    pcloud2d.resizeOutputContainers(clusters_soa.x.size());

    compute_histogram(tiles2d, pcloud2d);



    calculate_density3d(tiles2d, pcloud2d, clue3DParams_.dc, clue3DParams_.densitySiblingLayers);
    calculate_distanceToHigher3d(tiles2d, pcloud2d, clue3DParams_.outlierDeltaFactor, clue3DParams_.dc, clue3DParams_.densitySiblingLayers);
    auto total_clusters3d = findAndAssign_clusters3d(pcloud2d, clue3DParams_.outlierDeltaFactor, clue3DParams_.dc, clue3DParams_.rhoc);
    auto clusters3d = getClusters3d(total_clusters3d, pcloud2d);
    clusters3d_soa.load(clusters3d);
    
    clusters_tree.Fill();
    if (DEBUG) {
      dumpSoA(pcloud);
    }

    float total_energy_clustered = 0.f;
    for (auto const & cl : clusters) {
      auto pos = cl.position();
      if (DEBUG) {
        std::cout << std::get<0>(pos) << " "
          << std::get<1>(pos) << " "
          << std::get<2>(pos) << " "
          << cl.energy() << " "
          << std::endl;
      }
      total_energy_clustered += cl.energy();
    }


  }  // loop over entries
  f->cd();
  clusters_tree.Write();
  cout<<"Total events:"<<hgc_jentry<<endl;
  ///////////////////////////////////////////////////////////
  ///////  E N D     O F     E N T R Y     L O O P     //////
  ///////////////////////////////////////////////////////////

  cout << "Done :" << jentry << endl;
  f->Write();
}

