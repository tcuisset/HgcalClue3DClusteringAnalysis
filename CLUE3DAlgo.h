#ifndef CLUE3DAlgo_h
#define CLUE3DAlgo_h

#include <cmath>

#include "LayerTiles.h"
#include "PointsCloud.h"
#include "Cluster.h"

#define NLAYER1 1 ///< Index of first layer

struct Clue3DAlgoParameters
{
    std::array<float, 2> dc;
    std::array<float, 2> rhoc;
    float outlierDeltaFactor;
    int densitySiblingLayers; ///< define range of layers +- layer# of a point  

    friend std::ostream& operator<< (std::ostream& stream, const Clue3DAlgoParameters& p) {
        stream << "dc = " << p.dc[0] 
               << ", rhoc = " << p.rhoc[0]
               << ", outlierDeltaFactor = " << p.outlierDeltaFactor 
               << ", densitySiblingLayers = " << p.densitySiblingLayers;
        return stream;
    }
};

/**
 * Compute the distance between 2 points in PointsCloud in *2D* using (x;y) coordinates only (ignoring z completely)
 * *TODO* give a better name to this function
*/
inline float distance3d(PointsCloud &points, int i, int j) {
    const float dx = points.x[i] - points.x[j];
    const float dy = points.y[i] - points.y[j];
    return std::sqrt(dx * dx + dy * dy);
}

/**
 * Compute rho, the local energy density, for each 2D cluster (in 3D)
 * Note the way the distance is computed (ignoring layer completely in the calculation of distance) for computing rho
 * \param dc distance parameters (array as depends on layer)
 * \param densitySiblingLayers How many layers to consider before and beyond each cluster to compute rho
*/
void calculate_density3d(std::array<LayerTiles, NLAYERS> &d_hist,
			 PointsCloud &points, std::array<float, 2> dc, int densitySiblingLayers) {
  // loop over all 2D clusters
  for (unsigned int i = 0; i < points.n; i++) {
    int clayer = points.layer[i]; ///< Layer of 2D cluster
    int minLayer = 0;
    int maxLayer = NLAYERS;

    //Loop over layers that are within densitySiblingLayers of current layer
    minLayer = std::max(clayer - densitySiblingLayers, NLAYER1);
    maxLayer = std::min(clayer + densitySiblingLayers, maxLayer);
    for (int currentLayer = minLayer; currentLayer <= maxLayer; currentLayer++) {
      LayerTiles &lt = d_hist[currentLayer];
      float dc_effective = currentLayer < 41 ? dc[0] : dc[1];
      // get search box (2D)
      std::array<int, 4> search_box = lt.searchBox(
	     points.x[i] - dc_effective, points.x[i] + dc_effective,
	     points.y[i] - dc_effective, points.y[i] + dc_effective);
      
      // loop over bins in the search box
      for (int xBin = search_box[0]; xBin < search_box[1] + 1; ++xBin) {
        for (int yBin = search_box[2]; yBin < search_box[3] + 1; ++yBin) {
          
          // get the id of this bin
          int binId = lt.getGlobalBinByBin(xBin, yBin);
          // get the size of this bin
          int binSize = lt[binId].size();
          
          // iterate inside this bin
          for (int binIter = 0; binIter < binSize; binIter++) {
            unsigned int j = lt[binId][binIter];
            // query N_{dc_effective}(i)

            //This computes the 2D distance between the 2D clusters in x,y (ignoring layer)
            //ie 2 2D clusters will have the same distance as long as they have same x,y wether they are on the same layer or 2 layers apart
            // *TODO* Is this intended ?
            float dist_ij = distance3d(points, i, j);
            //    std::cout<<"dist_ij "<<dist_ij<<std::endl;
            
            if (dist_ij <= dc_effective) {
              // sum weights within N_{dc_effective}(i)
              points.rho[i] += (i == j ? 1.f : 0.5f) * points.weight[j];
            }
          }  // end of interate inside this bin
        }
      }  // end of loop over bins in search box
    }  // loop over layers for each point
  }    // end of loop over points
};

/**
 * Compute, for each 2D cluster, the distance to nearest higher (and set the ID of the nearest higher)
 * \param dc distance parameters (array as depends on layer)
 * \param outlierDeltaFactor multiplicative factor to dc to get distance to search for nearest higher
 * \param densitySiblingLayers How many layers to consider before and beyond each cluster to search for nearest higher
*/
void calculate_distanceToHigher3d(std::array<LayerTiles, NLAYERS> &d_hist,
                                PointsCloud &points, float outlierDeltaFactor,
				  std::array<float, 2> dc, int densitySiblingLayers) {
  // loop over all points
  for (unsigned int i = 0; i < points.n; i++) {
    float dc_effective = points.layer[i] < 41 ? dc[0] : dc[1];
    float dm = outlierDeltaFactor * dc_effective;
    // default values of delta and nearest higher for i
    float delta_i = std::numeric_limits<float>::max();
    int nearestHigher_i = -1;
    float xi = points.x[i];
    float yi = points.y[i];
    float rho_i = points.rho[i];

    int minLayer = 0;
    int maxLayer = NLAYERS;
    int clayer = points.layer[i];

    //Loop over layers that are within densitySiblingLayers of current layer
    minLayer = std::max(clayer - densitySiblingLayers, NLAYER1);
    maxLayer = std::min(clayer + densitySiblingLayers, maxLayer);
    for (int currentLayer = minLayer; currentLayer <= maxLayer; currentLayer++) {
      
      
      // get search box (2D)
      //      LayerTiles &lt = d_hist[points.layer[currentLayer]];
      LayerTiles &lt = d_hist[currentLayer];
      std::array<int, 4> search_box =
        lt.searchBox(xi - dm, xi + dm, yi - dm, yi + dm);
      
      // loop over all bins in the search box
      for (int xBin = search_box[0]; xBin < search_box[1] + 1; ++xBin) {
        for (int yBin = search_box[2]; yBin < search_box[3] + 1; ++yBin) {
          // get the id of this bin
          int binId = lt.getGlobalBinByBin(xBin, yBin);
          // get the size of this bin
          int binSize = lt[binId].size();
          
          // interate inside this bin
          for (int binIter = 0; binIter < binSize; binIter++) {
            unsigned int j = lt[binId][binIter];
            // query N'_{dm}(i)
            bool foundHigher = (points.rho[j] > rho_i);
            // in the rare case where rho is the same, use detid
            // foundHigher = foundHigher || ((points.rho[j] == rho_i) && (j > i));
            float dist_ij = distance3d(points, i, j);
            if (foundHigher && dist_ij <= dm) {  // definition of N'_{dm}(i)
              // find the nearest point within N'_{dm}(i)
              if (dist_ij < delta_i) {
                // update delta_i and nearestHigher_i
                delta_i = dist_ij;
                nearestHigher_i = j;
              }
            }
          }  // end of interate inside this bin
        }
      }  // end of loop over bins in search box
    } // end of loop over layers
    points.delta[i] = delta_i;
    points.nearestHigher[i] = nearestHigher_i;
  }  // end of loop over points
};

/**
 * For all 2D clusters, compute whether it is a seed, outlier, or follower (in this case register to the nearest higher)
 * Then expand clusters from seeds (setting point.clusterIndex for all points in each cluster)
 * \param dc distance parameters (array as depends on layer)
 * \param outlierDeltaFactor multiplicative factor to dc to get distance to search for nearest higher
 * \param rhoc critical energy density parameters (array as depends on layer)
 * \return number of 3D clusters created
*/
int findAndAssign_clusters3d(PointsCloud &points, float outlierDeltaFactor,
                           std::array<float, 2> dc,
                           std::array<float, 2> rhoc) {
  int nClusters = 0;

  // find cluster seeds and outlier
  std::vector<int> localStack;
  // loop over all points
  for (unsigned int i = 0; i < points.n; i++) {
    float dc_effective = points.layer[i] < 41 ? dc[0] : dc[1];
    float rhoc_effective = points.layer[i] < 41 ? rhoc[0] : rhoc[1];
    // initialize clusterIndex
    points.clusterIndex[i] = -1;

    float deltai = points.delta[i];
    float rhoi = points.rho[i];
    // determine seed or outlier
    bool isSeed = (deltai > dc_effective) && (rhoi >= rhoc_effective);
    bool isOutlier = (deltai > outlierDeltaFactor * dc_effective) && (rhoi < rhoc_effective);
    if (isSeed) {
      //std::cout<<"found seed"<<std::endl;  
    // set isSeed as 1
      points.pointType[i] = PointsCloud::SEED;
      // set cluster id
      points.clusterIndex[i] = nClusters;
      // increment number of clusters
      nClusters++;
      // add seed into local stack
      localStack.push_back(i);
    } else if (!isOutlier) {
      // register as follower at its nearest higher
      points.followers[points.nearestHigher[i]].push_back(i);
      points.pointType[i] = PointsCloud::FOLLOWER;
    } else {
      points.pointType[i] = PointsCloud::OUTLIER;
    }
  }

  // expend clusters from seeds
  while (!localStack.empty()) {
    int i = localStack.back();
    auto &followers = points.followers[i];
    localStack.pop_back();

    // loop over followers
    for (int j : followers) {
      // pass id from i to a i's follower
      points.clusterIndex[j] = points.clusterIndex[i];
      // push this follower to localStack
      localStack.push_back(j);
    }
  }
  //  std::cout << "Total clusters created: " <<  nClusters << std::endl;
  return nClusters;
};

/**
 * Compute position and total weight of a given cluster
 * Position is computed by weighted average of points positions (for x and y and z).
 * Weights are computed using an affine log scale of hits weights (ie hit energy).
*/
void calculatePosition3d(const PointsCloud & points, Cluster & cl) {
  constexpr float thresholdW0 = 2.9f;
  float total_weight = 0.f;
  float x = 0.f;
  float y = 0.f;
  float z = 0.f;

  unsigned int maxEnergyIndex = 0;
  float maxEnergyValue = 0.f;

  // loop over hits in cluster candidate
  // determining the maximum energy hit
  for (auto i : cl.hits()) {
    total_weight += points.weight[i];
    if (points.weight[i] > maxEnergyValue) {
      maxEnergyValue = points.weight[i];
      maxEnergyIndex = i;
    }
  }

  // TODO: this is recomputing everything twice and overwriting the position with log weighting position
  float total_weight_log = 0.f;
  float x_log = 0.f;
  float y_log = 0.f;
  float z_log = 0.f;
  for (auto i : cl.hits()) {
//    //for silicon only just use 1+6 cells = 1.3cm for all thicknesses
//    if (distance2(i, maxEnergyIndex, layerId, false) > positionDeltaRho2_)
//      continue;
    float rhEnergy = points.weight[i];
    float Wi = std::fmax(thresholdW0 + std::log(rhEnergy / total_weight), 0.);
    x_log += points.x[i] * Wi;
    y_log += points.y[i] * Wi;
    z_log += points.z[i] * Wi;
    total_weight_log += Wi;
  }

  total_weight = total_weight_log;
  x = x_log;
  y = y_log;
  z = z_log;
  if (total_weight != 0.) {
    float inv_tot_weight = 1.f / total_weight;
    cl.setPosition(x * inv_tot_weight,
                   y * inv_tot_weight,
                   z * inv_tot_weight);
  } else {
    return cl.setPosition(0.f, 0.f, 0.f);
  }
}

/**
 * Build all Cluster objects holding 3D cluster information
 * Cluster layer is undefined (since it is a 3D cluster)
 * \param totalClusters the nb of clusters generated by the algorithm (return value of \a findAndAssign_clusters3d )
 * \param points the PointsCloud ( *TODO* why is this passed by value ?)
*/
std::vector<Cluster> getClusters3d(int totalClusters, PointsCloud points) {
  std::vector<Cluster> clusters(totalClusters);

  // loop over all points
  for (unsigned int i = 0; i < points.n; i++) {
    // If it has a valid cluster index, assign that point to the corresponding cluster
    if (points.clusterIndex[i] != -1) {
      auto & thisCluster = clusters[points.clusterIndex[i]];
      thisCluster.addCell(i);
//      thisCluster.addEnergyAndRescale(points.weight[i], points.z[i]);
      thisCluster.addEnergy(points.weight[i]);
      //      thisCluster.setLayer(points.layer[i]);
      thisCluster.setLayer(-99);
    }
  }

  // Now that we have all clusters, compute their barycenter
  for (auto & cl : clusters) {
    calculatePosition3d(points, cl);
  }

  return clusters;
}

#endif
