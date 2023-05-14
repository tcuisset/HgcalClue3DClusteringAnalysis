#ifndef CLUE3DAlgo_h
#define CLUE3DAlgo_h

#include <cmath>

#include "CLUEAlgoParameters.h"
#include "LayerTiles.h"
#include "PointsCloud.h"
#include "Cluster.h"

#define NLAYER1 1 ///< Index of first layer


inline float square(float a) {
  return a * a;
}

/**
 * Compute the distance squared between 2 points in PointsCloud, in x-y plane
*/
inline float distance3d_squared(PointsCloud &points, int i, int j) {
  /*
  auto r = [&] (int k) -> float { // Compute r = sqrt(x^2 + y^2) (ie cylindrical coordinates)
    return std::sqrt( square(points.x[k]) + square(points.y[k]));
  };
  */
  /* This is the way it is done in CMSSW. Given our z is different than CMS z, using it is probably not a good idea (esp when z=0)
  Could probably be made to work by shifting z by the distance between detectir center and firts layer of HGCAL
  return square(points.z[i]) * square( r(i)/std::abs(points.z[i]) + r(j)/std::abs(points.z[j]) )
      + square(r(j))/square(points.z[j]) * square(std::atan2(points.y[j], points.x[j]) - std::atan2(points.y[i], points.x[i]));

  */
 /* Using transverse plane distance, probably a good approximation of the CMSSW way when considering a small subset of detector at high eta
 */
 return square(points.x[j] - points.x[i]) + square(points.y[j] - points.y[i]);
}

/**
 * Compute rho, the local energy density, for each 2D cluster (in 3D)
 * Note the way the distance is computed (ignoring layer completely in the calculation of distance) for computing rho
 * \param params CLUE3D params (used : dc, densitySiblingLayers, densityOnSameLayer, kernelDensityFactor)
*/
void calculate_density3d(std::array<LayerTilesClue3D, NLAYERS> &d_hist, PointsCloud &points, Clue3DAlgoParameters const& params) {
  // loop over all 2D clusters
  for (unsigned int i = 0; i < points.n; i++) {
    int clayer = points.layer[i]; ///< Layer of 2D cluster
    int minLayer = 0;
    int maxLayer = NLAYERS;

    //Loop over layers that are within densitySiblingLayers of current layer
    minLayer = std::max(clayer - params.densitySiblingLayers, NLAYER1);
    maxLayer = std::min(clayer + params.densitySiblingLayers, maxLayer);
    for (int currentLayer = minLayer; currentLayer <= maxLayer; currentLayer++) {
      LayerTilesClue3D &lt = d_hist[currentLayer];
      /*
      In CMSSW the search is done in the next two eta-phi bins, so here we emulate this by looking into the next 2 x-y bins
      */
      float dc_effective = 2 * LayerTilesClue3D::TilesConstants::tileSize;
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

            if (points.masked[j] == 1)
              continue; // Skip masked layer cluster

            //This computes the 2D distance between the 2D clusters in x,y (ignoring layer)
            //ie 2 2D clusters will have the same distance as long as they have same x,y wether they are on the same layer or 2 layers apart
            // *TODO* Is this intended ?
            float dist_ij = distance3d_squared(points, i, j);
            //    std::cout<<"dist_ij "<<dist_ij<<std::endl;

            if (i != j && !params.densityOnSameLayer && currentLayer == clayer) {
              // if densityOnSameLayer is false, then ignore all layer clusters on the same layer, except when considering itself
              continue;
            }
            if (dist_ij <= params.densityXYDistanceSqr) {
              // sum weights within N_{dc_effective}(i)
              points.rho[i] += (i == j ? 1.f : params.kernelDensityFactor) * points.weight[j];
            }
          }  // end of interate inside this bin
        }
      }  // end of loop over bins in search box
    }  // loop over layers for each point
  }    // end of loop over points
};

/**
 * Compute, for each 2D cluster, the distance to nearest higher (and set the ID of the nearest higher)
 * \param params CLUE3D parameters (used are : deltac, outlierDeltaFactor, densitySiblingLayers, nearestHigherOnSameLayer) 
*/
void calculate_distanceToHigher3d(std::array<LayerTilesClue3D, NLAYERS> &d_hist,
                                PointsCloud &points, Clue3DAlgoParameters const& params) {
  // loop over all points
  for (unsigned int i = 0; i < points.n; i++) {
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
    minLayer = std::max(clayer - params.densitySiblingLayers, NLAYER1);
    maxLayer = std::min(clayer + params.densitySiblingLayers, maxLayer);
    for (int currentLayer = minLayer; currentLayer <= maxLayer; currentLayer++) {
      if (!params.nearestHigherOnSameLayer && (currentLayer == clayer))
        continue;
      /*
      The way CMSSW CLUE3D works i by considering on eta bin and one phi bin on either side of the current bin
      and consider *all* layer clusters in these bins (there is no filtering on distance)
      So here we emulate it by using x-y bins of size approximately equal to the size of a phi-eta bin in CMSSW
      A computation in the front face of HGCAL around the middle of the radius gave 2.5 cm in eta-equivalent dimension,
      and 5.4 cm in phi-equivalent dimension.
      Therefore we approximate this by 3.7 cm * 3.7 cm square bins (to get approx the same area)

      We still use searchBox for convenience as it handles edge cases properly
      */
      int xyBinWindow = 1;
      float dm = xyBinWindow * LayerTilesClue3D::TilesConstants::tileSize;
      LayerTilesClue3D &lt = d_hist[currentLayer];
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

            if (points.masked[j] == 1)
              continue; //Skip masked layer cluster

            // query N'_{dm}(i)
            bool foundHigher = (points.rho[j] > rho_i);
            // in the rare case where rho is the same, use detid
            // foundHigher = foundHigher || ((points.rho[j] == rho_i) && (j > i));
            float dist_ij = std::sqrt(distance3d_squared(points, i, j));

            /* Important note : CMSSW CLUE and CLUE3D does not do the check dist_ij <= dm contrary to what standalone CLUE does
            (at least as of March 2023)
            */
            if (foundHigher /* && dist_ij <= dm*/) { 
              // find the nearest point within N'_{dm}(i)
              if (dist_ij <= delta_i) {
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
 * \param params Clue3D parameters (used are deltac, outlierDeltaFactor, criticalDensity and criticalZDistanceLyr)
 * \return number of 3D clusters created
*/
int findAndAssign_clusters3d(PointsCloud &points, Clue3DAlgoParameters const& params) {
  int nClusters = 0;

  // find cluster seeds and outlier
  std::vector<int> localStack;
  // loop over all points
  for (unsigned int i = 0; i < points.n; i++) {
    float critical_transverse_distance = params.criticalXYDistance; //useAbsoluteProjectiveScale_ ? criticalXYDistance_ : criticalEtaPhiDistance_;
    // initialize clusterIndex
    points.clusterIndex[i] = -1;

    float deltai = points.delta[i];
    float rhoi = points.rho[i];
    int distanceInLayersToNearestHigher = -1;
    if (points.nearestHigher[i] >= 0) {
      // Absolute difference for unsigned int
      distanceInLayersToNearestHigher = points.layer[i] > points.layer[points.nearestHigher[i]]
        ? points.layer[i] - points.layer[points.nearestHigher[i]] : points.layer[points.nearestHigher[i]] - points.layer[i];
    }
    // determine seed or outlier
    bool isSeed = (deltai > critical_transverse_distance || distanceInLayersToNearestHigher > params.criticalZDistanceLyr) 
                  && (rhoi >= params.criticalDensity)
                  && (points.weight[i] / rhoi > params.criticalSelfDensity);
    bool isOutlier = (deltai > params.outlierDeltaFactor * critical_transverse_distance) && (rhoi < params.criticalDensity);
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
      if (points.nearestHigher[i] >= 0) //Only register if it actually has a nearest higher
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
