#ifndef CLUEAlgo_h
#define CLUEAlgo_h

#include <cmath>
#include <array>

#include "CLUEAlgoParameters.h"
#include "LayerTiles.h"
#include "PointsCloud.h"
#include "Cluster.h"


///< 2-d distance on the layer
inline float distance(PointsCloud const& points, int i, int j) {
  if (points.layer[i] == points.layer[j]) {
    const float dx = points.x[i] - points.x[j];
    const float dy = points.y[i] - points.y[j];
    return std::sqrt(dx * dx + dy * dy);
  } else {
    return std::numeric_limits<float>::max();
  }
}

//Never used
inline void updateLayersAndEnergies(PointsCloud &points,
    const float z_boundaries[],
    const float mip2gev[]) {
  points.layer.resize(points.z.size(), 0);
  int counter = 0;
  for (auto z : points.z) {
    for (unsigned int l = 0; l < 50; ++l) {
      if (z_boundaries[l] > z) {
        points.layer[counter] = l+1;
        break;
      }
    }

    // Rescale energy to GeV, according to the layer
    if (points.layer[counter] < 29) {
      points.weight[counter] *= mip2gev[0];
    } else if (points.layer[counter] < 41) {
      points.weight[counter] *= mip2gev[1];
    } else {
      points.weight[counter] *= mip2gev[2];
    }
//    std::cout << "Got Layer " << points.layer[counter] << " for Z " << points.z[counter]
//      << " with calibrated energy: " << points.weight[counter] << std::endl;
    // Next Rechits
    ++counter;
  }
}

///< For all points in PointsCloud, add the index of the point into the correct bin of LayerTiles
template<typename LayerTiles>
void compute_histogram(std::array<LayerTiles, NLAYERS> &d_hist,
                       PointsCloud &points) {
  //    std::cout<<"points.n "<<points.n<<std::endl;

  for (unsigned int i = 0; i < points.n; i++) {
    // push index of points into tiles
    d_hist[points.layer[i]].fill(points.x[i], points.y[i], i);
  }
};

/**
 * Compute rho, the local energy density, for each point of the cloud (in 2D)
 * \param dc distance parameters (array as depends on layer)
*/
void calculate_density(std::array<LayerTilesClue, NLAYERS> &d_hist,
                       PointsCloud &points, std::array<float, 2> dc) {
  // loop over all points
  for (unsigned int i = 0; i < points.n; i++) {
    LayerTilesClue &lt = d_hist[points.layer[i]];

    float dc_effective = points.layer[i] < 41 ? dc[0] : dc[1];

    // get search box
    std::array<int, 4> search_box = lt.searchBox(
        points.x[i] - dc_effective, points.x[i] + dc_effective,
        points.y[i] - dc_effective, points.y[i] + dc_effective);

    // loop over bins in the search box
    for (int xBin = search_box[0]; xBin < search_box[1] + 1; ++xBin) {
      for (int yBin = search_box[2]; yBin < search_box[3] + 1; ++yBin) {
        // get the id of this bin
        int binId = lt.getGlobalBinByBin(xBin, yBin);
        // get the number of hits in this bin ("bin size")
        int binSize = lt[binId].size();

        // iterate inside this bin
        for (int binIter = 0; binIter < binSize; binIter++) {
          unsigned int j = lt[binId][binIter]; //Hit index
          // query N_{dc_effective}(i)
          float dist_ij = distance(points, i, j); //Distance between initial hit (i) and searched hit (j)
          if (dist_ij <= dc_effective) {
            // sum weights within N_{dc_effective}(i)
            points.rho[i] += (i == j ? 1.f : 0.5f) * points.weight[j];
          }
        }  // end of interate inside this bin
      }
    }  // end of loop over bins in search box
  }    // end of loop over points
};

/**
 * Compute, for each hit, the distance to nearest higher (and set the ID of the nearest higher)
 * \param dc distance parameters (array as depends on layer)
 * \param outlierDeltaFactor multiplicative factor to dc to get distance to search for nearest higher
*/
void calculate_distanceToHigher(std::array<LayerTilesClue, NLAYERS> &d_hist,
                                PointsCloud &points, float outlierDeltaFactor,
                                std::array<float, 2> dc) {
  // loop over all points
  for (unsigned int i = 0; i < points.n; i++) {
    float dc_effective = points.layer[i] < 41 ? dc[0] : dc[1];
    float dm = outlierDeltaFactor * dc_effective;
    // default values of delta and nearest higher for i
    float delta_i = std::numeric_limits<float>::max(); //Minimum value of distance yet found in the loop
    int nearestHigher_i = -1;
    float xi = points.x[i];
    float yi = points.y[i];
    float rho_i = points.rho[i];

    // get search box
    LayerTilesClue &lt = d_hist[points.layer[i]];
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
          float dist_ij = distance(points, i, j);
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

    points.delta[i] = delta_i;
    points.nearestHigher[i] = nearestHigher_i;
  }  // end of loop over points
};

/**
 * For all points, compute whether it is a seed, outlier, or follower (in this case register to the nearest higher)
 * Then expand clusters from seeds (setting point.clusterIndex for all points in each cluster)
 * \param dc distance parameters (array as depends on layer)
 * \param outlierDeltaFactor multiplicative factor to dc to get distance to search for nearest higher
 * \param rhoc critical energy density parameters (array as depends on layer)
 * \return number of clusters created
*/
int findAndAssign_clusters(PointsCloud &points, float outlierDeltaFactor,
                           std::array<float, 2> dc,
                           std::array<float, 2> rhoc) {
  int nClusters = 0;

  // find cluster seeds and outlier
  std::vector<int> localStack; ///< Stack holding points IDs, used for filling in clusterIndex. Initialized to all seeds
  // loop over all points
  for (unsigned int i = 0; i < points.n; i++) {
    float dc_effective = points.layer[i] < 41 ? dc[0] : dc[1];
    float rhoc_effective = points.layer[i] < 41 ? rhoc[0] : rhoc[1];
    // initialize clusterIndex
    points.clusterIndex[i] = -1;

    float deltai = points.delta[i];
    float rhoi = points.rho[i];

    // determine seed or outlier
    // Seed if : far enough from nearest higher && enough energy density
    bool isSeed = (deltai > dc_effective) && (rhoi >= rhoc_effective);
    // Outlier if : very far from nearest higher && low energy density
    bool isOutlier = (deltai > outlierDeltaFactor * dc_effective) && (rhoi < rhoc_effective);
    if (isSeed) {
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
  // std::cout << "Total clusters created: " <<  nClusters << std::endl;
  return nClusters;
};

/**
 * Compute position and total weight of a given cluster
 * Position is computed by weighted average of points positions (for x and y). For z just use the z of all hits (same since 2D cluster)
 * Weights are computed using an affine log scale of hits weights (ie hit energy).
 * \param positionDeltaRho2 ignore cells that have the distance squared to the highest energy cell in layer cluster greater than this squared distance
*/
void calculatePosition(const PointsCloud & points, Cluster & cl, float positionDeltaRho2) {
  constexpr float thresholdW0 = 2.9f;
  float total_weight = 0.f;
  float x = 0.f;
  float y = 0.f;

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
  float x_log = 0.f; ///< Sum of all x positions of clusters weighted 
  float y_log = 0.f;
  for (auto i : cl.hits()) {
//    //for silicon only just use 1+6 cells = 1.3cm for all thicknesses
    if (::distance(points, i, maxEnergyIndex) > std::sqrt(positionDeltaRho2)) // use ::distance to avoid colliding with std::distance
      continue;
    float rhEnergy = points.weight[i];
    float Wi = std::fmax(thresholdW0 + std::log(rhEnergy / total_weight), 0.);
    x_log += points.x[i] * Wi;
    y_log += points.y[i] * Wi;
    total_weight_log += Wi;
  }

  total_weight = total_weight_log;
  x = x_log;
  y = y_log;
  if (total_weight != 0.) {
    float inv_tot_weight = 1.f / total_weight;
    cl.setPosition(x * inv_tot_weight,
                   y * inv_tot_weight,
                   points.z[cl.hits()[0]]);
  } else {
    return cl.setPosition(0.f, 0.f, 0.f);
  }
}

/**
 * Build all Cluster objects holding cluster information
 * \param totalClusters the nb of clusters generated by the algorithm (return value of \a findAndAssign_clusters )
 * \param points the PointsCloud ( *TODO* why is this passed by value ?)
*/
std::vector<Cluster> getClusters(int totalClusters, PointsCloud points, ClueAlgoParameters const& params) {
  std::vector<Cluster> clusters(totalClusters);

  // loop over all points
  for (unsigned int i = 0; i < points.n; i++) {
    // If it has a valid cluster index, assign that point to the corresponding cluster
    if (points.clusterIndex[i] != -1) {
      auto & thisCluster = clusters[points.clusterIndex[i]];
      thisCluster.addCell(i);
//      thisCluster.addEnergyAndRescale(points.weight[i], points.z[i]);
      thisCluster.addEnergy(points.weight[i]);
      thisCluster.setLayer(points.layer[i]);
    }
  }

  // Now that we have all clusters, compute their barycenter
  for (auto & cl : clusters) {
    calculatePosition(points, cl, params.positionDeltaRho2);
  }

  return clusters;
}

#endif
