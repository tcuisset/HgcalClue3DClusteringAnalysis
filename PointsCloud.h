#ifndef Points_Cloud_h
#define Points_Cloud_h

#include <vector>
#include "Cluster.h"



/**
 * List of points, holding for each point its position, layer nb, etc 
 * as well as output of CLUE algorithm for each point (rho, delta, nearest higher, ...)
*/
struct PointsCloud {
  enum PointType : char {
    FOLLOWER = 0,
    SEED = 1,
    OUTLIER = 2
  };

  PointsCloud() = default;

  ///< 
  void resizeOutputContainers(unsigned int const& nPoints) {
    masked.resize(nPoints, 0);
    rho.resize(nPoints, 0.);
    delta.resize(nPoints, -1.);
    nearestHigher.resize(nPoints, -1);
    clusterIndex.resize(nPoints, -1);
    followers.resize(nPoints);
    pointType.resize(nPoints, PointsCloud::FOLLOWER);
    n = nPoints;
  }

  void inline reset() {
    x.clear();
    y.clear();
    z.clear();
    layer.clear();
    weight.clear();
    masked.clear();
    rho.clear();
    delta.clear();
    nearestHigher.clear();
    followers.clear();
    pointType.clear();
    clusterIndex.clear();
  }

  // Input variables
  std::vector<float> x;
  std::vector<float> y;
  std::vector<float> z;
  std::vector<unsigned int> layer;
  std::vector<float> weight; ///< Weights of all points (ie energy)

  std::vector<char> masked; ///< If 1, then point is masked (should be ignored). 0 otherwise

  // Output variables : vectors of size n
  std::vector<float> rho; ///< Local energy density
  std::vector<float> delta; ///< Distance to nearest highest
  std::vector<int> nearestHigher; ///< ID of nearest highest (-1 in case no nearest higher)
  std::vector<std::vector<int>> followers; ///< List of points that follow this point (ie points which are neither seeds nor outliers and for which we are the nearest higher)
  std::vector<char> pointType; ///< Is the point a follower, a seed or an outlier (use enum PointType for values)
  std::vector<int> clusterIndex; ///< ID of the cluster this point is member of

  unsigned int n; ///< Size of output variables vectors, usually same as input
};


/**
 * STructure holding all informations for all clusters (position of cluster, .., and all hits ID contained in each cluster)
*/
struct ClustersSoA {
  ClustersSoA() = default;

  /**
   * \param filterMinClusterSize mask clusters whose size is strictly less than filterMinClusterSize
  */
  void inline load(const std::vector<Cluster>& clusters, unsigned filterMinClusterSize=0) {
    auto total_clusters = clusters.size();
    x.clear();
    y.clear();
    z.clear();
    energy.clear();
    layer.clear();
    size.clear();
    hitidxs.clear();
    masked.clear();
    x.reserve(total_clusters);
    y.reserve(total_clusters);
    z.reserve(total_clusters);
    energy.reserve(total_clusters);
    layer.reserve(total_clusters);
    size.reserve(total_clusters);
    hitidxs.reserve(total_clusters); 
    masked.reserve(total_clusters);

    for (auto const & cl :  clusters) {
//      const auto [cl_x, cl_y, cl_z] = cl.position();
      if (cl.hits().size() >= filterMinClusterSize) {
        auto cl_pos = cl.position();
        x.push_back(std::get<0>(cl_pos));
        y.push_back(std::get<1>(cl_pos));
        z.push_back(std::get<2>(cl_pos));
        energy.push_back(cl.energy());
        layer.push_back(cl.layer());
        size.push_back(cl.hits().size());
        hitidxs.push_back(cl.hits());
        masked.push_back(cl.hits().size() >= filterMinClusterSize ? 0 : 1);
      }
    }
  };

  std::vector<float> x;
  std::vector<float> y;
  std::vector<float> z;
  std::vector<float> energy;
  std::vector<int> layer;
  std::vector<int> size; ///< Number of hits per cluster
  std::vector<std::vector<int>> hitidxs; ///< List of hits IDs per cluster
  std::vector<char> masked; ///< If 1 then cluster is masked

};

#endif
