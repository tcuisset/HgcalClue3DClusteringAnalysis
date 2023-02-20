#ifndef Points_Cloud_h
#define Points_Cloud_h

#include <vector>
#include "Cluster.h"


struct PointsCloud {
  PointsCloud() = default;

  void resizeOutputContainers(unsigned int const& nPoints) {
    rho.resize(nPoints, 0.);
    delta.resize(nPoints, -1.);
    nearestHigher.resize(nPoints, -1);
    clusterIndex.resize(nPoints, -1);
    followers.resize(nPoints);
    isSeed.resize(nPoints, 0);
    n = nPoints;
  }

  void inline reset() {
    x.clear();
    y.clear();
    z.clear();
    layer.clear();
    weight.clear();
    rho.clear();
    delta.clear();
    nearestHigher.clear();
    followers.clear();
    isSeed.clear();
    clusterIndex.clear();
  }

  // Input variables
  std::vector<float> x;
  std::vector<float> y;
  std::vector<float> z;
  std::vector<unsigned int> layer;
  std::vector<float> weight;

  // Output variables
  std::vector<float> rho;
  std::vector<float> delta;
  std::vector<int> nearestHigher;
  std::vector<std::vector<int>> followers;
  std::vector<int> isSeed;
  std::vector<int> clusterIndex;

  unsigned int n;
};

struct ClustersSoA {
  ClustersSoA() = default;

  void inline load(const std::vector<Cluster>& clusters) {
    auto total_clusters = clusters.size();
    x.resize(total_clusters, 0.f);
    y.resize(total_clusters, 0.f);
    z.resize(total_clusters, 0.f);
    energy.resize(total_clusters, 0.f);
    layer.resize(total_clusters, 0);
    size.resize(total_clusters, 0);
    hitidxs.resize(total_clusters); 

    int counter = 0;
    for (auto const & cl :  clusters) {
//      const auto [cl_x, cl_y, cl_z] = cl.position();
      auto cl_pos = cl.position();
      x[counter] = std::get<0>(cl_pos);
      y[counter] = std::get<1>(cl_pos);
      z[counter] = std::get<2>(cl_pos);
      energy[counter] = cl.energy();
      layer[counter] = cl.layer();
      size[counter] = cl.hits().size();
      hitidxs[counter] = cl.hits();
      counter++;
    }
  };

  std::vector<float> x;
  std::vector<float> y;
  std::vector<float> z;
  std::vector<float> energy;
  std::vector<int> layer;
  std::vector<int> size;
  std::vector<std::vector<int>> hitidxs;


};

#endif
