#ifndef LayerTiles_h
#define LayerTiles_h

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <memory>

#include "LayerTilesConstants.h"

/**
 * 2D layer of tiles, implementing binning
 * \tparam TilesConstants can be LayerTilesConstants (for CLUE3D) or TilesConstants (for CLUE)
*/
template<typename TilesConstantsT>
class LayerTiles {
 public:
 typedef TilesConstantsT TilesConstants;
  LayerTiles() {
    layerTiles_.resize(TilesConstants::nColumns *
                       TilesConstants::nRows);
  }

  void fill(const std::vector<float>& x, const std::vector<float>& y) {
    auto cellsSize = x.size();
    for (unsigned int i = 0; i < cellsSize; ++i) {
      layerTiles_[getGlobalBin(x[i], y[i])].push_back(i);
    }
  }

  ///< Add to the bin corresponding to coordinates (x;y) the index of hit i
  void fill(float x, float y, int i) {
    layerTiles_[getGlobalBin(x, y)].push_back(i);
  }

  int getXBin(float x) const {
    constexpr float xRange =
        TilesConstants::maxX - TilesConstants::minX;
    static_assert(xRange >= 0.);
    int xBin = (x - TilesConstants::minX) * TilesConstants::rX;
    xBin = std::min(xBin, TilesConstants::nColumns - 1);
    xBin = std::max(xBin, 0);
    return xBin;
  }

  int getYBin(float y) const {
    constexpr float yRange =
        TilesConstants::maxY - TilesConstants::minY;
    static_assert(yRange >= 0.);
    int yBin = (y - TilesConstants::minY) * TilesConstants::rY;
    yBin = std::min(yBin, TilesConstants::nRows - 1);
    yBin = std::max(yBin, 0);
    return yBin;
  }

  int getGlobalBin(float x, float y) const {
    return getXBin(x) + getYBin(y) * TilesConstants::nColumns;
  }

  ///< Map bin (x;y) coordinates to internal global bin number
  int getGlobalBinByBin(int xBin, int yBin) const {
    return xBin + yBin * TilesConstants::nColumns;
  }

  ///< Get box in bin numbers from box in space coordinates
  std::array<int, 4> searchBox(float xMin, float xMax, float yMin, float yMax) {
    int xBinMin = getXBin(xMin);
    int xBinMax = getXBin(xMax);
    int yBinMin = getYBin(yMin);
    int yBinMax = getYBin(yMax);
    return std::array<int, 4>({{xBinMin, xBinMax, yBinMin, yBinMax}});
  }

  void clear() {
    for (auto& t : layerTiles_) {
      t.clear();
    }
  }

  std::vector<int>& operator[](int globalBinId) {
    return layerTiles_[globalBinId];
  }

 private:
  /**
   * Outer vector : indexed by global bin number
   * Inner vector : list of hit ID in bin
  */
  std::vector<std::vector<int>> layerTiles_;
};

using LayerTilesClue = LayerTiles<LayerTilesConstantsClue>;
using LayerTilesClue3D = LayerTiles<LayerTilesConstantsClue3D>;

#endif  // LayerTiles_h
