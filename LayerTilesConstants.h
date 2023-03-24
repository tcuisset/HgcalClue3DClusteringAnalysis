#ifndef LayerTilesConstants_h
#define LayerTilesConstants_h

#include <array>
#include <cstdint>

#define NLAYERS 50

static constexpr int32_t my_ceil(float num) {
return (static_cast<float>(static_cast<int32_t>(num)) == num)
          ? static_cast<int32_t>(num)
          : static_cast<int32_t>(num) + ((num > 0) ? 1 : 0);
}

struct LayerTilesConstantsClue {
  static constexpr float minX = -10.f;
  static constexpr float maxX = 10.f;
  static constexpr float minY = -10.f;
  static constexpr float maxY = 10.f;

  static constexpr float tileSize = 5.f; ///< Copied from CMSSW (has to be the same, this has physical consequences in findNearestHigher)
  static constexpr int nColumns = my_ceil((maxX - minX) / tileSize);
  static constexpr int nRows = my_ceil((maxY - minY) / tileSize);
  static constexpr int maxTileDepth = 40;

  static constexpr float rX = nColumns / (maxX - minX);
  static constexpr float rY = nRows / (maxY - minY);
};


struct LayerTilesConstantsClue3D {
  static constexpr float minX = -10.f;
  static constexpr float maxX = 10.f;
  static constexpr float minY = -10.f;
  static constexpr float maxY = 10.f;
  /**
   * The tile size is important in calculate_distanceToHigher in both CLUE and CLUE3D
   * This value was chosen so that our square x-y bins have approximately the same area as eta-phi bins in CLUE3D
   * Though it is a somewhat arbitray choice since the area of bins in CMSSW CLUE3D varies depending on where you are in the detector
  */
  static constexpr float tileSize = 3.7f; 
  static constexpr int nColumns = my_ceil((maxX - minX) / tileSize);
  static constexpr int nRows = my_ceil((maxY - minY) / tileSize);
  static constexpr int maxTileDepth = 40;

  static constexpr float rX = nColumns / (maxX - minX);
  static constexpr float rY = nRows / (maxY - minY);
};

#endif  // LayerTilesConstants_h
