#ifndef CLUEAlgoParameters_h
#define CLUEAlgoParameters_h

#include <ostream>
#include <array>
#include <Rtypes.h>


struct ClueAlgoParameters
{
    std::array<float, 2> deltac; ///< Critical distance parameters, index 0 for HGCAL, index 1 for AHCAL
    std::array<float, 2> rhoc;
    float outlierDeltaFactor;
    /**
     * When computing layer cluster position,
     *  ignore cells that have the distance squared to the highest energy cell in layer cluster greater than this squared distance
    */
    float positionDeltaRho2;

    friend std::ostream& operator<< (std::ostream& stream, const ClueAlgoParameters& p) {
        stream << "deltac = " << p.deltac[0] 
               << ", rhoc = " << p.rhoc[0]
               << ", outlierDeltaFactor = " << p.outlierDeltaFactor
               << ", positionDeltaRho2 = " << p.positionDeltaRho2;
        return stream;
    }

    ClassDef(ClueAlgoParameters, 1)
};



struct Clue3DAlgoParameters
{
    float densityXYDistanceSqr; ///< in cm^2, 2.6*2.6, distance squared on the transverse plane to consider for local density
    float criticalXYDistance; ///< Minimal distance in cm on the XY plane from nearestHigher to become a seed
    int criticalZDistanceLyr;  ///< Minimal distance in layers along the Z axis from nearestHigher to become a seed
    float criticalDensity;  ///< Critical density parameters, index 0 for HGCAL, index 1 for AHCAL (called criticalDensity in CMSSW)
    float outlierDeltaFactor; ///< multiplicative factor to deltac to get distance to search for nearest higher
    int densitySiblingLayers; ///< define range of layers +- layer# of a point  
    bool densityOnSameLayer;  ///< Consider layer clusters on the same layer when computing local energy density
    float kernelDensityFactor; ///< Kernel factor to be applied to other LC while computing the local density
    bool nearestHigherOnSameLayer; ///< Allow the nearestHigher to be located on the same layer
    float criticalSelfDensity; ///< Minimum ratio of self_energy/local_density to become a seed. (roughly 1/(densitySiblingLayers+1) )

    friend std::ostream& operator<< (std::ostream& stream, const Clue3DAlgoParameters& p) {
        stream << "densityXYDistanceSqr = " << p.densityXYDistanceSqr
               << ", criticalXYDistance = " << p.criticalXYDistance
               << ", criticalZDistanceLyr = " << p.criticalZDistanceLyr
               << ", criticalDensity = " << p.criticalDensity
               << ", outlierDeltaFactor = " << p.outlierDeltaFactor 
               << ", densitySiblingLayers = " << p.densitySiblingLayers
               << ", densityOnSameLayer = " << p.densityOnSameLayer
               << ", kernelDensityFactor = " << p.kernelDensityFactor
               << ", nearestHigherOnSameLayer = " << p.nearestHigherOnSameLayer
               << ", criticalSelfDensity = " << p.criticalSelfDensity;
        return stream;
    }

    ClassDef(Clue3DAlgoParameters, 1)
};

#endif