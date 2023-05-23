#define Runclustering_cxx
#include "Runclustering.h"

#include <TF1.h>
#include <math.h>

#include <cstring>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>

#include "CLUEAlgo.h"
#include "CLUE3DAlgo.h"
#include "Math/SMatrix.h"
#include "Math/SVector.h"

#include "OptionParser.h"

namespace option = ROOT::option;
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
      << " pointType = " << points.pointType[c]
      << " clusterIdx = " << points.clusterIndex[c]
      << " nearestHigher = " << points.nearestHigher[c]
      << std::endl;
  }
}

/**
 * Opens file fileListPath and return a vector containing all lines (except empty lines)
*/
std::vector<std::string> readFileList(const char* fileListPath)
{
  std::vector<std::string> listOfFiles;
  std::ifstream file(fileListPath);
  if (!file.is_open()) {
    std::cerr << "** ERROR: Can't open '" << fileListPath << "' for input"
              << std::endl;
    throw std::runtime_error("File list opening failed");
  }

  std::string str; 
  while (std::getline(file, str))
  {
    if(str.find_first_not_of(' ') != std::string::npos)
    {
      //This line does not have only whitespace
      listOfFiles.push_back(str);
    }
  }
  return listOfFiles;
}

struct Arg: public option::Arg
{
  static void printError(const char* msg1, const option::Option& opt, const char* msg2)
  {
    fprintf(stderr, "ERROR: %s", msg1);
    fwrite(opt.name, opt.namelen, 1, stderr);
    fprintf(stderr, "%s", msg2);
  }
/*  
  static option::ArgStatus Unknown(const option::Option& option, bool msg)
  {
    if (msg) printError("Unknown option '", option, "'\n");
    return option::ARG_ILLEGAL;
  }
 
  static option::ArgStatus Required(const option::Option& option, bool msg)
  {
    if (option.arg != 0)
      return option::ARG_OK;
 
    if (msg) printError("Option '", option, "' requires an argument\n");
    return option::ARG_ILLEGAL;
  }
*/
  static option::ArgStatus RequiredDatatype(const option::Option& option, bool msg)
  {
    if (option.arg != 0) {
      std::string arg_str(option.arg);
      if (arg_str == "data" || arg_str == "simulation")
        return option::ARG_OK;
      else if (msg)
        printError("Option '", option, "' allows as arguments one of : data, simulation\n");
      return option::ARG_ILLEGAL;
    }
    if (msg) printError("Option '", option, "' requires an argument\n");
    return option::ARG_ILLEGAL;
  }
 
  static option::ArgStatus NonEmpty(const option::Option& option, bool msg)
  {
    if (option.arg != 0 && option.arg[0] != 0)
      return option::ARG_OK;
 
    if (msg) printError("Option '", option, "' requires a non-empty argument\n");
    return option::ARG_ILLEGAL;
  } 
 
  static option::ArgStatus Numeric(const option::Option& option, bool msg)
  {
    char* endptr = 0;
    if (option.arg != 0 && strtol(option.arg, &endptr, 10)){};
    if (endptr != option.arg && *endptr == 0)
      return option::ARG_OK;
 
    if (msg) printError("Option '", option, "' requires a numeric argument\n");
    return option::ARG_ILLEGAL;
  }

  static option::ArgStatus RequiredFloat(const option::Option& option, bool msg)
  {
    if (option.arg != 0){
      char* endptr = 0;
      std::strtof(option.arg, &endptr);
      if (endptr != option.arg) {
        //String to float conversion was successful
        return option::ARG_OK;
      }
    }
    if (msg) printError("Option '", option, "' requires a floating point argument\n");
    return option::ARG_ILLEGAL;
  }
};


enum  optionIndex { UNKNOWN, HELP, OUTPUT_FILE, INPUT_FILE_LIST, DATATYPE, SHIFT_RECHITS, FILTER_DWC,
  CLUE_DC, CLUE_RHOC, CLUE_OUTLIER_DELTA_FACTOR, CLUE_POSITION_DELTA_RHO2,
  CLUE3D_DENSITY_XY_DISTANCE, CLUE3D_CRITICAL_XY_DISTANCE, CLUE3D_CRITICAL_DENSITY, CLUE3D_OUTLIER_DELTA_FACTOR, CLUE3D_DENSITY_SIBLING_LAYERS, CLUE3D_NEAREST_HIGHER_SAME_LAYER, 
    CLUE3D_CRITICAL_Z_DISTANCE, CLUE3D_CRITICAL_SELF_DENSITY, CLUE3D_DENSITY_SAME_LAYER, CLUE3D_KERNEL_DENSITY,
    CLUE3D_MIN_CLUSTER_SIZE};
enum optionToggle { ENABLE, DISABLE};
const option::Descriptor usage[] =
{// index type shortopt longopt check_arg help
 {UNKNOWN, 0,"" , ""    ,option::Arg::None, "USAGE: runclustering [options]\n\n"
                                            "Non-options : input files (can also be specified by -f parameter)\n"
                                            "Options:" },
 {HELP,    0,"" , "help",option::Arg::None, "  --help  \tPrint usage and exit." },
 {OUTPUT_FILE, 0,"o", "output-file",Arg::NonEmpty, "-o, --output-file=  \tLocation of output file" },
 {INPUT_FILE_LIST, 0, "f", "input-file-list", Arg::NonEmpty, "-f, --input-file-list= \t Path to a file holding the paths of all the input files to read"},
 
 {DATATYPE, 0, "", "datatype", Arg::RequiredDatatype, "--datatype : specify data or simulation (used for DWC filtering)"},
 {SHIFT_RECHITS, ENABLE, "", "shift-rechits", Arg::None, "--shift-rechits : use ce_clean_x_shifted and impactX_unshifted columns (suitable for data). (auto-set if datatype=data)"},
 {SHIFT_RECHITS, DISABLE, "", "no-shift-rechits", Arg::None, "--no-shift-rechits :  use ce_clean_x_unshifted and impactX_unshifted columns (suitable only for Monte Carlo) (auto-set if datatype=simulation)"},
 {FILTER_DWC, ENABLE, "", "filter-dwc", Arg::None, "--filter-dwc : Filter events passing DWC cuts (Track Chisquare and 2x2cm impact position)"},
 {FILTER_DWC, DISABLE, "", "no-filter-dwc", Arg::None, "--no-filter-dwc : Do not filter events passing DWC cuts (Track Chisquare and 2x2cm impact position)"},
 
 {CLUE_DC, 0, "", "clue-deltac", Arg::RequiredFloat, "--clue-deltac= \t CLUE2D critical distance parameter"},
 {CLUE_RHOC, 0, "", "clue-rhoc", Arg::RequiredFloat, "--clue-rhoc= \t CLUE2D critical density parameter"},
 {CLUE_OUTLIER_DELTA_FACTOR, 0, "", "clue-outlier-factor", Arg::RequiredFloat, "--clue-outlier-factor= \t CLUE2D outlier delta factor"},
 {CLUE_POSITION_DELTA_RHO2, 0, "", "clue-position-delta-rho2", Arg::RequiredFloat, "--clue-position-delta-rho2= \t CLUE2D max distance squared to look for cells away from highest energy cell when computing cluster position"},
 
 {CLUE3D_DENSITY_XY_DISTANCE, 0, "", "clue3d-densityXYDistaeeceSqr", Arg::RequiredFloat, "--clue3d-densityXYDistanceSqr= \t CLUE3D distance squared (in cm^2) on the transverse plane to consider for local density"},
 {CLUE3D_CRITICAL_XY_DISTANCE, 0, "", "clue3d-criticalXYDistance", Arg::RequiredFloat, "--clue3d-criticalXYDistance= \t CLUE3D Minimal distance in cm on the XY plane from nearestHigher to become a seed"},
 {CLUE3D_CRITICAL_DENSITY, 0, "", "clue3d-criticalDensity", Arg::RequiredFloat, "--clue3d-criticalDensity= \t CLUE3D critical density parameter"},
 {CLUE3D_OUTLIER_DELTA_FACTOR, 0, "", "clue3d-outlier-factor", Arg::RequiredFloat, "--clue3d-outlier-factor= \t CLUE3D outlier delta factor"},
 {CLUE3D_DENSITY_SIBLING_LAYERS, 0, "", "clue3d-density-sibling-layers", Arg::Numeric, "--clue3d-density-sibling-layers= \t CLUE3D density sibling layers parameters, define range of layers +- layer# to look for another 2D cluster"},
 {CLUE3D_DENSITY_SAME_LAYER, ENABLE, "", "densityOnSameLayer", Arg::None, "--densityOnSameLayer \t CLUE3D : Consider layer clusters on the same layer when computing local energy density"},
 {CLUE3D_DENSITY_SAME_LAYER, DISABLE, "", "no-densityOnSameLayer", Arg::None, "--no-densityOnSameLayer \t CLUE3D : Do not consider layer clusters on the same layer when computing local energy density"},
 {CLUE3D_NEAREST_HIGHER_SAME_LAYER, ENABLE, "", "nearestHigherOnSameLayer", Arg::None, "--nearestHigherOnSameLayer \t CLUE3D : Allow the nearestHigher to be located on the same layer"},
 {CLUE3D_NEAREST_HIGHER_SAME_LAYER, DISABLE, "", "no-nearestHigherOnSameLayer", Arg::None, "--no-nearestHigherOnSameLayer \t CLUE3D : Do not allow the nearestHigher to be located on the same layer"},
 {CLUE3D_CRITICAL_Z_DISTANCE, 0, "", "clue3d-criticalZDistanceLyr", Arg::Numeric, "--clue3d-criticalZDistanceLyr= \t CLUE3D Minimal distance in layers along the Z axis from nearestHigher to become a seed"},
 {CLUE3D_CRITICAL_SELF_DENSITY, 0, "", "clue3d-criticalSelfDensity", Arg::RequiredFloat, "--clue3d-criticalSelfDensity= \t CLUE3D Minimum ratio of self_energy/local_density to become a seed."},
 {CLUE3D_KERNEL_DENSITY, 0, "", "clue3d-kernelDensityFactor", Arg::RequiredFloat, "--clue3d-kernelDensityFactor= \t CLUE3D Kernel factor to be applied to other layer clusters while computing the local density"},

 // This is not a parameter for CLUE3D itself, but a replication of the filtering of layer clusters in CMSSW before running CLUE3D
 {CLUE3D_MIN_CLUSTER_SIZE, 0, "", "clue3d-minLayerClusterSize", Arg::Numeric, "--clue3d-minLayerClusterSize= \t CLUE3D : Minimum size of layer clusters to consider (inclusive)"},

 {UNKNOWN, 0,"" ,  ""   ,option::Arg::None, "\nExamples:\n"
                                            "  ./runclustering -f files-single.txt -o ./CLUE_clusters.root\n"
                                            "  ./runclustering --clue-rhoc=2. -f files-single.txt -o ./CLUE_clusters.root\n" },
 {0,0,0,0,0,0}
};


int main(int argc, char *argv[]) 
{
  argc-=(argc>0); argv+=(argc>0); // skip program name argv[0] if present
  option::Stats  stats(usage, argc, argv);
  std::vector<option::Option> options(stats.options_max);
  std::vector<option::Option> buffer(stats.buffer_max);
  option::Parser parse(usage, argc, argv, &options[0], &buffer[0]);
 
  if (parse.error())
    return 1;
  
  if (options[HELP] || argc == 0 || options[UNKNOWN]) {
    option::printUsage(std::cout, usage);
    return 0;
  }
  if (!options[FILTER_DWC]) {
    cerr << "Either --filter-dwc or --no-filter-dwc must be specified" << endl;
    option::printUsage(std::cout, usage);
    return 1;
  }

  cout << "Datatype : " << options[DATATYPE].arg << endl;

  ClueAlgoParameters clueParameters;
  if (options[CLUE_DC])
    clueParameters.deltac = {std::stof(options[CLUE_DC].arg), -1.};
  else
    clueParameters.deltac = {1.3f, 3.f * sqrt(2.f) + 0.1f}; // (cm) No idea where AHCAL value comes from (CE-E value is from CMSSW)
  
  constexpr float MIP2GeV[3] = {0.0105, 0.0812, 0.12508};
  if (options[CLUE_RHOC])
    clueParameters.rhoc = {std::stof(options[CLUE_RHOC].arg), -1.};
  else {
    /* In test beam : 1 sigma noise is ~ 5 to 7 ADC HG counts
    From DN19-019, figure 3, 4000 HG ADC counts is 100 MIPs
    so with kappa=9, 9 sigma noise is 9sigma * 6 ADC counst *  100 MIP / 4000 ADC  = 1.35 MIPs
    (about 13 MeV)
    (values for AHCAL were just copied over from CE-E)
    */
    clueParameters.rhoc = {1.35f * MIP2GeV[0], 1.35f * MIP2GeV[2]};
    // Old values : 
    //clueParameters.rhoc = {4.f * MIP2GeV[0], 4.f * MIP2GeV[2]};
  }
  
  if (options[CLUE_OUTLIER_DELTA_FACTOR])
    clueParameters.outlierDeltaFactor = std::stof(options[CLUE_OUTLIER_DELTA_FACTOR].arg);
  else
    clueParameters.outlierDeltaFactor = 2.f;
  
  if (options[CLUE_POSITION_DELTA_RHO2])
    clueParameters.positionDeltaRho2 = std::stof(options[CLUE_POSITION_DELTA_RHO2].arg);
  else
    clueParameters.positionDeltaRho2 = 1.69f;
  
  cout << "Using CLUE parameters : " << clueParameters << endl;

  Clue3DAlgoParameters clue3DParameters;
  if (options[CLUE3D_DENSITY_XY_DISTANCE])
    clue3DParameters.densityXYDistanceSqr = std::stof(options[CLUE3D_DENSITY_XY_DISTANCE].arg);
  else
    clue3DParameters.densityXYDistanceSqr = 3.24f; // = 2.6 cm * 2.6 cm, copied from CMSSW
  
  if (options[CLUE3D_CRITICAL_XY_DISTANCE])
    clue3DParameters.criticalXYDistance = std::stof(options[CLUE3D_CRITICAL_XY_DISTANCE].arg);
  else
    clue3DParameters.criticalXYDistance = 1.8f; //cm, copied from CMSSW

  if (options[CLUE3D_CRITICAL_Z_DISTANCE])
    clue3DParameters.criticalZDistanceLyr = std::stoi(options[CLUE3D_CRITICAL_Z_DISTANCE].arg);
  else
    clue3DParameters.criticalZDistanceLyr = 5; //Layer count, copied from CMSSW
  
  if (options[CLUE3D_CRITICAL_DENSITY])
    clue3DParameters.criticalDensity = std::stof(options[CLUE3D_CRITICAL_DENSITY].arg);
  else
    clue3DParameters.criticalDensity = 0.6f; // From CMSSW RecoHGCal/TICL/python/CLUE3DHighStep_cff.py
  
  if (options[CLUE3D_OUTLIER_DELTA_FACTOR])
    clue3DParameters.outlierDeltaFactor = std::stof(options[CLUE3D_OUTLIER_DELTA_FACTOR].arg);
  else
    clue3DParameters.outlierDeltaFactor = 2.f; // from CMSSW (called outlierMultiplier there)

  if (options[CLUE3D_DENSITY_SIBLING_LAYERS]) {
    clue3DParameters.densitySiblingLayers = std::stoi(options[CLUE3D_DENSITY_SIBLING_LAYERS].arg, nullptr, 10);
  } else
    clue3DParameters.densitySiblingLayers = 3; // from CMSSW
  
  if (options[CLUE3D_DENSITY_SAME_LAYER])
    clue3DParameters.densityOnSameLayer = options[CLUE3D_DENSITY_SAME_LAYER].last()->type() == ENABLE;
  else
    clue3DParameters.densityOnSameLayer = false; // from CMSSW

  if (options[CLUE3D_NEAREST_HIGHER_SAME_LAYER])
    clue3DParameters.nearestHigherOnSameLayer = options[CLUE3D_NEAREST_HIGHER_SAME_LAYER].last()->type() == ENABLE;
  else
    clue3DParameters.nearestHigherOnSameLayer = false; // from CMSSW
  
  if (options[CLUE3D_CRITICAL_SELF_DENSITY])
    clue3DParameters.criticalSelfDensity = std::stof(options[CLUE3D_CRITICAL_SELF_DENSITY].arg);
  else
    clue3DParameters.criticalSelfDensity = 0.15f; // from CMSSW
  
  if (options[CLUE3D_KERNEL_DENSITY])
    clue3DParameters.kernelDensityFactor = std::stof(options[CLUE3D_KERNEL_DENSITY].arg);
  else
    clue3DParameters.kernelDensityFactor = 0.2f; // from CMSSW
  
  cout << "Using CLUE3D parameters : " << clue3DParameters << endl;

  // Shifting rechits
  std::string datatype(options[DATATYPE].arg);
  bool shiftRechits;
  if (options[SHIFT_RECHITS]) {
    shiftRechits = (options[SHIFT_RECHITS].last()->type() == ENABLE);
  }
  else {
    shiftRechits = (datatype == "data");
  }
  if (shiftRechits)
    cout << "Shifting rechits positions (suitable for data only)" << endl;
  else
    cout << "Not shifting rechits positions (suitable for simulation only)" << endl;
  
  if (options[FILTER_DWC].last()->type() == ENABLE)
    cout << "Filtering events passing DWC cuts (Track Chisquare and 2x2cm impact position)" << endl;
  else
    cout << "Not filtering events passing DWC cuts (Track Chisquare and 2x2cm impact position)" << endl;

  cout << endl;

  int filterMinLayerClusterSize = 2; // From CMSSW : RecoHGCal/TICL/python/CLUE3DHighStep_cff.py
  if (options[CLUE3D_MIN_CLUSTER_SIZE])
    filterMinLayerClusterSize = std::stoul(options[CLUE3D_MIN_CLUSTER_SIZE].arg);

  // Dealing with input files
  if (parse.nonOptionsCount() > 0 && options[INPUT_FILE_LIST]) {
    cerr << "You cannot pass an input file list as a file at the same time as input files on the command line ! Aborting" << endl;
    return 1;
  }
  std::vector<std::string> listOfFilePaths;
  for (int i = 0; i < parse.nonOptionsCount(); i++) {
    listOfFilePaths.push_back(parse.nonOption(i));
  }
  if (options[INPUT_FILE_LIST]) {
    listOfFilePaths = readFileList(options[INPUT_FILE_LIST].arg);
  }

  Runclustering tbCLUS(std::move(listOfFilePaths), options[OUTPUT_FILE].arg, clueParameters, clue3DParameters,
    datatype, shiftRechits, options[FILTER_DWC].last()->type() == ENABLE);
  tbCLUS.EventLoop(filterMinLayerClusterSize);
  return 0;
}

/**
 * \param filterMinLayerClusterSize only consider for CLUE3D layer clusters whose size is greater or equal than filterMinLayerClusterSize
*/
void Runclustering::EventLoop(unsigned filterMinLayerClusterSize) {

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
  output_file_->WriteObject<ClueAlgoParameters>(&clueParams_, "clueParams");
  output_file_->WriteObject<Clue3DAlgoParameters>(&clue3DParams_, "clue3DParams");

  TTree clusters_tree("clusters", "clusters");

  bool DWC_passes_cuts = 0;

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
  clusters_tree.Branch("event", &event);
  //clusters_tree.Branch("NRechits", &NRechits);

  if (fChain->GetBranch("trueBeamEnergy")) { //If branch exists in input dataset, copy it to output
    clusters_tree.Branch("trueBeamEnergy", &trueBeamEnergy);
  }

  // Delay wire chambers data
  clusters_tree.Branch("impactX", &impactX); //These are vector<float> of size 40 (nb of layers)
  clusters_tree.Branch("impactY", &impactY);
  clusters_tree.Branch("DWC_passes_cuts", &DWC_passes_cuts); // If true, event DWC value passes cuts

  clusters_tree.Branch("rechits_x", &pcloud.x);
  clusters_tree.Branch("rechits_y", &pcloud.y);
  clusters_tree.Branch("rechits_z", &pcloud.z);
  clusters_tree.Branch("rechits_energy", &pcloud.weight);
  clusters_tree.Branch("rechits_energy_MIP", &ce_clean_energy_MIP);
  clusters_tree.Branch("rechits_layer", &pcloud.layer);
  clusters_tree.Branch("rechits_rho", &pcloud.rho);
  clusters_tree.Branch("rechits_delta", &pcloud.delta);
  clusters_tree.Branch("rechits_nearestHigher", &pcloud.nearestHigher);
  clusters_tree.Branch("rechits_pointType", &pcloud.pointType);

  clusters_tree.Branch("clus2D_x", &clusters_soa.x);
  clusters_tree.Branch("clus2D_y", &clusters_soa.y);
  clusters_tree.Branch("clus2D_z", &clusters_soa.z);
  clusters_tree.Branch("clus2D_energy", &clusters_soa.energy);
  clusters_tree.Branch("clus2D_layer", &clusters_soa.layer);
  clusters_tree.Branch("clus2D_size", &clusters_soa.size);
  clusters_tree.Branch("clus2D_idxs", &clusters_soa.hitidxs);
  clusters_tree.Branch("clus2D_rho", &pcloud2d.rho);
  clusters_tree.Branch("clus2D_delta", &pcloud2d.delta);
  clusters_tree.Branch("clus2D_nearestHigher", &pcloud2d.nearestHigher);
  clusters_tree.Branch("clus2D_pointType", &pcloud2d.pointType);

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

    // ---- Delay Wire Chamber filtering
    DWC_passes_cuts = DWC_trackChi2_X < 10 && DWC_trackChi2_Y < 10;
    /* These selection values were taken from Matteo Bonanomi's code (personal communication, as the version on GitHub is older and outdated ?)
    the b_x and b_y branches were never changed in previous steps so no need to mirror/unmirror 
    */
    if (datatype == "data") {
      DWC_passes_cuts = DWC_passes_cuts && (std::abs(DWC_b_x + 2.7)<1.) && (std::abs(DWC_b_y - 1.)<1.);
    } else if (datatype == "simulation") {
      float xcorr = - 3.6 + 2.7; //x MC
      float ycorr = + 2.6 - 1.0 ; //y MC
      // Note the minus sign
      DWC_passes_cuts = DWC_passes_cuts && (std::abs(-DWC_b_x + xcorr)<1.) && (std::abs(-DWC_b_y + ycorr)<1.);
    } else {
      assert(false && "datatype should be either data or simulation");
    }

    if (filterDwc && !DWC_passes_cuts)
      continue;

    // ----- rescaling energy from MeV to GeV
    vector<float> scaled_energy;
    scaled_energy = *ce_clean_energy_MeV;
    double en_scale{0.001};
    std::transform(scaled_energy.begin(), scaled_energy.end(), scaled_energy.begin(), [&en_scale](auto& c){return c*en_scale;});



    // Compute clusters using CLUE
    //NLAYERS defined in LayersTilesConstants
    std::array<LayerTilesClue, NLAYERS> tiles; ///< Array of LayerTiles (2D layer of tiles) for each layer
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
    calculate_density(tiles, pcloud, clueParams_.deltac);
    // Calculate nearest higher density  point
    calculate_distanceToHigher(tiles, pcloud, clueParams_.outlierDeltaFactor, clueParams_.deltac);
    // get seeds and followers
    auto total_clusters = findAndAssign_clusters(pcloud, clueParams_.outlierDeltaFactor, clueParams_.deltac, clueParams_.rhoc);
    std::vector<Cluster> clusters = getClusters(total_clusters, pcloud, clueParams_);
    // Fill in the clusters_SoA and masking for CLUE3D
    clusters_soa.load(clusters, filterMinLayerClusterSize);

    // Compute clusters using CLUE3D [similar sequence as the 2D part, starting pointcloud here composed of the 2D clusters just made above]
    std::array<LayerTilesClue3D, NLAYERS> tiles2d;
    pcloud2d.x = clusters_soa.x ;
    pcloud2d.y = clusters_soa.y ;
    pcloud2d.z = clusters_soa.z ;
    vector<int> clusters_soao ; ///< Layer of each 2D cluster, indexed by cluster ID
    clusters_soao = clusters_soa.layer;
    std::vector<unsigned int> clusters_soau(std::begin(clusters_soao), std::end(clusters_soao)); // convert to unsigned int
    pcloud2d.layer = clusters_soau ;
    pcloud2d.weight = clusters_soa.energy ;
    pcloud2d.masked = clusters_soa.masked;
    pcloud2d.resizeOutputContainers(clusters_soa.x.size());

    compute_histogram(tiles2d, pcloud2d);



    calculate_density3d(tiles2d, pcloud2d, clue3DParams_);
    calculate_distanceToHigher3d(tiles2d, pcloud2d, clue3DParams_);
    auto total_clusters3d = findAndAssign_clusters3d(pcloud2d, clue3DParams_);
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

