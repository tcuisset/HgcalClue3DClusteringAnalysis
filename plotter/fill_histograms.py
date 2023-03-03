import os
import pickle
import argparse
import sys

import uproot
import awkward as ak
import hist
import pandas as pd


parser = argparse.ArgumentParser(description="Fill histograms from CLUE_clusters.root files for plotting. Writes pickled boost-histograms")
parser.add_argument("--input", dest="input_file", default='ClusteringAnalysis/CLUE_clusters_single.root',
    help="Complete path to CLUE output, usually named CLUE_clusters.root")
parser.add_argument("--output", dest="output_file",
    default='/home/llr/cms/cuisset/hgcal/testbeam18/clue3d-dev/src/plots/cache/hists.pkl',
    help="Path to output histograms (as a pickle file), usually names hists.pkl")
args = parser.parse_args()


#fname='ClusteringAnalysis/CLUE_clusters_single.root' ## with CLUS->LC->rechit mapping
#outputCacheHistsFolder = '/home/llr/cms/cuisset/hgcal/testbeam18/clue3d-dev/src/plots/cache'


beamEnergies = [20, 50, 80, 100, 120, 150, 200, 250, 300]
beamEnergiesAxis = hist.axis.IntCategory(beamEnergies, name="beamEnergy", label="Beam energy (GeV)")
layerAxis = hist.axis.Integer(start=0, stop=30, name="layer", label="Layer number")

rho_axis = hist.axis.Regular(bins=100, start=0, stop=10., name="rho", label="Energy density")
delta_axis = hist.axis.Regular(bins=100, start=0, stop=3., name="delta", label="Distance to nearest higher")
seed_axis = hist.axis.IntCategory([0, 1], name="isSeed") #0: not a seed, 1: is a seed

cluster2dEnergy_axis = hist.axis.Regular(100, 0, 1, name="cluster2dEnergy", label="2D cluster energy")

cluster3dSizeAxis = hist.axis.Integer(1, 40, name="cluster3dSize", label="3D cluster size")
cluster3dEnergyAxis = hist.axis.Regular(100, 0, 20, name="cluster3dEnergy", label="3D cluster energy")

xPosition_axis = hist.axis.Regular(bins=50, start=-10., stop=10., name="x", label="x")
yPosition_axis = hist.axis.Regular(bins=50, start=-10., stop=10., name="y", label="y")
zPosition_axis = hist.axis.Regular(bins=50, start=0., stop=60., name="z", label="z")
# An axis for difference in position (for example cluster spatial resolution)
diffX_axis = hist.axis.Regular(bins=50, start=-5., stop=5., name="diffX", label="x position difference")
diffY_axis = hist.axis.Regular(bins=50, start=-5., stop=5., name="diffY", label="y position difference")

energy_axis = hist.axis.Regular(bins=100, start=0., stop=10., name="energy", label="Reconstructed hit energy")

totalRecHitEnergy_axis = hist.axis.Regular(bins=500, start=0, stop=300, name="totalRecHitEnergy", label="Total RecHit energy per event (GeV)")

hist_dict = {}

hist_dict["rechits_position"] = hist.Hist(beamEnergiesAxis, layerAxis, energy_axis, xPosition_axis, yPosition_axis)
hist_dict["rechits_prop"] = hist.Hist(beamEnergiesAxis, layerAxis, seed_axis, rho_axis, delta_axis)

hist_dict["clusters2d_position"] = hist.Hist(beamEnergiesAxis, layerAxis, seed_axis, xPosition_axis, yPosition_axis)
hist_dict["clusters2d_prop"] = hist.Hist(beamEnergiesAxis, layerAxis, seed_axis, rho_axis, delta_axis)

hist_dict["clusters3d_position_xy"] = hist.Hist(beamEnergiesAxis, cluster3dEnergyAxis, xPosition_axis, yPosition_axis)
hist_dict["clusters3d_position_z"] = hist.Hist(beamEnergiesAxis, cluster3dEnergyAxis, zPosition_axis)
hist_dict["clusters3d_prop"] = hist.Hist(beamEnergiesAxis, cluster3dEnergyAxis, cluster3dSizeAxis)

hist_dict["event_prop"] = hist.Hist(beamEnergiesAxis, totalRecHitEnergy_axis)

# Spatial resolution (ie 2D cluster position minus particle impact position estimated from delay wire chambers) for all 3D clusters
# considering all 3D clusters equally
hist_dict["clue3d_clusters_spatial_res_count"] = hist.Hist(beamEnergiesAxis, layerAxis, cluster2dEnergy_axis, diffX_axis, diffY_axis)
#considering only highest energy 3D cluster per event
hist_dict["clue3d_clusters_spatial_res_largest_3D_cluster"] = hist.Hist(beamEnergiesAxis, layerAxis, cluster2dEnergy_axis, diffX_axis, diffY_axis)

#Total clustered energy per layer and per 3D cluster
hist_dict["clue3d_total_clustered_energy"] = hist.Hist(beamEnergiesAxis, layerAxis, cluster2dEnergy_axis)
#Same but only consider for each event the 3D cluster with the highest energy
hist_dict["clue3d_total_clustered_energy_by_largest_cluster"] = hist.Hist(beamEnergiesAxis, layerAxis, cluster2dEnergy_axis)

try:
    for (array, report) in uproot.iterate(args.input_file + ":clusters", step_size="100MB", library="ak", report=True):
        print("Processing events [" + str(report.start) + ", " + str(report.stop) + "[")

        impact = ak.to_dataframe(array[
            ["impactX", "impactY"]],
            levelname=lambda i : {0 : "event", 1:"layer"}[i])

        rechits = ak.to_dataframe(array[
            ["beamEnergy", "NRechits", "rechits_x", "rechits_y", "rechits_z", "rechits_energy", "rechits_layer",
            "rechits_rho", "rechits_delta", "rechits_isSeed"]], 
            levelname=lambda i : {0 : "event", 1:"rechit_id"}[i])

        hist_dict["rechits_position"].fill(rechits.beamEnergy, rechits.rechits_layer, rechits.rechits_energy, rechits.rechits_x, rechits.rechits_y)
        hist_dict["rechits_prop"].fill(rechits.beamEnergy, rechits.rechits_layer, rechits.rechits_isSeed, rechits.rechits_rho, rechits.rechits_delta)

        event_group = rechits.groupby(level=0) # group by event number (first level)
        hist_dict["event_prop"].fill(
            event_group.beamEnergy.first(), #beam energy is the same for all rows of an event, so take the first hit
            event_group.rechits_energy.sum()) #Sum rechits_energy for all rechits of an event


        clusters_2d = ak.to_dataframe(array[
            ["beamEnergy", "NRechits", "clus2D_x", "clus2D_y", "clus2D_z", "clus2D_energy", "clus2D_layer",
            "clus2D_rho", "clus2D_delta", "clus2D_idxs", "clus2D_isSeed"]
            ], 
            levelname=lambda i : {0 : "event", 1:"clus2D_id", 2:"hit_id"}[i])
        
        clusters2D_slice = clusters_2d.loc[(slice(None), slice(None), 0)].reset_index(drop=True) #Slice to get row of first hit (we do not care about the hits)
        hist_dict["clusters2d_prop"].fill(clusters2D_slice.beamEnergy, clusters2D_slice.clus2D_layer, clusters2D_slice.clus2D_isSeed, clusters2D_slice.clus2D_rho, clusters2D_slice.clus2D_delta)


        clusters_3d = ak.to_dataframe(array[
            ["beamEnergy", "NRechits", "clus3D_x", "clus3D_y", "clus3D_z", "clus3D_energy", "clus3D_size", "clus3D_idxs"]
            ], 
            levelname=lambda i : {0 : "event", 1:"clus3D_id", 2:"hit_id"}[i])
        
        clusters3D_slice = clusters_3d.loc[(slice(None), slice(None), 0)].reset_index(drop=True) #Slice to get row of first 2D cluster (we do not care about the hits)
        hist_dict["clusters3d_prop"].fill(clusters3D_slice.beamEnergy, clusters3D_slice.clus3D_energy, clusters3D_slice.clus3D_size)
        hist_dict["clusters3d_position_xy"].fill(clusters3D_slice.beamEnergy, clusters3D_slice.clus3D_energy, clusters3D_slice.clus3D_x, clusters3D_slice.clus3D_y)
        hist_dict["clusters3d_position_z"].fill(clusters3D_slice.beamEnergy, clusters3D_slice.clus3D_energy, clusters3D_slice.clus3D_z)

        #############  Merge clusters3D and clusters2D
        clusters_3d_2d = pd.merge(
            #Left : clusters3d
            # Reset multiindex so its columns can be kept in joined dataframe (otherwise clus3D_id column disappears)
            clusters_3d.reset_index(level=("clus3D_id", "hit_id"), names=["event", "clus3D_id", "clus3D_hit_id"]),

            #Right : clusters_2d
            # We don't care about rechits so we slice the df by taking the row of first rechit of 2D cluster
            #                  event     cluster2d_id  index of hit inside cluster2D
            # Also reset index for same reason
            clusters_2d.loc[(slice(None), slice(None),   0                       )].reset_index(level="event"),
            how='inner',

            # Map event on both sides
            # Map clus3D_idxs to clus2D_id
            left_on=('event', 'clus3D_idxs'),
            right_on=('event', 'clus2D_id'),

            suffixes=('', '_clus2D') # This is to avoid beamEnergy column (which exists on both sides) to get renamed. We just keep the one from the left
        )

        #Now merge impact to get impact position of beam on each layer for each event
        clusters_3d_2d_layer = pd.merge(
            # Left : previously merged dataframe
            clusters_3d_2d,

            #Right : impact df (indexed by event and layer)
            impact, 

            # Map event on both sides
            # Map layer of 2D cluster with layer of impact computation
            left_on=("event", "clus2D_layer"),
            right_on=("event", "layer")
        )

        #### Make an index that gets us the highest energy 3D cluster per event
        # Slice clusters3d to remove 2D cluster rows (just take the first row)
        clusters3d_slice = clusters_3d.loc[(slice(None), slice(None), 0)]
        # Build an index that selects for each event the 3D cluster with highest energy
        index_largest_3D_cluster = clusters3d_slice.groupby(["event"])['clus3D_energy'].transform(max) == clusters3d_slice["clus3D_energy"]


        ########## 2D cluster positions wrt incident particle
        #Now we can compute the difference between 2D cluster position and impact of trajectory per layer :
        clusters_3d_2d_layer["clus2D_diff_impact_x"] = clusters_3d_2d_layer["clus2D_x"] - clusters_3d_2d_layer["impactX"]
        clusters_3d_2d_layer["clus2D_diff_impact_y"] = clusters_3d_2d_layer["clus2D_y"] - clusters_3d_2d_layer["impactY"]
        
        hist_dict["clue3d_clusters_spatial_res_count"].fill(
            clusters_3d_2d_layer.beamEnergy, clusters_3d_2d_layer.clus2D_layer, 
            clusters_3d_2d_layer.clus2D_energy, clusters_3d_2d_layer.clus2D_diff_impact_x, clusters_3d_2d_layer.clus2D_diff_impact_y)
        
        ########## TODO fix this :
        #clusters_3d_2d_layer.set_index(["event", "clus3D_id"])
        #clusters_3d_2d_largest_3d_cluster = clusters_3d_2d_layer[index_largest_3D_cluster]
        #hist_dict["clue3d_clusters_spatial_res_largest_3D_cluster"].fill(
        #    clusters_3d_2d_largest_3d_cluster.beamEnergy, clusters_3d_2d_largest_3d_cluster.clus2D_layer, 
        #    clusters_3d_2d_largest_3d_cluster.clus2D_energy, clusters_3d_2d_largest_3d_cluster.clus2D_diff_impact_y, clusters_3d_2d_largest_3d_cluster.clus2D_diff_impact_y)

        # Compute total clustered energy per layer and per 3D cluster
        # Start from df with all 2D clusters per 3D clusters
        # Group by event, 3D cluster and layer
        # Select 2D cluster energy and sum (this will sum over all 2D clusters on the same layer)
        # Also keep beamEnergy
        total_clustered_energy_per_layer_clue3D = (clusters_3d_2d
            .groupby(["event", "clus3D_id", "clus2D_layer"])
            .agg({"clus2D_energy":"sum", "beamEnergy":"first"})
            .reset_index("clus2D_layer")) #so we can use clus2D_layer as a column
        
        hist_dict["clue3d_total_clustered_energy"].fill(total_clustered_energy_per_layer_clue3D.beamEnergy,
            total_clustered_energy_per_layer_clue3D.clus2D_layer, total_clustered_energy_per_layer_clue3D.clus2D_energy)
        
        # Compute total clustered energy per layer for the 3D cluster with highest energy
        # Apply the index to the df with total clustered energy per layer and per 3D cluster, selecting only the correct 3D cluster
        total_clustered_energy_per_layer_highest_clue3D = total_clustered_energy_per_layer_clue3D[index_largest_3D_cluster]
        hist_dict["clue3d_total_clustered_energy_by_largest_cluster"].fill(total_clustered_energy_per_layer_highest_clue3D.beamEnergy,
            total_clustered_energy_per_layer_highest_clue3D.clus2D_layer, total_clustered_energy_per_layer_highest_clue3D.clus2D_energy)

except IndexError as e:
    print("WARNING : an IndexError exception ocurred. This can happen for improperly closed ROOT files.")
    print("WARNING : the last 100MB batch of entries may not have been processed, but the histograms will be written anyway")
    print("The exception was : ")
    print(e)

if hist_dict["rechits_position"].empty():
    print("Result histogram is empty")
    sys.exit(-1)
else:
    os.makedirs(os.path.dirname(args.output_file), exist_ok=True)
    with open(args.output_file, "wb") as f:
        pickle.dump(hist_dict, f)

