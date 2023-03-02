import os
import pickle

import uproot
import awkward as ak
import numpy as np
import boost_histogram as bh
import hist
import pandas as pd

fname='ClusteringAnalysis/CLUE_clusters_single.root' ## with CLUS->LC->rechit mapping
#outPlotFolder = '/eos/user/t/tcuisset/www/testbeam18/clue3d'
outPlotFolder = '/home/llr/cms/cuisset/hgcal/testbeam18/clue3d-dev/src/plots'
outputCacheHistsFolder = '/home/llr/cms/cuisset/hgcal/testbeam18/clue3d-dev/src/plots/cache'


beamEnergies = [20, 50, 80, 100, 120, 150, 200, 250, 300]
beamEnergiesAxis = hist.axis.IntCategory(beamEnergies, name="beamEnergy", label="Beam energy (GeV)")
layerAxis = hist.axis.Integer(start=0, stop=30, name="layer", label="Layer number")

rho_axis = hist.axis.Regular(bins=100, start=0, stop=10., name="rho", label="Energy density")
delta_axis = hist.axis.Regular(bins=100, start=0, stop=3., name="delta", label="Distance to nearest higher")
seed_axis = hist.axis.IntCategory([0, 1], name="isSeed") #0: not a seed, 1: is a seed

cluster2dEnergy_axis = hist.axis.Regular(100, 0, 1, name="cluster2dEnergy", label="2D cluster energy")

cluster3dSizeAxis = hist.axis.Integer(1, 40, name="cluster3dSize", label="3D cluster size")
cluster3dEnergyAxis = hist.axis.Regular(100, 0, 20, name="cluster3dEnergy", label="3D cluster energy")

xPosition_axis = hist.axis.Regular(bins=10, start=-10., stop=10., name="x", label="x")
yPosition_axis = hist.axis.Regular(bins=10, start=-10., stop=10., name="y", label="y")
zPosition_axis = hist.axis.Regular(bins=10, start=0., stop=60., name="z", label="z")
energy_axis = hist.axis.Regular(bins=100, start=0., stop=10., name="energy", label="Reconstructed hit energy")

totalRecHitEnergy_axis = hist.axis.Regular(bins=500, start=0, stop=300, name="totalRecHitEnergy", label="Total RecHit energy per event (GeV)")

hist_dict = {}

hist_dict["rechits_position"] = hist.Hist(beamEnergiesAxis, layerAxis, energy_axis, xPosition_axis, yPosition_axis)
hist_dict["rechits_prop"] = hist.Hist(beamEnergiesAxis, layerAxis, seed_axis, rho_axis, delta_axis)

hist_dict["clusters2d_position"] = hist.Hist(beamEnergiesAxis, layerAxis, seed_axis, xPosition_axis, yPosition_axis)
hist_dict["clusters2d_prop"] = hist.Hist(beamEnergiesAxis, layerAxis, seed_axis, rho_axis, delta_axis)

hist_dict["clusters3d_position"] = hist.Hist(beamEnergiesAxis, cluster3dEnergyAxis, xPosition_axis, yPosition_axis, zPosition_axis)
hist_dict["clusters3d_prop"] = hist.Hist(beamEnergiesAxis, cluster3dEnergyAxis, cluster3dSizeAxis)

hist_dict["event_prop"] = hist.Hist(beamEnergiesAxis, totalRecHitEnergy_axis)
hist_dict["clue3d_clusters_spatial_res"] = hist.Hist(beamEnergiesAxis, layerAxis, cluster2dEnergy_axis, xPosition_axis, yPosition_axis)

for array in uproot.iterate(fname + ":clusters", step_size="100MB", library="ak"):
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
    hist_dict["clusters3d_position"].fill(clusters3D_slice.beamEnergy, clusters3D_slice.clus3D_energy, clusters3D_slice.clus3D_x, clusters3D_slice.clus3D_y, clusters3D_slice.clus3D_z)


    #Merge clusters3D and clusters2D
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

    hist_dict["clue3d_clusters_spatial_res"].fill(clusters_3d_2d_layer.beamEnergy, clusters_3d_2d_layer.clus2D_layer, 
        clusters_3d_2d_layer.clus2D_energy, clusters_3d_2d_layer.clus2D_x, clusters_3d_2d_layer.clus2D_y)

os.makedirs(outputCacheHistsFolder, exist_ok=True)
with open(os.path.join(outputCacheHistsFolder, "hists.pkl"), "wb") as f:
    pickle.dump(hist_dict, f)

