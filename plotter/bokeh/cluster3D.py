import hist

from bokeh.layouts import layout, column
from bokeh.plotting import curdoc

from utils import *

args = parseArgs()
histStore = HistogramStore(args.hist_folder, args.single_file)

datatype_selector = DatatypeSelector()
clueParamSelector = ClueParamsSelector(histStore.getPossibleClueParameters())
layerSelector = makeLayerSelector()
beamEnergySelector = makeBeamEnergySelector()

h_x = MultiBokehHistogram(
    HistogramProvider(histStore, "clusters3d_position_xy", ["beamEnergy", "x"], clueParamSelector, datatype_selector),
    "x", [beamEnergySelector, datatype_selector, clueParamSelector],
    title="3D cluster x position (cluster count)"
)

h_y = MultiBokehHistogram(
    HistogramProvider(histStore, "clusters3d_position_xy", ["beamEnergy", "y"], clueParamSelector, datatype_selector),
    "y", [beamEnergySelector, datatype_selector, clueParamSelector],
    title="3D cluster y position (cluster count)"
)

h_xy = MultiBokehHistogram2D(
    HistogramProvider(histStore, "clusters3d_position_xy", ["beamEnergy", "x", "y"], clueParamSelector, datatype_selector),
    "x", "y", [beamEnergySelector, datatype_selector, clueParamSelector],
    title="3D cluster x-y position (cluster count)"
)

h_spatial_resolution_xy = MultiBokehHistogram2D(
    HistogramProvider(histStore, "clue3d_clusters_spatial_res_count", ["beamEnergy", "layer", "diffX", "diffY"], clueParamSelector, datatype_selector),
    "diffX", "diffY", [beamEnergySelector, datatype_selector, layerSelector, clueParamSelector],
    title="2D cluster resolution of all 2D clusters clustered inside a 3D cluster (clustered energy is projected away)"
)

h_spatial_resolution_xy_profiled = MultiBokehHistogram2D(
    HistogramProvider(histStore, "clue3d_clusters_spatial_res_count", ["beamEnergy", "layer", "diffX", "diffY", "cluster2dEnergy"],
        clueParamSelector, datatype_selector, profileOn="cluster2dEnergy"),
    "diffX", "diffY", [beamEnergySelector, datatype_selector, layerSelector, clueParamSelector],
    title="2D cluster resolution of all 2D clusters clustered inside a 3D cluster (projected on 2D cluster energy)"
)

curdoc().add_root(layout([
    [column(clueParamSelector.widget, datatype_selector.widget, beamEnergySelector.widget, layerSelector.widget),
    h_x.figure, h_y.figure],
    [h_xy.figure, h_spatial_resolution_xy.figure, h_spatial_resolution_xy_profiled.figure]
]))


 