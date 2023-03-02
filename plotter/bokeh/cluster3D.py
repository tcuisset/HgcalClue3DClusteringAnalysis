import hist

from bokeh.layouts import layout
from bokeh.plotting import curdoc

from utils import *


hist_per_datatype = {"data":getHistDict()}

datatype_selector = DatatypeSelector()
layerSelector = makeLayerSelector()
beamEnergySelector = makeBeamEnergySelector()

h_x = MultiBokehHistogram(
    HistogramProvider(hist_per_datatype, "clusters3d_position", ["beamEnergy", "x"], datatype_selector),
    "x", [beamEnergySelector, datatype_selector]
)

h_y = MultiBokehHistogram(
    HistogramProvider(hist_per_datatype, "clusters3d_position", ["beamEnergy", "y"], datatype_selector),
    "y", [beamEnergySelector, datatype_selector]
)

h_xy = MultiBokehHistogram2D(
    HistogramProvider(hist_per_datatype, "clusters3d_position", ["beamEnergy", "x", "y"], datatype_selector),
    "x", "y", [beamEnergySelector, datatype_selector]
)

h_spatial_resolution_xy = MultiBokehHistogram2D(
    HistogramProvider(hist_per_datatype, "clue3d_clusters_spatial_res", ["beamEnergy", "layer", "x", "y", "cluster2dEnergy"], datatype_selector),
    "x", "y", [beamEnergySelector, datatype_selector, layerSelector]
)

curdoc().add_root(layout([
    [[datatype_selector.widget, beamEnergySelector.widget, layerSelector.widget], h_x.figure, h_y.figure, h_xy.figure, h_spatial_resolution_xy.figure]
]))


 