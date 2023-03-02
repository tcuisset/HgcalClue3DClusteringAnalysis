import hist

from bokeh.layouts import layout
from bokeh.plotting import curdoc

from utils import *


hist_per_datatype = {"data":getHistDict()}

datatype_selector = DatatypeSelector()
layerSelector = makeLayerSelector()
beamEnergySelector = makeBeamEnergySelector()

h_energy = MultiBokehHistogram(
    HistogramProvider(hist_per_datatype, "rechits_position", ["beamEnergy", "energy", "layer"], datatype_selector),
    "energy", [layerSelector, beamEnergySelector, datatype_selector]
)

h_rec_energy_per_event = MultiBokehHistogram(
    HistogramProvider(hist_per_datatype, "event_prop", [], datatype_selector),
    "totalRecHitEnergy", [beamEnergySelector, datatype_selector]
)

h_xy = MultiBokehHistogram2D(
    HistogramProvider(hist_per_datatype, "rechits_position", ["beamEnergy", "layer", "x", "y"], datatype_selector),
    "x", "y", [layerSelector, beamEnergySelector, datatype_selector]
)

curdoc().add_root(layout([
    [[datatype_selector.widget, beamEnergySelector.widget, layerSelector.widget], 
    h_energy.figure,
    h_rec_energy_per_event.figure,
    h_xy.figure]
]))


 