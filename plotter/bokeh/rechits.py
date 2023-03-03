import hist

from bokeh.layouts import layout
from bokeh.plotting import curdoc

from utils import *


args = parseArgs()
histStore = HistogramStore(args.hist_folder, args.single_file)

datatype_selector = DatatypeSelector()
clueParamSelector = ClueParamsSelector(histStore.getPossibleClueParameters())
layerSelector = makeLayerSelector()
beamEnergySelector = makeBeamEnergySelector()

h_energy = MultiBokehHistogram(
    HistogramProvider(histStore, "rechits_position", ["beamEnergy", "energy", "layer"], clueParamSelector, datatype_selector),
    "energy", [layerSelector, beamEnergySelector, datatype_selector]
)

h_rec_energy_per_event = MultiBokehHistogram(
    HistogramProvider(histStore, "event_prop", [], clueParamSelector, datatype_selector),
    "totalRecHitEnergy", [beamEnergySelector, datatype_selector]
)

h_xy = MultiBokehHistogram2D(
    HistogramProvider(histStore, "rechits_position", ["beamEnergy", "layer", "x", "y"], clueParamSelector, datatype_selector),
    "x", "y", [layerSelector, beamEnergySelector, datatype_selector]
)

curdoc().add_root(layout([
    [[clueParamSelector.widget, datatype_selector.widget, beamEnergySelector.widget, layerSelector.widget], 
    h_energy.figure,
    h_rec_energy_per_event.figure,
    h_xy.figure]
]))


 