import pickle
import os

import bokeh.plotting
import bokeh.models
from bokeh.models import ColumnDataSource
import hist

outputCacheHistsFolder = '/home/llr/cms/cuisset/hgcal/testbeam18/clue3d-dev/src/plots/cache'
beamEnergies = [20, 50, 80, 100, 120, 150, 200, 250, 300]


def getHistDict():
    with open(os.path.join(outputCacheHistsFolder, "hists.pkl"), "rb") as f:
        hist_dict = pickle.load(f)
    return hist_dict


class AxisSelector:
    def getSlice(self):
        return ()
    def registerCallback(self, callback):
        self.widget.on_change('value', callback)

class RangeAxisSelector(AxisSelector):
    def __init__(self, axisName:str, **kwargs) -> None:
        self.axisName = axisName
        self.widget = bokeh.models.RangeSlider(**kwargs)

    def getSlice(self):
        (min, max) = self.widget.value
        return (self.axisName, slice(hist.loc(min), hist.loc(max+1), sum))

class MultiSelectAxisSelector(AxisSelector):
    def __init__(self, axisName:str, **kwargs) -> None:
        self.axisName = axisName
        self.widget = bokeh.models.MultiSelect(**kwargs)
    
    def getSlice(self):
        # !!! does not project (needs to be projected)
        # for now boost-histogram using list slicing does not fill overflow bins
        # we need this behaviour for the project call in MultiBokehHistogram.getHistogramProjection to return the correct projection on selected beam energies
        # when boos-histogram is updated you should write probably [hist.loc(int(energy)) for energy in self.widget.value]:sum or something
        # see  https://github.com/scikit-hep/boost-histogram/issues/296
        return (self.axisName, [hist.loc(int(energy)) for energy in self.widget.value])

class PlaceholderAxisSelector(AxisSelector):
    def __init__(self, slice) -> None:
        self.slice = slice
    def getSlice(self):
        return self.slice
    def registerCallback(self, callback):
        pass

class DatatypeSelector(AxisSelector):
    def __init__(self) -> None:
        self.labels = ["data", "sim_proton", "sim_noproton"]
        self.widget = bokeh.models.RadioButtonGroup(
            name="datatype",
            labels=self.labels,
            active=0
        )

    def getDatatype(self):
        # self.radio.value is the button number that is pressed -> map it to label
        return self.labels[self.widget.active]
    
    def registerCallback(self, callback):
        self.widget.on_change('active', callback)

class PlaceholderDatatypeSelector(AxisSelector):
    def getDatatype(self):
        return "data"
    def registerCallback(self, callback):
        pass

def makeLayerSelector():
    return RangeAxisSelector("layer",
        title="Layer selection",
        start=0,
        end=30,
        step=1,
        value=(1,28)
    )

def makeBeamEnergySelector():
    return MultiSelectAxisSelector("beamEnergy",
        title="Beam energy",
        options=[(str(value), str(value) + " GeV") for value in beamEnergies],
        value=[str(energy) for energy in beamEnergies],
        height=200
    )

class HistogramProvider:
    def __init__(self, hist_per_datatype, hist_name, projections, datatypeSelector:DatatypeSelector) -> None:
        self.hist_per_datatype = hist_per_datatype
        self.hist_name = hist_name
        self.projections = projections
        self.datatypeSelector = datatypeSelector

    def getHist(self):
        h = self.hist_per_datatype[self.datatypeSelector.getDatatype()][self.hist_name]

        if self.projections:
            # Need to expand self.projections which is a list into arguments for project 
            # https://docs.python.org/3/tutorial/controlflow.html#unpacking-argument-lists
            return h.project(*self.projections)
        else:
            return h


class MultiBokehHistogram:
    def __init__(self, histProvider:HistogramProvider, plotAxis:str, axisSelectors) -> None:
        self.figure = bokeh.plotting.figure()

        self.axisSelectors = axisSelectors
        self.histProvider = histProvider
        self.plotAxis = plotAxis

        for axisSelector in axisSelectors:
            axisSelector.registerCallback(self.update)

        self.figure.xaxis.axis_label = self.getHistogramProjection().axes[0].label
        self.figure.yaxis.axis_label = "Count"

        self.source = ColumnDataSource()
        self.update(0, 0, 0)
        self.figure.quad(bottom=0, source=self.source)
    
    def getAllSlices(self):
        slice_dict = {}
        for axisSelector in self.axisSelectors:
            slice = axisSelector.getSlice()
            if slice != ():
                slice_dict[slice[0]] = slice[1]
        return slice_dict
    
    def getHistogramProjection(self):
        #We do an additional projection on plotAxis as required by MultiSelectAxisSelector (see MultiSelectAxisSelector.getSlice)
        return self.histProvider.getHist()[self.getAllSlices()].project(self.plotAxis)
    
    def update(self, attr, old, new):
        h_proj = self.getHistogramProjection()
        self.source.data = {"top":h_proj.view(), "left":h_proj.axes[0].edges[:-1], "right":h_proj.axes[0].edges[1:]}


class MultiBokehHistogram2D:
    def __init__(self, histProvider:HistogramProvider, xAxis:str, yAxis:str, axisSelectors) -> None:
        self.figure = bokeh.plotting.figure()

        self.axisSelectors = axisSelectors
        self.histProvider = histProvider
        self.xAxis = xAxis
        self.yAxis = yAxis

        for axisSelector in axisSelectors:
            axisSelector.registerCallback(self.update)

        h_proj = self.getHistogramProjection()
        self.figure.xaxis.axis_label = h_proj.axes[0].label
        self.figure.yaxis.axis_label = h_proj.axes[1].label

        self.source = ColumnDataSource()
        self.update(0, 0, 0)
        self.figure.image(image='histogram_2D_view', x=h_proj.axes[0].edges[0], y=h_proj.axes[1].edges[0],
            dw=h_proj.axes[0].edges[-1]-h_proj.axes[0].edges[0], dh=h_proj.axes[1].edges[-1]-h_proj.axes[1].edges[0],
            palette="Spectral11",
            source=self.source)
    
    def getAllSlices(self):
        slice_dict = {}
        for axisSelector in self.axisSelectors:
            slice = axisSelector.getSlice()
            if slice != ():
                slice_dict[slice[0]] = slice[1]
        return slice_dict
    
    def getHistogramProjection(self):
        #We do an additional projection on plotAxis as required by MultiSelectAxisSelector (see MultiSelectAxisSelector.getSlice)
        return self.histProvider.getHist()[self.getAllSlices()].project(self.xAxis, self.yAxis)
    
    def update(self, attr, old, new):
        h_proj = self.getHistogramProjection()
        self.source.data = {"histogram_2D_view":[h_proj.view()]}