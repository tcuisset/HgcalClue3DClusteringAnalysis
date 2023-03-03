import pickle
import os
import argparse

import numpy as np
import bokeh.plotting
import bokeh.models
from bokeh.models import ColumnDataSource
import hist
import boost_histogram as bh

outputCacheHistsFolder = '/home/llr/cms/cuisset/hgcal/testbeam18/clue3d-dev/src/plots/cache'
beamEnergies = [20, 50, 80, 100, 120, 150, 200, 250, 300]
datatypes = ['data', 'sim_proton', 'sim_noproton']

def parseArgs():
    parser = argparse.ArgumentParser(description="Plotting code to be run using Bokeh server, use bokeh serve SCRIPT.py --args ARGS")
    parser.add_argument("--hist-folder", dest="hist_folder",
        help="path to folder holding all histograms. Will load recursively all clueparams and datatypes inside this folder")
    parser.add_argument("--single-file", dest='single_file',
        help="Only load a single pickle file (for debugging), given by this full path to the file")
    return parser.parse_args()


class HistogramStore:
    # {
    #     ("default", "sim_proton") : {
    #         "rechits_position" : hist.Hist object,
    #         more hists ...
    #     },
    #     ("other_clue3D_param", "data") : {
    #         hists...
    #     }
    # }
    hist_dict = {} # dict of dicts. 
    def __init__(self, hist_folder=None, single_file=None) -> None:
        if hist_folder is not None:
            self._load_recursive(hist_folder)
        elif single_file is not None:
            self._load_single(single_file)
        else:
            raise ValueError("Need to specify either --hist-folder or --single-file argument to the script, use bokeh serve --args --histfolder=...")

    def _load_single(self, single_file):
        with open(os.path.join(single_file, 'hists.pkl'), "rb") as f:
            self.hist_dict[(None, None)] = pickle.load(f)

    def _load_recursive(self, hist_folder):
        (clue_param_path, clue_param_names, _) = next(os.walk(hist_folder))
        for clue_param_name in clue_param_names:
            (datatype_path, datatype_names, filenames) = next(os.walk(os.path.join(clue_param_path, clue_param_name)))
            for datatype_name in datatype_names:
                try:
                    with open(os.path.join(datatype_path, datatype_name, 'hists.pkl'), "rb") as f:
                        self.hist_dict[(clue_param_name, datatype_name)] = pickle.load(f)
                except FileNotFoundError:
                    pass
                    
    def getHistogram(self, clue_param_name, datatype, hist_name) -> hist.Hist:
        return self.hist_dict[(clue_param_name, datatype)][hist_name]
    
    def getPossibleParameters(self):
        return self.hist_dict.keys()
    def getPossibleClueParameters(self):
        return list(set([t[0] for t in self.hist_dict.keys()]))




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
    
    def registerCallback(self, callback):
        #Use value_throttled so that it updates only on mouse release to avoid recomputing all the time when dragging
        self.widget.on_change('value_throttled', callback)

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

class ClueParamsSelector(AxisSelector):
    def __init__(self, clueParamList) -> None:
        self.widget = bokeh.models.RadioButtonGroup(
            name="clue_params",
            labels=clueParamList,
            active=0
        )
    
    def getClueParam(self):
        return self.widget.labels[self.widget.active]
    def registerCallback(self, callback):
        self.widget.on_change('active', callback)

class PlaceholderClueParamsSelector(AxisSelector):
    def getClueParam(self):
        return "default"
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
    def __init__(self, histStore:HistogramStore, hist_name, projections, clue3DSelector:ClueParamsSelector, datatypeSelector:DatatypeSelector,
        profileOn=None) -> None:
        self.histStore = histStore
        self.hist_name = hist_name
        self.projections = projections
        self.clue3DSelector = clue3DSelector
        self.datatypeSelector = datatypeSelector
        self.profileOn = profileOn

    def getHist(self):
        h = self.histStore.getHistogram(
            self.clue3DSelector.getClueParam(),
            self.datatypeSelector.getDatatype(),
            self.hist_name
        )
        h_projected = h
        if self.projections:
            # Need to expand self.projections which is a list into arguments for project 
            # https://docs.python.org/3/tutorial/controlflow.html#unpacking-argument-lists
            h_projected = h.project(*self.projections)
        
        if self.profileOn is not None:
            return h_projected.profile(self.profileOn)
        else:
            return h_projected


class MultiBokehHistogram:
    def __init__(self, histProvider:HistogramProvider, plotAxis:str, axisSelectors, **kwargs) -> None:
        self.figure = bokeh.plotting.figure(**kwargs)

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
    def __init__(self, histProvider:HistogramProvider, xAxis:str, yAxis:str, axisSelectors, **kwargs) -> None:
        self.figure = bokeh.plotting.figure(**kwargs)

        self.axisSelectors = axisSelectors
        self.histProvider = histProvider
        self.xAxis = xAxis
        self.yAxis = yAxis

        for axisSelector in axisSelectors:
            axisSelector.registerCallback(self.update)

        h_proj = self.getHistogramProjection()
        self.figure.xaxis.axis_label = h_proj.axes[xAxis].label
        self.figure.yaxis.axis_label = h_proj.axes[yAxis].label

        self.source = ColumnDataSource()
        self.update(0, 0, 0)
        self.figure.image(image='histogram_2D_view', x=h_proj.axes[xAxis].edges[0], y=h_proj.axes[yAxis].edges[0],
            dw=h_proj.axes[xAxis].edges[-1]-h_proj.axes[xAxis].edges[0], dh=h_proj.axes[yAxis].edges[-1]-h_proj.axes[yAxis].edges[0],
            palette="Spectral11",
            source=self.source)
    
    def getAllSlices(self):
        slice_dict = {}
        for axisSelector in self.axisSelectors:
            slice = axisSelector.getSlice()
            if slice != ():
                slice_dict[slice[0]] = slice[1]
        return slice_dict
    
    def getHistogramProjection(self) -> hist.Hist:
        #We do an additional projection on plotAxis as required by MultiSelectAxisSelector (see MultiSelectAxisSelector.getSlice)
        return self.histProvider.getHist()[self.getAllSlices()].project(self.xAxis, self.yAxis)
    
    def update(self, attr, old, new):
        h_proj = self.getHistogramProjection()
        
        if h_proj.kind is bh.Kind.COUNT:
            # regular histogram
            view = h_proj.view()
        elif h_proj.kind is bh.Kind.MEAN:
            # profile histogram : we want the mean
            view = h_proj.view().value
        else:
            raise ValueError("Histogram kind not supported : " + str(h_proj.kind) + ". Only COUNT and MEAN histograms are supported.")

        # It would seem that the histogram x, y view is the transpose of what is expected by bokeh, though it needs to be checked
        self.source.data = {"histogram_2D_view":[np.transpose(view)]}