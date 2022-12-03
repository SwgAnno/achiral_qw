from simulator import *
import numpy as np
import multiprocessing as mp
from Graph import QWGraph as qgw
from Graph import *

def unit_list(bounds, unit):

    b_0 = max(0, bounds[0]// unit.distance())

    return np.arange(b_0, bounds[1]// unit.distance())

class QWGraphCollection(object) :

    def __init__( self, analyzer : Analyzer = None ) :

        if analyzer == None:
            self._analyzer = Analyzer()

        self._gr_list : list[QWGraph] = []

    def add(self, gr : QWGraph) :
        self._gr_list.append( gr )

    def set_analyzer( self, aAnalyzer):
        self._analyzer = aAnalyzer

    def get_analyzer( self):
        return self._analyzer

class CollectionBuilder(object) :

    def from_list( gr_list : list[QWGraph], analyzer : Analyzer = None) -> QWGraphCollection :

        cb = QWGraphCollection( analyzer=analyzer)

        for gr in gr_list :
            cb.add( gr )

        return cb

    def P_progression( bounds = None, select : list[int] = None, analyzer : Analyzer = None) :

        cb = QWGraphCollection( analyzer=analyzer)

        assert bounds or select

        if select != None:
            drange = select
        else :
            drange = np.arange(bounds[0], bounds[1], bounds[2])

        for d in drange :
            cb.add( qgw.Line(d))

        return cb

    def C_progression( bounds = None, select : list[int] = None, analyzer : Analyzer = None, **kwargs) :

        cb = QWGraphCollection( analyzer=analyzer)

        assert bounds or select

        if select != None:
            drange = select
        else :
            drange = np.arange(bounds[0], bounds[1], bounds[2])

        for d in drange :
            cb.add( qgw.Ring(d, **kwargs))

        return cb

    def chain_progression( gr_unit, bounds = None, select : list[int] = None, analyzer: Analyzer = None, **kwargs) :

        cb = QWGraphCollection( analyzer=analyzer)

        assert bounds or select

        if select != None:
            drange = select
        else :
            drange = unit_list( bounds, gr_unit)

        for d in drange :
            cb.add( qgw.chain( gr_unit, rep = d, **kwargs))

        return cb