from simulator import *
import numpy as np
import multiprocessing  as mp
import os, copy
from Graph import QWGraph as qgw
from Graph import *
import scipy.stats as stats


def unit_list(bounds, unit):

    b_0 = max(0, bounds[0]// unit.distance())

    return np.arange(b_0, bounds[1]// unit.distance())

class QWGraphCollection(object) :

    def __init__( self, analyzer : Analyzer = None, name = None ) :

        if analyzer == None:
            analyzer = Analyzer()
        self._analyzer = analyzer
        self._name = name

        self._gr_list : list[QWGraph] = []

    def add(self, gr : QWGraph) :
        self._gr_list.append( gr )

    def get_data( self, target = "p", diag = True, opt_mode = None, opt_phi = None, x_mode = "dist"):

        graphs = self.get_list()
        tester = self.get_analyzer()
        N = len(graphs)
        data = np.empty( N)


        prog_label = self.get_name() +" collection " + target + " data:  {:2.1%}"

        def evaluate( tester):

            best_phi = 0

            if opt_mode == "smart" :
                best_phi = tester.optimum_phase_smart()[0]
            elif opt_mode == "fix" :
                assert opt_phi != None
                best_phi = opt_phi
            else :
                best_phi = tester.optimum_phase_minimize(diag = diag)[0]

            #print(best_phi, opt_mode)
            target_t = (target != "p")

            #todo : t is boolean???????
            if diag:
                return tester.performance_diag(best_phi, t = target_t)
            else:
                return tester.performance(best_phi, t = target_t)

        for i in range(N):
            print (prog_label.format(1/N), end = "\r")
            tester.set_gr(graphs[i])

            data[i] = evaluate(tester)

        print(prog_label.format(1))

        # n_proc = os.cpu_count()*2  
        # n_proc = 1 
        # print("Starting pool evaluation with {} process".format(n_proc))
        # with mp.Pool( n_proc) as pool:
        #     print( pool.map(evaluate, testers))

        #     pool.close()
        #     pool.join()

        x = self.get_x_vec(x_mode = x_mode)

        return x , data

    def get_x_vec(self, x_mode = "dist"):
        out = []

        for gr in self.get_list():
            if x_mode == "dist":
                out.append( gr.distance())
            elif x_mode == "size" :
                out.append( gr.N)
            else :
                raise ValueError("get_x_vec mode not supported")
        
        return out


    #todo make it return lm object and not just the two coefficients
    def transport_time_lm( self, mode = "dist"):
        """
        Construct a linear model of the best transport time from the graph collection
        """
        
        x, data = self.get_data( target = "t")

        #print(data)
        out = stats.linregress(x, data)

        print("m: ", out.slope, " +- ", out.stderr)
        print("q: ", out.intercept, " +- ", out.intercept_stderr)
        print("r: ", out.rvalue)

        return out.slope, out.intercept
        
        
    def get_list(self):

        return self._gr_list

    def set_analyzer( self, aAnalyzer):
        self._analyzer = aAnalyzer

    def get_analyzer( self):
        return self._analyzer

    def get_name( self):
        return self.get_list()[0].code if self._name == None else self._name

class CollectionBuilder(object) :

    def from_list( gr_list : list[QWGraph], analyzer : Analyzer = None) -> QWGraphCollection :

        cb = QWGraphCollection( analyzer=analyzer)

        for gr in gr_list :
            cb.add( gr )

        return cb

    def P_progression( bounds = None, step = 1, select : list[int] = None, analyzer : Analyzer = None) :

        cb = QWGraphCollection( analyzer=analyzer)

        assert bounds or select

        if select != None:
            drange = select
        else :
            drange = np.arange(bounds[0], bounds[1], step)

        for d in drange :
            cb.add( qgw.Line(d))

        return cb

    def C_progression( bounds = None, step = 1, select : list[int] = None, analyzer : Analyzer = None, **kwargs) :

        cb = QWGraphCollection( analyzer=analyzer)

        assert bounds or select

        if select != None:
            drange = select
        else :
            drange = np.arange(bounds[0], bounds[1], step)

        for d in drange :
            cb.add( qgw.Ring(d, **kwargs))

        return cb

    def base_progression( g_type, **kwargs):
        """
        wrapper for the 3 basic standard progression
        """

        if g_type == "P":
            return CollectionBuilder.P_progression(**kwargs)
        if g_type == "C":
            return CollectionBuilder.C_progression(**kwargs)
        if g_type == "Ch":
            return CollectionBuilder.C_progression(HANDLES = True, **kwargs)

    def chain_progression( gr_unit, bounds = None, step = 1, select : list[int] = None, analyzer: Analyzer = None, **kwargs) :

        cb = QWGraphCollection( analyzer=analyzer)

        assert bounds or select

        if select != None:
            drange = select
        else :
            drange = unit_list( bounds, gr_unit)

        for d in drange :
            cb.add( qgw.chain( gr_unit, rep = d, **kwargs))

        return cb


#####################################à
#utils methods

def get_line_data(bounds = (3,10), target = "p", x_mode = "dist", **kwargs):
    """
    Simple wrapper for L graphs progression references   
    """

    an = Analyzer( mode = "first")
    line_collection = CollectionBuilder.P_progression( bounds, analyzer=an)

    return line_collection.get_data( target = target, x_mode = x_mode)