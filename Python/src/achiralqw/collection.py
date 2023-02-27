from achiralqw.simulator import Analyzer
import numpy as np
import achiralqw.istarmap
import multiprocessing  as mp
import os, copy
from itertools import repeat
from achiralqw.graph import QWGraph, QWGraphBuilder
import scipy.stats as stats
from scipy.optimize import curve_fit
import tqdm

MULTIPROCESS = True

def get_gr_t_data( an : Analyzer,):
    return an.evaluate(target="t")

def get_gr_p_data( an : Analyzer,):
    return an.evaluate(target = "p")

def set_graph( an : Analyzer, graph : QWGraph):

    an = copy.deepcopy(an)
    an.set_gr(graph)
    return an

def create_c(dist):
    return QWGraphBuilder.Ring(dist)

def create_ch(dist):
    return QWGraphBuilder.Ring(dist, HANDLES = True)

def unit_list_bounds(bounds, unit):

    b_0 = max(0, bounds[0]// unit.distance())

    return np.arange(b_0, bounds[1]// unit.distance())

def unit_list( select, unit):

    num_unit = (select -2)// unit.distance()

    return num_unit

class QWGraphCollection(object) :

    def __init__( self, analyzer : Analyzer = None, name = None ) :

        if analyzer == None:
            analyzer = Analyzer()
        self._analyzer = analyzer
        self._name = name

        global MULTIPROCESS

        if MULTIPROCESS:
            self.get_data = self.get_data_multiprocess
        else :
            self.get_data = self.get_data_singleprocess

        self._gr_list : list[QWGraph] = []

    def add(self, gr : QWGraph) :
        self._gr_list.append( gr )

    def get_data_singleprocess( self, target = "p", diag = True, opt_mode = None, opt_phi = None, x_mode = "dist"):

        graphs = self.get_list()
        tester = self.get_analyzer()
        N = len(graphs)
        data = np.empty( N)

        prog_label = self.get_name() +" collection " + target + " data:  {:2.1%}"

        for i in range(N):
            print (prog_label.format(i/N), end = "\r")
            tester.set_gr(graphs[i])

            data[i] = tester.evaluate(target=target)

        print(prog_label.format(1))

        x = self.get_x_vec(x_mode = x_mode)

        return x , data

    def get_data_multiprocess( self, target = "p", x_mode = "dist"):

        graphs = self.get_list()
        tester = self.get_analyzer()
        N = len(graphs)
        out = []

        global get_gr_t_data
        global get_gr_p_data
        global set_graph

        if target == "t":
            get_gr_data = get_gr_t_data
        else:
            get_gr_data = get_gr_p_data

        n_proc = os.cpu_count()*2  

        greeting_string = self.get_name() +" collection: Starting pool evaluation with {} process" 
        print(greeting_string.format(n_proc))
        with mp.Pool( n_proc) as pool:

            input_vec = zip([tester]* len(graphs), graphs)
            testers = []

            print("Data Setup")
            for _ in tqdm.tqdm(pool.istarmap(set_graph, input_vec), total=len(testers)):
                testers.append(_)

            print("Evaluation")
            for _ in tqdm.tqdm(pool.imap(get_gr_data, testers), total=len(testers)):
                out.append(_)
            data = np.array(out)

            pool.close()
            pool.join()

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

    def transport_prob_model( self, mode = "poly"):
        """
        Very specific analysis tool: it extracts the slope of the probability progression on a loglog scale,
        which represent the power law exponent of its decay
        """
        x, data = self.get_data( target = "p")

        #test over y = ax^2 + bx + c
        if mode == "poly" :
            param = np.polyfit(np.log(x), np.log(data), deg = 2)
            
            print("ax^2: ", param[0])
            print("bx: ", param[1])
            print("c: ", param[2])

            return param

        elif mode == "custom": 
            # test over y = ax + b + c*1/x   
            def model(x,a,b,c):
                return a*x + b + c/x

            #first guess
            p0 = [-.5,2.5,1]
            param, cov_param = curve_fit(model, np.log(x), np.log(data), p0 = p0)

            print("ax: ", param[0], " +- ", cov_param[0,0])
            print("b: ", param[1], " +- ", cov_param[1,1])
            print("c*1/x: ", param[2], " +- ", cov_param[2,2])

            return param


        return None

    def transport_prob_loglog_lm( self, mode = "dist"):
        """
        Very specific analysis tool: it extracts the slope of the probability progression on a loglog scale,
        which represent the power law exponent of its decay
        """
        
        x, data = self.get_data( target = "p")

        #print(data)
        out = stats.linregress(np.log(x), np.log(data))

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

    def __init__(self):

        global MULTIPROCESS

        if MULTIPROCESS:
            self.chain_progression = self.chain_progression_multiprocess
            self.P_progression = self.P_progression_multiprocess
            self.C_progression = self.C_progression_multiprocess
        else :
            self.chain_progression = self.chain_progression_singleprocess
            self.P_progression = self.P_progression_singleprocess
            self.C_progression = self.C_progression_singleprocess


    def from_list(self, gr_list , analyzer : Analyzer = None) -> QWGraphCollection :

        collection = QWGraphCollection( analyzer=analyzer)

        for gr in gr_list:
            collection.add( gr )

        return collection

    def P_progression_singleprocess(self, bounds = None, step = 1, select = None, analyzer : Analyzer = None) :

        collection = QWGraphCollection( analyzer=analyzer)
        collection.get_analyzer().set_opt_mode("none")

        assert bounds or np.any(select)

        if np.any(select):
            drange = select
        else :
            drange = np.arange(bounds[0], bounds[1], step)

        for d in drange :
            collection.add( QWGraphBuilder.Line(d))

        return collection

    def P_progression_multiprocess(self, bounds = None, step = 1, select = None, analyzer : Analyzer = None, **kwargs) :

        collection = QWGraphCollection( analyzer=analyzer)
        collection.get_analyzer().set_opt_mode("none")

        assert bounds or np.any(select)

        if np.any(select):
            drange = select
        else :
            drange = np.arange(bounds[0], bounds[1], step)

        n_proc = os.cpu_count()*2  
        greeting_string = "P progression: Starting pool creation with {} process" 
        print(greeting_string.format(n_proc))

        with mp.Pool( n_proc) as pool:
            for _ in tqdm.tqdm(pool.imap(QWGraphBuilder.Line, drange ), total=len(drange)):
                collection.add(_)

            pool.close()
            pool.join()

        return collection



    def C_progression_singleprocess(self, bounds = None, step = 1, odd = False, select = None, analyzer : Analyzer = None, HANDLES = False, **kwargs) :

        collection = QWGraphCollection( analyzer=analyzer)

        assert bounds or np.any(select)

        #start with an odd cycle
        offset = 1 if odd else 0
        
        #account fo extra distance
        if HANDLES:
            offset -= 4 


        if np.any(select):
            drange = select
        else :
            drange = np.arange(bounds[0]*2 + offset, bounds[1]*2 + offset, step)

        for d in drange :
            collection.add( QWGraphBuilder.Ring(d, HANDLES = HANDLES, **kwargs))

        return collection

    def C_progression_multiprocess(self, bounds = None, step = 1, odd = False, select = None, analyzer : Analyzer = None, HANDLES = False, 
     **kwargs) :

        collection = QWGraphCollection( analyzer=analyzer)

        assert bounds or np.any(select)
        
        #start with an odd cycle
        offset = 1 if odd else 0

        #account fo extra distance
        if HANDLES:
            offset -= 4 


        if np.any(select):
            drange = select
        else :
            drange = np.arange(bounds[0]*2 + offset, bounds[1]*2 + offset, step)

        n_proc = os.cpu_count()*2  
        greeting_string = "C progression: Starting pool creation with {} process" 
        print(greeting_string.format(n_proc))

        with mp.Pool( n_proc) as pool:
            if HANDLES :
                for _ in tqdm.tqdm(pool.imap(create_ch, drange ), total=len(drange)):
                    collection.add(_)
            else :
                for _ in tqdm.tqdm(pool.imap(create_c, drange ), total=len(drange)):
                    collection.add(_)

            pool.close()
            pool.join()

        return collection

    def chain_progression_singleprocess(self, gr_unit, bounds = None, step = 1, select = None, analyzer: Analyzer = None, **kwargs) :

        collection = QWGraphCollection( analyzer=analyzer)

        assert bounds or np.any(select)

        if np.any(select):
            drange = unit_list(select, gr_unit)
        else :
            drange = unit_list_bounds( bounds, gr_unit)

        for d in drange :
            collection.add( gr_unit.chain( rep = d, **kwargs))

        return collection

    def chain_progression_multiprocess(self, gr_unit, bounds = None, step = 1, select = None, analyzer = None, **kwargs):
        collection = QWGraphCollection( analyzer=analyzer)

        assert bounds or np.any(select)

        if np.any(select):
            drange = unit_list(select, gr_unit)
        else :
            drange = unit_list_bounds( bounds, gr_unit)


        n_proc = os.cpu_count()*2  
        greeting_string = gr_unit.code +" chain progression: Starting pool creation with {} process" 
        print(greeting_string.format(n_proc))

        input_vec = [(gr_unit, rep) for rep in drange]

        with mp.Pool( n_proc) as pool:
            for _ in tqdm.tqdm(pool.istarmap(QWGraph.chain, input_vec ), total=len(drange)):
                collection.add(_)

            pool.close()
            pool.join()

        return collection

    def base_progression(self, g_type, **kwargs):
        """
        wrapper for the 3 basic standard progression
        """

        if g_type == "P":
            return self.P_progression(**kwargs)
        if g_type == "C":
            return self.C_progression(**kwargs)
        if g_type == "Ch":
            return self.C_progression(HANDLES = True, **kwargs)
        else:
            raise ValueError("g_type not supported in base_progression")

    def log_progression( self, g_type, bounds, points = 10, **kwargs):
        """
        get a standard progression evenly spread out on a log scale
        """

        select = np.geomspace(*bounds, num = points, dtype=int)
        select = set(select)
        select = [x for x in select]
        select.sort()
        select = np.array(select)

        return CollectionBuilder.base_progression(self, g_type, select = select, **kwargs)

    def log_chain_progression( self, gr_unit, bounds, points = 10, **kwargs):
        """
        get a standard progression evenly spread out on a log scale
        """

        select = np.geomspace(*bounds, num = points, dtype=int)
        select = set(select)
        select = [x for x in select]
        select.sort()
        select = np.array(select)

        return self.chain_progression(gr_unit, select = select, **kwargs)


#####################################à
#utils methods

def get_line_data(bounds = (3,10), target = "p", x_mode = "dist", **kwargs):
    """
    Simple wrapper for L graphs progression references   
    """

    an = Analyzer( mode = "first", opt_mode= "none")
    line_collection = CollectionBuilder().P_progression( bounds, analyzer=an)

    return line_collection.get_data( target = target, x_mode = x_mode)
