from abc import abstractmethod
from typing import List, Tuple
from achiralqw.analyze import TransportParameters, performance_best
from achiralqw.graph import QWGraph, QWGraphBuilder
import achiralqw.istarmap

import numpy as np
from numpy.typing import NDArray
import scipy.stats as stats
from scipy.optimize import curve_fit

import multiprocessing  as mp
import tqdm
import os, copy
import json

##############################
#global function that can be pickled and shared between processes

MULTIPROCESS = True

def get_perf(args):
    """
    This global function expects an array with three items in it:
    args[0] : QWGraph   (the relevant graph)
    args[1] : str       (the target of performance())
    args[2] : TransportParameters
    """
    return performance_best(args[0], target = args[1], tp = args[2])

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

##############################

class QWGraphCollection(object) :

    def __init__( self, tp : TransportParameters = TransportParameters(), name = None ) :

        self._tp = tp
        self._name = name


    @staticmethod
    def _get_data(gr_list, tp : TransportParameters = TransportParameters(), target = "p", x_mode = "dist", name = "QWgraphCollection") -> Tuple[List[int], NDArray]:
        """
        Router method for get_data() according to the global variable MULTIPROCESS
        """
        global MULTIPROCESS

        if MULTIPROCESS:
            return QWGraphCollection._get_data_multiprocess( gr_list = gr_list,      \
                                                            tp = tp,                \
                                                            target = target,        \
                                                            x_mode = x_mode,        \
                                                            name = name             )
        else :
            return QWGraphCollection._get_data_singleprocess( gr_list = gr_list,      \
                                                            tp = tp,                \
                                                            target = target,        \
                                                            x_mode = x_mode,        \
                                                            name = name             )

    @staticmethod
    def _get_data_singleprocess( gr_list, tp : TransportParameters = TransportParameters(), target = "p", x_mode = "dist", name = "QWgraphCollection") -> Tuple[List[int], NDArray] :

        N = len(gr_list)
        data = np.empty( N)

        prog_label = name +" " + target + " data:  {:2.1%}"

        for i in range(N):
            print (prog_label.format(i/N), end = "\r")

            data[i] = performance_best(gr_list[i], target=target, tp = tp)

        print(prog_label.format(1))

        x = QWGraphCollection.get_list_x(gr_list, x_mode = x_mode)

        return x , data

    @staticmethod
    def _get_data_multiprocess( gr_list, tp : TransportParameters = TransportParameters(), target = "p", x_mode = "dist", name = "QWgraphCollection") -> Tuple[List[int], NDArray]:

        N = len(gr_list)
        out = []

        def get_gr_data(gr : QWGraph, target = target, tp = tp):
            return performance_best(gr, target=target, tp = tp)

        n_proc = os.cpu_count()*2  

        greeting_string = name + " : Starting pool evaluation with {} process" 
        print(greeting_string.format(n_proc))

        with mp.Pool( n_proc) as pool:
            
            #todo: second and third argument are constant and read only, might as well share them
            zipped = zip(gr_list, [target]*N, [tp]*N)
            print("Evaluation")
            for _ in tqdm.tqdm(pool.imap(get_perf, zipped), total=N):
                out.append(_)
            data = np.array(out)

            pool.close()
            pool.join()

        x = QWGraphCollection.get_list_x(gr_list, x_mode = x_mode)

        return x , data

    @staticmethod
    def get_list_x(gr_list, x_mode = "dist") -> List[int]:
        out = []

        for gr in gr_list:
            if x_mode == "dist":
                out.append( gr.distance())
            elif x_mode == "size" :
                out.append( gr.N)
            else :
                raise ValueError("get_list_x mode not supported")
        
        return out

    @abstractmethod  
    def evaluate( self, select, target = "p", x_mode = "dist")  -> Tuple[List[int], NDArray]:
        """
        Evaluate the transport performance of the whole collection on a specific selection of elements
        """
        pass

    #todo make it return lm object and not just the two coefficients
    def transport_time_lm( self, select = None, x_mode = "dist") -> Tuple[float, float]:
        """
        Construct a linear model of the best transport time from the graph collection
        """
        
        x, data = self.evaluate(select, target = "t", x_mode= x_mode)

        #print(data)
        out = stats.linregress(x, data)

        print("m: ", out.slope, " +- ", out.stderr)
        print("q: ", out.intercept, " +- ", "????")
        print("r: ", out.rvalue)

        return out.slope, out.intercept

    def transport_prob_model( self, select = None, mode = "poly", x_mode = "dist") -> NDArray:
        """
        Very specific analysis tool: it extracts the slope of the probability progression on a loglog scale,
        which represent the power law exponent of its decay

        mode : str ["poly" , "banchi1", "banchi2", "banchi3", "banchi4", "banchi4log", "custom"]
        """

        modes = ["poly" , "banchi1", "banchi2", "banchi3", "banchi4", "banchi4log", "custom"]
        assert mode in modes, "Transport prob model mode not recognized"

        print("######### Transpor prob model: ", mode)

        x, data = self.evaluate(select, target = "p", x_mode = x_mode)

        #test over y = ax^2 + bx + c
        if mode == "poly" :
            param = np.polyfit(np.log(x), np.log(data), deg = 2)
            
            print("ax^2: ", param[0])
            print("bx: ", param[1])
            print("c: ", param[2])

            return param

        elif mode == "banchi1":
            # test over y = a(x)^-2/3 
            def model(x,a,b):
                return a*np.power(x, -2/3)

            #first with stegun bessel expansion
            p0 = [0.455469,-0.147317]
            param, cov_param = curve_fit(model, x, data, p0 = p0)

            print("a*x^-2/3: ", param[0], " +- ", cov_param[0,0])
            
            return param

        elif mode == "banchi2":
            # test over y = a(x)^-2/3 + b(x)^-4/3   
            def model(x,a,b):
                return a*np.power(x, -2/3) + b*np.power(x, -4/3)

            #first with stegun bessel expansion
            p0 = [0.455469,-0.147317]
            param, cov_param = curve_fit(model, x, data, p0 = p0)

            print("a*x^-2/3: ", param[0], " +- ", cov_param[0,0])
            print("b*x^-4/3: ", param[1], " +- ", cov_param[1,1])

            return param

        elif mode == "banchi3":
            # banchi model with 3 terms 
            def model(x,a,b,c):
                return a*np.power(x, -2/3) + b*np.power(x, -4/3) + c*np.power(x, -2)

            #first with stegun bessel expansion
            p0 = [0.455469,-0.147317, .25]
            param, cov_param = curve_fit(model, x, data, p0 = p0)

            print("a*x^-2/3: ", param[0], " +- ", cov_param[0,0])
            print("b*x^-4/3: ", param[1], " +- ", cov_param[1,1])
            print("c*x^-2: ", param[2], " +- ", cov_param[2,2])

            return param

        elif mode == "banchi4":
            # banchi model with 3 terms 
            def model(x,a,b,c,d):
                return a*np.power(x, -2/3) + b*np.power(x, -4/3) + c*np.power(x, -2.0) + d*np.power(x, -8/3)

            #first with stegun bessel expansion
            p0 = [0.455469,-0.147317, .25, -.5]
            param, cov_param = curve_fit(model, x, data, p0 = p0)

            print("a*x^-2/3: ", param[0], " +- ", cov_param[0,0])
            print("b*x^-4/3: ", param[1], " +- ", cov_param[1,1])
            print("c*x^-2: ", param[2], " +- ", cov_param[2,2])
            print("d*x^-8/3: ", param[3], " +- ", cov_param[3,3])

            return param

        elif mode == "banchi4log":
            # banchi model with 3 terms 
            def model(x,a,b,c,d):
                return np.log(a*np.power(x, -2/3) + b*np.power(x, -4/3) + c*np.power(x, -2.0) + d*np.power(x, -8/3))

            #first with stegun bessel expansion
            p0 = [0.455469,-0.147317, .25, -.5]
            param, cov_param = curve_fit(model, x, np.log(data), p0 = p0)

            print("a*x^-2/3: ", param[0], " +- ", cov_param[0,0])
            print("b*x^-4/3: ", param[1], " +- ", cov_param[1,1])
            print("c*x^-2: ", param[2], " +- ", cov_param[2,2])
            print("d*x^-8/3: ", param[3], " +- ", cov_param[3,3])

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

    def transport_prob_loglog_lm( self, select = None, x_mode = "dist") -> Tuple[float, float]:
        """
        Very specific analysis tool: it extracts the slope of the probability progression on a loglog scale,
        which represent the power law exponent of its decay
        """
        
        x, data = self.evaluate(select, target = "p", x_mode = x_mode)

        #print(data)
        out = stats.linregress(np.log(x), np.log(data))

        print("m: ", out.slope, " +- ", out.stderr)
        print("q: ", out.intercept, " +- ", "???")
        print("r: ", out.rvalue)

        return out.slope, out.intercept
        
    def set_transport_params( self, tp : TransportParameters ):
        self._tp = tp

    def get_transport_params( self) -> TransportParameters:
        return self._tp

    def get_name( self) -> str:
        return "noname QWgraphCollection" if self._name is None else self._name

class QWGraphList(QWGraphCollection):
    """
    QWGraphList supports QWgraphCollection primitives and wraps around a generic list of QWGraphs
    """    
    def __init__( self, **kwargs):
        super().__init__(**kwargs)
        self._gr_list : list[QWGraph] = []

    def evaluate( self, select = None, target = "p", x_mode = "dist") -> Tuple[List[int], NDArray]:
        """
        Evaluate the transport performance of the whole collection on a specific selection of elements
        """

        gr_list = self.get_list()

        if select is not None:
            gr_list = [gr_list[i] for i in select] 

        return QWGraphCollection._get_data(gr_list, tp = self.get_transport_params(), target = target, x_mode = x_mode, name = self.get_name())


    def __getitem__(self, key : int) -> QWGraph :
        return self.get_list()[key]
    
    def __setitem__(self, key : int, value : QWGraph) -> None:
        self.get_list()[key] = value

    def append(self, gr : QWGraph) -> None:
        self.get_list().append( gr )

    def get_list(self) -> List[QWGraph]:
        return self._gr_list
    
    def get_name(self) -> str:
        return self.get_list()[0].code if self._name == None else self._name

class CachedQWGraphCollection(QWGraphCollection):

    """
    Spinn-off of QWGraphCollection class with an internal json object 
    which you can load from or dump into a file.
    Stores previous computations
    """

    def __init__(self, create_func, filename : str, tp : TransportParameters = TransportParameters(), target = "p", x_mode = "dist", **kwargs):

        super().__init__(tp = tp, name = filename,  **kwargs)
        
        #todo: smart regex file name
        filename = filename # +.json

        if os.path.isfile("./{}.json".format(filename)):
            with open("{}.json".format(filename), "r") as file :
                self._data = json.load(file)

                #creating from an existing file: loading previous transport parameters
                tp = TransportParameters(   mode        = self._data["tp"]["mode"],
                                            opt_mode    = self._data["tp"]["opt_mode"],
                                            TC          = self._data["tp"]["TC"],
                                            diag        = self._data["tp"]["diag"])
                
                tp.fix_phi = self._data["tp"]["fix_phi"]

                self._tp = tp

        else :
            #initialize the new data dicionary, the various option infos are not going to be ever changed

            self._data = dict()
            self._data["name"] = filename

            self._data["tp"] = dict()
            self._data["tp"]["diag"] = tp.diag
            self._data["tp"]["TC"] = tp.TIME_CONSTANT
            self._data["tp"]["mode"] = tp.evt_mode
            self._data["tp"]["opt_mode"] = tp.opt_mode
            self._data["tp"]["fix_phi"] = tp.fix_phi

            self._data["options"] = dict()
            self._data["options"]["target"] = target
            self._data["options"]["x_mode"] = x_mode

            self._data["data"] = dict()

        self._create = create_func

    def offload(self, filename = None):

        if filename == None:
            filename = self._name

        filename += ".json"
        with open(filename, "w") as out :

            json.dump(self._data,out)

    def offload_data(self, select = None, filename = None):

        if filename == None:
            filename = self._name

        filename += ".data"
        with open(filename, "w") as out :

            if select is None:
                select = self._data["data"].keys()

            for k in select:
                p = self._data["data"].get(str(k), -1)
                out.write(str(k)+ " " + str(p) + "\n")

            out.close()

    def get_select_x(self, select, x_mode = "dist") -> List[int]:
        #estimate a linear behavour and pray

        x1 = 8
        x0 = 6
        gr1 = self._create(x1)
        gr0 = self._create(x0)


        if x_mode == "dist":
            y1 = gr1.distance()
            y0 = gr0.distance()
        elif x_mode == "size" :
            y1 = gr1.N
            y0 = gr0.N
        else :
            raise ValueError("get_select_x mode not supported")

        m = (y1 - y0)/(x1 - x0)
        q = y1 - m*x1
        out = [ m*x +q for x in select]

        return out


    def evaluate( self, select, target = "p", x_mode = "dist") -> Tuple[List[int], NDArray]:
        """
        Evaluate the transport performance of the whole collection on a specific selection of elements

        Retrive all the useful cached data and compute the missing ones
        target and x_mode appear as argument for interface uniformity but are going to be ignored
        """

        out_x = self.get_select_x(select = select, x_mode = x_mode)
        out_data = [-1]* len(select)

        missing_ids = []
        missing = []
        for i, id in enumerate(select) :
            perf = self._data["data"].get(str(int(id)))

            if perf is None  :
                missing_ids.append(id)
                missing.append(i)
            else:
                out_data[i] = perf

        #print(missing)

        if len(missing) > 0 :

            gr_list = CollectionBuilder.build_gr_list( self._create, missing_ids)
            missing_x, missing_data = QWGraphCollection._get_data(  gr_list, tp = self.get_transport_params(),   \
                                                                    target = self._data["options"]["target"],    \
                                                                    x_mode = self._data["options"]["x_mode"],    \
                                                                    name = self.get_name()                       )
       
            #update cache and out_data with new new values
            for i, perf in zip(missing, missing_data):
                self._data["data"][ str(int(select[i]))] = perf
                out_data[i] = perf

        return out_x, out_data
    
    def update(self, select) -> Tuple[List[int], NDArray]:
        """
        Update a selection of the data cache entries
        """

        gr_list = []
        for id in select:
            gr_list.append( self._create(id))
       
        x, data = QWGraphCollection._get_data(  gr_list, tp = self.get_transport_params(),   \
                                                target = self._data["options"]["target"],    \
                                                x_mode = self._data["options"]["x_mode"],    \
                                                name = self.get_name()                       )
        
        for id, perf in zip(x,data):
            self._data["data"][str(int(id))] = perf

        return x, data

    def get_name(self) -> str:
        return self._data["name"]




#todo:  why not make everything static
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

    @staticmethod
    def build_gr_list( create_func, input_vec) -> List[QWGraph]:
        """
            Helper which returns a list of QWGraph build with the function passed as input,
            possibly employing multiple processes
        """

        global MULTIPROCESS
        gr_list = []

        if MULTIPROCESS:
            n_proc = os.cpu_count()*2  
            greeting_string = "P progression: Starting pool creation with {} process" 
            print(greeting_string.format(n_proc))

            with mp.Pool( n_proc) as pool:
                for _ in tqdm.tqdm(pool.imap(create_func, input_vec ), total=len(input_vec)):
                    gr_list.append(_)

                pool.close()
                pool.join()
        else :
            for l in input_vec:
                gr_list.append( create_func(l))

        return gr_list


    def from_list(self, gr_list , tp : TransportParameters = TransportParameters()) -> QWGraphList :

        collection = QWGraphList( tp = tp)

        for gr in gr_list:
            collection.append( gr )

        return collection

    def P_progression_singleprocess(self, bounds = None, step = 1, select = None, tp : TransportParameters = TransportParameters()) -> QWGraphList:

        collection = QWGraphList( tp = tp)
        collection.get_transport_params().opt_mode = "none"

        assert bounds or np.any(select)

        if np.any(select):
            drange = select
        else :
            drange = np.arange(bounds[0], bounds[1], step)

        for d in drange :
            collection.append( QWGraphBuilder.Line(d))

        return collection

    def P_progression_multiprocess(self, bounds = None, step = 1, select = None, tp : TransportParameters = TransportParameters(), **kwargs) -> QWGraphList :

        collection = QWGraphList( tp=tp)
        collection.get_transport_params().opt_mode = "none"

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
                collection.append(_)

            pool.close()
            pool.join()

        return collection



    def C_progression_singleprocess(self, bounds = None, step = 1, odd = False, select = None, tp : TransportParameters = TransportParameters(), HANDLES = False, **kwargs) -> QWGraphList :

        collection = QWGraphList( tp=tp)

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
            collection.append( QWGraphBuilder.Ring(d, HANDLES = HANDLES, **kwargs))

        return collection

    def C_progression_multiprocess(self, bounds = None, step = 1, odd = False, select = None, tp : TransportParameters = TransportParameters(), HANDLES = False, 
     **kwargs) -> QWGraphList:

        collection = QWGraphList( tp=tp)

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
                    collection.append(_)
            else :
                for _ in tqdm.tqdm(pool.imap(create_c, drange ), total=len(drange)):
                    collection.append(_)

            pool.close()
            pool.join()

        return collection

    def chain_progression_singleprocess(self, gr_unit : QWGraph, bounds = None, step = 1, select = None, tp : TransportParameters = TransportParameters(), **kwargs) -> QWGraphList:

        collection = QWGraphList( tp = tp)

        assert bounds or np.any(select)

        if np.any(select):
            drange = unit_list(select, gr_unit)
        else :
            drange = unit_list_bounds( bounds, gr_unit)

        for d in drange :
            collection.append( gr_unit.chain( rep = d, **kwargs))

        return collection

    def chain_progression_multiprocess(self, gr_unit : QWGraph, bounds = None, step = 1, select = None, tp : TransportParameters = TransportParameters(), **kwargs) -> QWGraphList:
        collection = QWGraphList( tp = tp)

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
                collection.append(_)

            pool.close()
            pool.join()

        return collection

    def base_progression(self, g_type, **kwargs) -> QWGraphList:
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

    @staticmethod
    def log_selection(bounds, points = 10) -> List[int]:
        """
        return a integer geomspace selection with different indexes
        """
        select = np.geomspace(*bounds, num = points, dtype=int)
        select = set(select)
        select = [x for x in select]
        select.sort()
        select = np.array(select)

        return select

    def log_progression( self, g_type, bounds, points = 10, **kwargs) -> QWGraphList:
        """
        get a standard progression evenly spread out on a log scale
        """

        select = CollectionBuilder.log_selection(bounds = bounds, points = points)
        return CollectionBuilder.base_progression(self, g_type, select = select, **kwargs)

    def log_chain_progression( self, gr_unit, bounds, points = 10, **kwargs) -> QWGraphList:
        """
        get a standard progression evenly spread out on a log scale
        """

        select = np.geomspace(*bounds, num = points, dtype=int)
        select = set(select)
        select = [x/gr_unit.distance() for x in select]
        select.sort()
        select = np.array(select)

        return self.chain_progression(gr_unit, select = select, **kwargs)


#####################################à
#utils methods


#todo remove get_line data in favour of CachedQWgraphCollection


def get_line_data(bounds = (3,10), target = "p", x_mode = "dist", **kwargs) -> Tuple[List[int], NDArray]:
    """
    Simple wrapper for L graphs progression references   
    """

    tp = TransportParameters( evt_mode= "first", opt_mode="none")
    line_collection = CollectionBuilder().P_progression( bounds, tp = tp)

    return line_collection.evaluate( target = target, x_mode = x_mode)
