from achiralqw.collection import CollectionBuilder, QWGraphList, CachedQWGraphCollection
from achiralqw.simulator import Analyzer
from achiralqw.graph import QWGraphBuilder

import numpy as np

create_path = QWGraphBuilder.Line

def C3_chain(id):
    return QWGraphBuilder.Ring(3).chain(id)
def C4_chain(id):
    return QWGraphBuilder.Ring(4).chain(id)
def SC_chain(id):
    return  QWGraphBuilder.SquareCut().chain(id)

def C_odd(id):
    return QWGraphBuilder.Ring(2*id + 1)
def C_even(id):
    return QWGraphBuilder.Ring(2*id )

def load_P_cache():
    return  CachedQWGraphCollection( create_func = create_path,   filename = "P_first_fast")

def load_C3_cache():
    return  CachedQWGraphCollection( create_func = C3_chain,   filename = "C3_chain_first_fast")

def load_C4_cache():
    return  CachedQWGraphCollection( create_func = C4_chain,   filename = "C4_chain_first_fast")

def load_DiC_cache():
    return  CachedQWGraphCollection( create_func = SC_chain,   filename = "DiC4_chain_first_fast")

def load_C_odd_cache():
    return CachedQWGraphCollection( create_func = C_odd,   filename = "C_odd_first_fast")

def load_C_even_cache():
    return CachedQWGraphCollection( create_func = C_even,   filename = "C_even_first_fast")

def build_chain_cache():
    base_analyzer =  Analyzer( solver_mode= "eigen", mode = "first", diag = False)

    fast_analyzer_C3 =     Analyzer( solver_mode= "eigen", mode = "first", opt_mode = "fix", diag = True, gr = QWGraphBuilder.Ring(3))
    fast_analyzer_C3.set_fix_phi( opt_mode = "smart")
    fast_analyzer_C4 =     Analyzer( solver_mode= "eigen", mode = "first", opt_mode = "fix", diag = True, gr = QWGraphBuilder.Ring(4))
    fast_analyzer_C4.set_fix_phi( opt_mode = "smart")
    fast_analyzer_SquareCut =     Analyzer( solver_mode= "eigen", mode = "first", opt_mode = "fix", diag = True, gr = QWGraphBuilder.SquareCut())
    fast_analyzer_SquareCut.set_fix_phi( opt_mode = "smart")


    create_path = QWGraphBuilder.Line

    cached_line = CachedQWGraphCollection( create_func = create_path,   filename = "P_first_fast",          analyzer = base_analyzer)
    cached_C3   = CachedQWGraphCollection( create_func = C3_chain,      filename = "C3_chain_first_fast",   analyzer = fast_analyzer_C3)
    cached_C4   = CachedQWGraphCollection( create_func = C4_chain,      filename = "C4_chain_first_fast",   analyzer = fast_analyzer_C4)
    cached_DiC   = CachedQWGraphCollection( create_func = SC_chain,     filename = "DiC4_chain_first_fast", analyzer = fast_analyzer_SquareCut)

    selection = np.arange(2, 50)

    cached_line.evaluate(selection)
    cached_C3.evaluate(selection)
    cached_C4.evaluate(selection)
    cached_DiC.evaluate(selection)

    #fast evaluation with cached data!!!
    lm_test_range = [ np.arange(i*10, i*20) for i  in range(1, 3)]

    for selection in lm_test_range:
        cached_line.transport_time_lm(selection)
        cached_C3.transport_time_lm(selection)
        cached_C4.transport_time_lm(selection)
        cached_DiC.transport_time_lm(selection)

    cached_line.offload()
    cached_C3.offload()
    cached_C4.offload()
    cached_DiC.offload()

def build_C_cache():

    fast_analyzer_C3 =     Analyzer( solver_mode= "eigen", mode = "first", opt_mode = "fix", diag = True, gr = QWGraphBuilder.Ring(3))
    fast_analyzer_C3.set_fix_phi( opt_mode = "smart")
    fast_analyzer_C4 =     Analyzer( solver_mode= "eigen", mode = "first", opt_mode = "fix", diag = True, gr = QWGraphBuilder.Ring(4))
    fast_analyzer_C4.set_fix_phi( opt_mode = "smart")

    cached_odd = CachedQWGraphCollection( create_func = C_odd,   filename = "C_odd_first_fast",  analyzer = fast_analyzer_C3)
    cached_even = CachedQWGraphCollection( create_func = C_even,   filename = "C_even_first_fast",  analyzer = fast_analyzer_C4)

    selection = np.arange(2, 50)

    cached_odd.evaluate(selection)
    cached_even.evaluate(selection)

    cached_odd.offload()
    cached_even.offload()


if __name__ == "__main__" :

    #build_chain_cache()
    #build_C_cache()

    
    select = np.arange(1,100)
    fast_analyzer_C4 =     Analyzer( solver_mode= "eigen", mode = "first", opt_mode = "fix", diag = True, gr = QWGraphBuilder.Ring(4))
    fast_analyzer_C4.set_fix_phi(0)
    cached_C4   = CachedQWGraphCollection( create_func = C4_chain,      filename = "C4_chain_first_fast",   analyzer = fast_analyzer_C4)
    x, data = cached_C4.evaluate(select = select)

    cached_C4.offload()