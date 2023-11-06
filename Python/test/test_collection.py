from achiralqw.collection import CollectionBuilder, QWGraphList, CachedQWGraphCollection
from achiralqw.simulator import Analyzer
from achiralqw.graph import QWGraphBuilder

import numpy as np

if __name__ == "__main__" :

    global MULTIPROCESS
    cb = CollectionBuilder()

    #test1 = cb.P_progression((4,50))
    #test1.evaluate()


    aList = QWGraphList()
    for i in range(4,10) :
        aList.append( QWGraphBuilder.Ring(i, HANDLES= True))

    aList.transport_time_lm()

    base_analyzer =  Analyzer( solver_mode= "eigen", mode = "first", diag = False)
    general_analyzer =  Analyzer( solver_mode= "eigen", mode = "first", diag = False)
    fast_analyzer =     Analyzer( solver_mode= "eigen", mode = "first", opt_mode = "fix", diag = True, gr = QWGraphBuilder.Ring(3))
    fast_analyzer.set_fix_phi( opt_mode = "smart")


    create_path = QWGraphBuilder.Line
    test_line = CachedQWGraphCollection( create_func = create_path, filename = "P_line", analyzer = base_analyzer)

    selection = np.arange(2, 100)
    test_line.update(selection)

    #fast evaluation with cached data!!!
    lm_test_range = [ np.arange(i*10, i*20) for i  in range(1, 6)]

    for selection in lm_test_range:
        test_line.transport_time_lm(selection)

    test_line.offload()

    