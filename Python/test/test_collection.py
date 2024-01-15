from achiralqw.analyze import TransportParameters
from achiralqw.collection import CollectionBuilder, QWGraphList, CachedQWGraphCollection
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

    base_params =  TransportParameters( solver_mode= "eigen", evt_mode = "first", diag = False)
    general_params =  TransportParameters( solver_mode= "eigen", evt_mode = "first", diag = False)
    fast_params =     TransportParameters( solver_mode= "eigen", evt_mode = "first", opt_mode = "smart", diag = True)
    fast_params.fix_phase( gr = QWGraphBuilder.Ring(3))


    create_path = QWGraphBuilder.Line
    test_line = CachedQWGraphCollection( create_func = create_path, filename = "P_line", tp = base_params)

    selection = np.arange(2, 100)
    test_line.update(selection)

    #fast evaluation with cached data!!!
    lm_test_range = [ np.arange(i*10, i*20) for i  in range(1, 6)]

    for selection in lm_test_range:
        test_line.transport_time_lm(selection)

    test_line.offload()

    