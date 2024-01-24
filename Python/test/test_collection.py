from achiralqw.analyze import TransportParameters
import achiralqw.collection as collection
from achiralqw.collection import CollectionBuilder, QWGraphCollection, QWGraphList, CachedQWGraphCollection
from achiralqw.graph import QWGraph, QWGraphBuilder

import numpy as np
import pytest 
import pathlib
import time

def ring_plus_tail( i : int) -> QWGraph:
    return QWGraphBuilder.Ring(i) + QWGraphBuilder.Line(8)

def C3_chain(i: int) -> QWGraph:
    unit = QWGraphBuilder.Ring(3)
    return unit.chain(i)

def test_create_collection():

    gr_list = [QWGraphBuilder.Ring(x) for x in range(3,30, 2)]

    coll = QWGraphList()
    for gr in gr_list:
        coll.append(gr)

    assert coll.get_name() == gr_list[0].code
    assert coll.get_list_x(gr_list, x_mode = "size") == [x for x in range(3,30,2)]
    assert coll.get_list_x(gr_list, x_mode = "dist") == [x for x in range(1,15,1)]

def test_collection_consistency():

    cb = CollectionBuilder()

    test_a = cb.P_progression( bounds=(5,15))
    test_b = cb.P_progression( select=[5,6,7,8,9,10,11,12,13,14,15])
    test_c = cb.P_progression_singleprocess(bounds= (5,15))

    if collection.MULTIPROCESS:
        test_d = cb.P_progression_multiprocess(bounds= (5,15))
    else:
        test_d = test_c
    for g1, g2, g3, g4 in zip(test_a.get_list(), test_b.get_list(),test_c.get_list(),test_d.get_list()):
        assert g1.code == g2.code
        assert g1.code == g3.code
        assert g1.code == g4.code

def test_collection_builder():
    cb = CollectionBuilder()

    #two base topologies
    test = cb.C_progression( bounds=(5,15))
    test = cb.C_progression( bounds=(5,15), HANDLES=True)
    
    #the three chains
    unit = QWGraphBuilder.Ring(3)
    test = cb.chain_progression(unit, select= [3,4,5,6,7,8,9]) 
    unit = QWGraphBuilder.Ring(4)
    test = cb.chain_progression(unit, bounds=(3,10)) 
    unit = QWGraphBuilder.SquareCut()
    test = cb.chain_progression(unit, select= [3,4,5,6,7,8,9]) 

    #custom creation
    select = [x for x in range(4,12,2)]
    l = cb.build_gr_list(create_func= ring_plus_tail, input_vec= select)
    test = cb.from_list(l)

def test_evaluate_collection():
    
    cb = CollectionBuilder()
    
    tp = TransportParameters(evt_mode = "first", opt_mode= "none")
    test = cb.P_progression(bounds=(5,10), tp = tp)
    test.evaluate()

    tp = TransportParameters(evt_mode = "first", opt_mode= "smart")
    test = cb.C_progression(bounds=(5,10), tp = tp)
    test.evaluate()

def test_linear_model():
    cb = CollectionBuilder()
    
    tp = TransportParameters(evt_mode = "first", opt_mode= "none")
    test = cb.P_progression(bounds=(5,10), tp = tp)
    params = test.transport_time_lm()

    print(params)
    assert params[0] == pytest.approx(0.5499988, rel = 1e-6)
    assert params[1] == pytest.approx(1.1884961, rel = 1e-6)

def test_prob_model():
    cb = CollectionBuilder()
    
    tp = TransportParameters(evt_mode = "first", opt_mode= "none")
    test = cb.C_progression(bounds=(5,10), step = 2, tp = tp)
    print( [gr.code for gr in test.get_list()])
    params = test.transport_prob_model(mode = "banchi2")

    print(params)
    assert params[0] == pytest.approx(1.8112682,  rel = 1e-6)
    assert params[1] == pytest.approx(-0.5107786, rel = 1e-6)

def test_cached_collection( tmp_path : pathlib.Path):

    tp = TransportParameters(evt_mode="first", opt_mode="smart")
    unit = QWGraphBuilder.Ring(3)
    tp.fix_phase(unit)

    collection_file = tmp_path / "test_collection"

    test = CachedQWGraphCollection(create_func=C3_chain, filename= collection_file, tp = tp)
    test.evaluate( [x for x in range(1,8)])
    test.offload()

    test = CachedQWGraphCollection(create_func=C3_chain, filename= collection_file, tp = tp)

    t0 = time.time()
    test.evaluate( [x for x in range(1,8)])
    t1 = time.time()
    test.update(   [x for x in range(1,8)])

    #did it run fast?
    assert t1 == pytest.approx(t0, abs=1e-1)

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

    