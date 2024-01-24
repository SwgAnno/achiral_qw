from achiralqw.analyze import *
from achiralqw.graph import QWGraph, QWGraphBuilder
from achiralqw.article import plot_size_progression_multi

import pytest as pt
import numpy as np

from achiralqw.plotter import *


def test_label():

    qwgb = QWGraphBuilder()

    test_gr = qwgb.Ring(6, HANDLES = True)
    tp = TransportParameters(evt_mode= "TC", TC = 123456, opt_mode= "yolo")
    assert tp.get_label(test_gr) == "h(C6) {}=123456".format(chr(957))

    tp.TIME_CONSTANT = 5
    assert tp.get_label(test_gr) == "h(C6) {}=5".format(chr(957))
    
    test_gr = qwgb.Ring(7, HANDLES = False)
    tp = TransportParameters(evt_mode= "TC", TC = 9, opt_mode= "yolo")
    assert tp.get_label(test_gr) == "C7 {}=9".format(chr(957))

    assert tp.get_label(qwgb.Line(123)) == "P123 {}=9".format(chr(957))

def test_locate_max():
    qwgb = QWGraphBuilder()

    test_gr = qwgb.Ring(6, HANDLES = True, COMPUTE_EIGEN=True)

    tp = TransportParameters(evt_mode= "TC", TC = 20, opt_mode= "yolo")

    assert performance(test_gr, np.pi, tp = tp) == pt.approx(0)
    assert performance(test_gr, 0, tp = tp) == pt.approx(0.8907076163598934)

    tp.TIME_CONSTANT = 1
    assert performance(test_gr, np.pi, tp = tp) == pt.approx(0)
    assert performance(test_gr, 0, tp = tp) == pt.approx(0.7997209692969043)

    test_gr = qwgb.Ring(7, HANDLES = True,  COMPUTE_EIGEN=True)
    tp = TransportParameters(evt_mode= "TC", TC = 20, opt_mode= "yolo")

    assert performance(test_gr, np.pi/2, tp = tp) == pt.approx(0.7506993687215633)
    assert performance(test_gr, -np.pi/2, tp = tp) == pt.approx(0.9430288692063485)
    assert performance(test_gr, 0, tp = tp) == pt.approx(0.7075783036942821)

    tp.TIME_CONSTANT = 1
    assert performance(test_gr, np.pi/2, tp = tp) == pt.approx(0.750699368721558)
    assert performance(test_gr, -np.pi/2, tp = tp) == pt.approx(0.041007030503443904)
    assert performance(test_gr, 0, tp = tp) == pt.approx(0.3698988609910833)

def test_optimum_phase():
    
    gr = QWGraphBuilder.Ring(3)
    tp = TransportParameters(evt_mode="first")
    assert optimum_phase_minimize(gr, tp = tp)[0] == pt.approx(np.pi/2)
    assert optimum_phase_smart(gr, tp = tp)[0] == pt.approx(np.pi/2)
    assert optimum_phase_yolo(gr, tp = tp)[0] == pt.approx(np.pi/2, abs = 1e-3)

    gr = QWGraphBuilder.Ring(4)
    tp = TransportParameters(evt_mode="first")
    assert (optimum_phase_minimize(gr, tp = tp)[0] + np.pi/2) % (2*np.pi)== pt.approx(np.pi/2)
    assert (optimum_phase_smart(gr, tp = tp)[0] + np.pi/2) % (2*np.pi)== pt.approx(np.pi/2)
    assert (optimum_phase_yolo(gr, tp = tp)[0] + np.pi/2) % (2*np.pi)== pt.approx(np.pi/2)

    gr = QWGraphBuilder.Ring(4).chain(4)
    tp = TransportParameters(evt_mode = "TC", TC = 5)
    assert (optimum_phase_minimize(gr, tp = tp)[0] + np.pi/2) % (2*np.pi)== pt.approx(np.pi/2)
    assert (optimum_phase_smart(gr, tp = tp)[0] + np.pi/2) % (2*np.pi)== pt.approx(np.pi/2)
    assert (optimum_phase_yolo(gr, tp = tp)[0] + np.pi/2) % (2*np.pi)== pt.approx(np.pi/2)

def test_performance_reference():

    tp = TransportParameters(evt_mode = "first", opt_mode="smart")
    gr = QWGraphBuilder.SquareCut()
    assert performance_best(gr, tp = tp) == pt.approx(1)
    gr = QWGraphBuilder.Ring(3)
    assert performance_best(gr, tp = tp) == pt.approx(1)
    gr = QWGraphBuilder.Ring(4)
    assert performance_best(gr, tp = tp) == pt.approx(1)

    tp = TransportParameters(evt_mode = "first", opt_mode="smart", diag = True)
    gr = QWGraphBuilder.SquareCut().chain(10, HANDLES = True)
    assert performance_best(gr, tp = tp) == pt.approx(0.9188742)
    gr = QWGraphBuilder.Ring(3).chain(9, HANDLES = True)
    assert performance_best(gr, tp = tp) == pt.approx(0.96146705)
    gr = QWGraphBuilder.Ring(4).chain(8, HANDLES = True)
    assert performance_best(gr, tp = tp) == pt.approx(0.9392034)



# plitting helper
#confront locate_max results with actual evolution profile
def check_locate_max( test_gr = QWGraphBuilder.Ring(6), mode = None, qutip = True, TC = 1):

    fix, ax = plot_evo_vs_derivative(test_gr, show = False)

    tp = TransportParameters(solver_mode= "eigen" , TC = TC)
    points = np.zeros( (4,2) )

    tp.evt_mode = "first"
    print("############First###############")
    points[0,:] = locate_max(test_gr, tp = tp)
    tp.evt_mode = "TC"
    print("############TC###############")
    points[1,:] = locate_max(test_gr, tp = tp)
    
    if qutip:
        tp.solver_mode = "qutip"

        print("############qutip First###############")
        tp.evt_mode = "first"
        points[2,:] = locate_max(test_gr, tp = tp)
        print("############qutip TC###############")
        tp.evt_mode = "TC"
        points[3,:] = locate_max(test_gr, tp = tp)

    print(points)

    plt.scatter(points[:,0], points[:,1], color= "green")

    plt.show()

if __name__ == "__main__":
    
    qwgb = QWGraphBuilder()

    test1 = qwgb.Ring(5, COMPUTE_EIGEN= True)

    fig, axx = plt.subplots(3,1, figsize = (5,10))

    tp = TransportParameters(evt_mode="TC", TC = 10)
    plot_performance(test1,tp = tp,  ax = axx[0])
    tp = TransportParameters(evt_mode="TC", TC = 5)
    plot_performance(test1,evt_mode = "TC", TC = 5,  ax = axx[1])
    tp = TransportParameters(evt_mode="first")
    plot_performance(test1,evt_mode = "first",  ax = axx[2])

    plot_performance(test1|qwgb.Ring(8))

    #?????
    tp = TransportParameters(evt_mode= "first")
    plot_size_progression_multi( loglog = True, bounds = (4,100), tp = tp)

    #my_solver = QutipSESolver()
    #plot_evo_vs_derivative( test1, TC = 4, solver = my_solver, ax = axx[1])

    plt.show()