from achiralqw.graph import QWGraph, QWGraphBuilder
import numpy as np


def test_builder():
    gb = QWGraphBuilder()

    a = gb.Line(3, COMPUTE_EIGEN= True)
    a = gb.Line(30, COMPUTE_EIGEN= True)
    a = gb.Line(300, COMPUTE_EIGEN= True)

    b = gb.Ring(3, COMPUTE_EIGEN= True)
    b = gb.Ring(4, COMPUTE_EIGEN= True)
    b = gb.Ring(30, COMPUTE_EIGEN= True)
    b = gb.Ring(40, COMPUTE_EIGEN= True)
    b = gb.Ring(300, COMPUTE_EIGEN= True)
    b = gb.Ring(400, COMPUTE_EIGEN= True)

    c = gb.SquareCut(COMPUTE_EIGEN= True)

    d = gb.Parallel(paths= 5, p_len= 2, COMPUTE_EIGEN= True)
    d = gb.Parallel(paths= 10, p_len= 3, COMPUTE_EIGEN= True)
    d = gb.Parallel(paths= 100, p_len=  4, COMPUTE_EIGEN= True)

def test_line():

    gb = QWGraphBuilder()

    l = gb.Line(5)

    assert l.distance() == 4
    assert l.code == "P5"
    assert l.N == 5
    assert l.re_coord == []

    mat = [[0,-1,0,0,0 ],
           [-1,0,-1,0,0],
           [0,-1,0,-1,0],
           [0,0,-1,0,-1],
           [0,0,0,-1,0 ]]
    mat = np.array(mat)

    assert np.sum(np.abs(l.mat - mat)) < 1e-10

def test_ring():

    gb = QWGraphBuilder()
    t = gb.Ring(3, E = 2)

    assert t.distance() == 1
    assert t.code == "C3"
    assert t.N == 3
    assert t.re_coord == [(1,2)]
    assert (t.start, t.target) == (0,2)

    mat = [[2,-1,-1],
           [-1,2,-1],
           [-1,-1,2]]
    mat = np.array(mat)
    assert np.sum(np.abs(t.mat - mat)) < 1e-10

    t = gb.Ring(4, E = 1)

    assert t.distance() == 2
    assert t.code == "C4"
    assert t.N == 4
    assert t.re_coord == [(2,3)]
    assert (t.start, t.target) == (0,3)

    #Third and fourth siets are swapped in order to have the target site as last
    mat = [[1,-1,-1,0],
           [-1,1,0,-1],
           [-1,0,1,-1],
           [0,-1,-1,1]]
    mat = np.array(mat)
    assert np.sum(np.abs(t.mat - mat)) < 1e-10

def test_squarecut():
    gb = QWGraphBuilder()
    t = gb.SquareCut(E = 3)

    assert t.distance() == 1
    assert t.code == "DiC4"
    assert t.N == 4
    assert t.re_coord == [(1,3),(2,3)]
    assert (t.start, t.target) == (0,3)

    mat = [[3,-1,-1,-1],
           [-1,3,0,-1],
           [-1,0,3,-1],
           [-1,-1,-1,3]]
    mat = np.array(mat)
    assert np.sum(np.abs(t.mat - mat)) < 1e-10

def test_graph_manipulation():
    gr = QWGraphBuilder

    a = gr.Ring(3)
    b = gr.Line(1)
    c = gr.Ring(4)

    t = a|b
    assert t.distance() == 1
    assert t.code == "C3|P1"
    assert t.N == 3
    assert t.re_coord == [(1,2)]
    assert (t.start, t.target) == (0,2)

    mat = [[0,-1,-1],
           [-1,0,-1],
           [-1,-1,0]]
    mat = np.array(mat)
    assert np.sum(np.abs(t.mat - mat)) < 1e-10

    t = a+c
    assert t.distance() == 4
    assert t.code == "C3+C4"
    assert t.N == 7
    assert t.re_coord == [(1,2), (5,6)]

    mat = [[0,-1,-1,0,0,0,0],
           [-1,0,-1,0,0,0,0],
           [-1,-1,0,-1,0,0,0],
           [0,0,-1,0,-1,-1,0],
           [0,0,0,-1,0,0,-1],
           [0,0,0,0-1,0,0,-1],
           [0,0,0,0,-1,-1,0],
           ]
    mat = np.array(mat)
    assert np.sum(np.abs(t.mat - mat)) < 1e-10

    t = a|c
    assert t.distance() == 3
    assert t.code == "C3|C4"
    assert t.N == 6
    assert t.re_coord == [(1,2), (4,5)]

    mat = [[0,-1,-1,0,0,0],
           [-1,0,-1,0,0,0],
           [-1,-1,0,-1,-1,0],
           [0,0,-1,0,0,-1],
           [0,0,0-1,0,0,-1],
           [0,0,0,-1,-1,0],
           ]
    mat = np.array(mat)
    assert np.sum(np.abs(t.mat - mat)) < 1e-10

    t = a.add_handles(1)
    assert t.distance() == 3
    assert t.code == "h(C3)"
    assert t.N == 5
    assert t.re_coord == [(2,3)]

    mat = [[0,-1, 0, 0,0],
           [-1,0,-1,-1,0],
           [0,-1,0,-1, 0],
           [0,-1,-1,0,-1],
           [0, 0,0,-1, 0]]
    mat = np.array(mat)
    assert np.sum(np.abs(t.mat - mat)) < 1e-10

    t = a.chain(10, HANDLES = True)
    assert t.distance() == 12
    assert t.N == 23
    assert t.re_coord == [ (2+2*i, 3+2*i) for i in range(10)]

def test_rephase():

    gb = QWGraphBuilder()

    t = gb.Ring(3)
    t.rephase(  np.pi/2)
    t.rephase( [np.pi/2])
    t.rephase( np.array([np.pi/2]))

    mat = [[0,-1,-1],
           [-1,0,-1j],
           [-1,1j,0]]
    mat = np.array(mat)
    assert np.sum(np.abs(t.mat - mat)) < 1e-10

    t = gb.SquareCut()
    t.rephase( [np.pi/2, -np.pi/2])
    t.rephase( np.array([np.pi/2, -np.pi/2]))

    mat = [[ 0,-1,-1,-1],
           [-1,0,0,-1j],
           [-1,0,0,1j],
           [-1,1j,-1j,0]]
    mat = np.array(mat)
    assert np.sum(np.abs(t.mat - mat)) < 1e-10