import matplotlib.pyplot as plt


from achiralqw.graph import QWGraph, QWGraphBuilder
from achiralqw.simulator import EigenSESolver,QutipSESolver
from achiralqw.plotter import *




if __name__ == "__main__":
    
    qwgb = QWGraphBuilder()

    test1 = qwgb.Ring(5, COMPUTE_EIGEN= True)

    fig, axx = plt.subplots(3,1, figsize = (5,10))

    my_solver = EigenSESolver()
    plot_evo_vs_derivative( test1, TC = 4, solver = my_solver, ax = axx[0])
    plot_evo_mat(test1, TC = 4, solver = my_solver, ax = axx[1])
    plot_evo_mat_heatmap(test1, TC = 4, solver = my_solver,fig = fig, ax = axx[2])

    #my_solver = QutipSESolver()
    #plot_evo_vs_derivative( test1, TC = 4, solver = my_solver, ax = axx[1])

    plt.show()