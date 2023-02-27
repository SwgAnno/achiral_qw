from achiralqw.graph import QWGraph, QWGraphBuilder
from achiralqw.article import plot_base_progression

from achiralqw.plotter import *


if __name__ == "__main__":
    
    qwgb = QWGraphBuilder()

    test1 = qwgb.Ring(5, COMPUTE_EIGEN= True)

    fig, axx = plt.subplots(3,1, figsize = (5,10))

    plot_performance(test1,an_mode = "TC", TC = 10,  ax = axx[0])
    plot_performance(test1,an_mode = "TC", TC = 5,  ax = axx[1])
    plot_performance(test1,an_mode = "first",  ax = axx[2])

    #?????
    plot_base_progression( g_type= "P", bounds = (2,100))

    #my_solver = QutipSESolver()
    #plot_evo_vs_derivative( test1, TC = 4, solver = my_solver, ax = axx[1])

    plt.show()