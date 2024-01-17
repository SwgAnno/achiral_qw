from achiralqw.analyze import performance
from achiralqw.graph import QWGraph, QWGraphBuilder
from achiralqw.graph import *
from achiralqw.plotter import *
from achiralqw.simulator import *
from achiralqw.trends import *
from achiralqw.article import *
from achiralqw.collection import *
import sys
import dill

import matplotlib.pyplot as plt

#plotting helper
#graphic comparison of optimum phase result
def check_optimum_phase( test_gr = QWGraphBuilder.Ring(6), mode = None, evt_mode = "first", qutip = True):

    if mode == "diag":
        plot_performance(test_gr, mode = mode, evt_mode = evt_mode, show = False)
    else :
        plot_performance(test_gr, evt_mode = evt_mode, show = False)

    if qutip:
        solver_mode = "qutip"
    else:
        solver_mode = "eigen"
    tp = TransportParameters(evt_mode = evt_mode, qutip = qutip, solver_mode=solver_mode)

    print("running optimization algorithm")
    if mode == "diag":
        res = optimum_phase_minimize(test_gr, diag = True, tp = tp)
    elif mode == "yolo":
        res = optimum_phase_yolo(test_gr, tp = tp)
    elif mode == "smart":
        res = optimum_phase_smart(test_gr, tp = tp)
    else:
        res = optimum_phase_minimize(test_gr, tp = tp)

    plt.scatter(res[0], performance(res[0]))

    plt.show()

def inspect_chain_first_maxima( gr_unit, bounds = (1,10), by = .1):

    sample = 4
    aTC = 4

    fig, axx = plt.subplots(1,4, figsize = (18,6))
    t_chain_progression_phases(gr_unit, bounds = bounds, sample_step = sample, ax = axx[0])
    axx[0].legend()

    chain1 = QWGraph.chain(gr_unit, int(bounds[1]/4))
    chain2 = QWGraph.chain(gr_unit, int(bounds[1]/2))
    chain3 = QWGraph.chain(gr_unit, bounds[1])

    plot_evo_vs_phase(chain1, TC = aTC, ax = axx[1], phase_by = by,plot_max= True)
    plot_evo_vs_phase(chain2, TC = aTC, ax = axx[2], phase_by = by, plot_max= True )
    plot_evo_vs_phase(chain3, TC = aTC, ax = axx[3], phase_by = by, plot_max= True)

    title_temp = "size = {}, d = {}"
    axx[1].set_title(title_temp.format(int(bounds[1]/4), chain1.distance()))
    axx[2].set_title(title_temp.format(int(bounds[1]/2), chain2.distance()))
    axx[3].set_title(title_temp.format(int(bounds[1]), chain3.distance()))

    prob_map = probability_colorbar_map()
    fig.subplots_adjust(.05, hspace = .25)
    cbar_ax = fig.add_axes([0.92, 0.1, 0.025, 0.8])
    fig.colorbar(prob_map, cbar_ax)

    #######################

    fig2, axx2 = plt.subplots(1, sample, figsize = (15,6))

    #phi vec stolen from  t_chain_progression_multi
    if gr_unit.code[0] == "C" and gr_unit.N % 2 == 0:
        phase_vec = np.arange(0, sample) * np.pi/sample
    else :
        phase_vec = np.arange(0, sample) * 2*np.pi/sample
    print(phase_vec)

    for phase, ax in zip(phase_vec, axx2):
        chain1.rephase( np.repeat( phase, chain1.get_phase_n() ))
        plot_evo_mat(chain1, TC = 1, filter = "target", ax = ax)
        ax.set_title( "phase = pi*%.2f" % (np.imag(np.log(phase))/np.pi))

    plt.show()

#todo: cleanup
def random():
    #a  = QWGraph.chain(QWGraphBuilder.Ring(3), 10)
    #a = QWGraphBuilder.Paralle l(3,2)
    #a = QWGraph.chain(a,2)
    ##a = a+a
    a = QWGraphBuilder.Ring(6)


    plot_performance(a)

    plot_evo_vs_derivative(a)

def general_test():
    a = QWGraphBuilder.Ring(7, HANDLES = False)
    b = QWGraphBuilder.Ring(6, HANDLES = True)

    plot_performance(a)
    plot_performance(b)

    #simple evolution
    plot_evo_mat(QWGraphBuilder.Ring(8), end = 10)


    #create a chain and compute krylov basis
    b = QWGraphBuilder.Ring(3)
    b = QWGraph.chain(b, 30)
    b.re_coord[0] = (2,1)
    plot_performance(b, mode = "diag")
    b.rephase( np.repeat(np.pi*-.5, len(b.re_coord)))
    print(b.mat)

    plot_krylov_basis(b)
    plot_krylov_couplings(b)
    plot_evo_mat_heatmap(b, TC = 4)

    plt.show()




#######################################


if __name__ == "__main__" :

    #plot_speedup_performance_multi(bounds = (4,24), target = "p", show = True)
    #plot_speedup_performance_multi_chain(QWGraphBuilder.Ring(3), bounds = (4,20), target = "p", show = True
    #plot_performance(QWGraph.chain(QWGraphBuilder.Ring(3), 10, speedup = 1), sample_step= 1000, mode = "diag")
    #plot_evo_mat(QWGraph.chain(QWGraphBuilder.Ring(3), 5, speedup = 10), end = 50)
    #plot_evo_mat(QWGraphBuilder.Line(6, speedup = 20), end = 20)

    # su_vec = [1,2,4,5]
    # a = QWGraphBuilder.Ring(3)
    # plot_performance(a, evt_mode = "TC", TC = 2)
    # chain_performance_multi_speedup(a, su_vec = su_vec,rep = 10, sample_step = 1000)
    #plot_evo_vs_phase( QWGraph.chain(a, 10, speedup = 5), by = .05, TC = 5, phase_by = .01)

    #a = QWGraph.chain(QWGraphBuilder.Ring(4), 1, speedup = 20)
    #a.plot()
    #b = QWGraphBuilder.Line(4, speedup = 20)
    #b.plot()

    #a = QWGraphBuilder.Ring(9)
    #a.plot()
    #plot_evo_vs_phase(a, TC = 2, phase_by = .01)
    #plot_performance(a,TC = 2, target = "t", sample_step = 500)

    a = QWGraphBuilder.Ring(7, HANDLES = False)
    b = QWGraphBuilder.Ring(6, HANDLES = True)
    #plot_performance_1_multi(a, TC_vec = [5,20], first = True, step = 1000)
    #plot_performance_odd_even(step = 1000)
    #comp_evo_vs_phase(a, [2,12], phase_by = .02, by = .05)


    #comp_size_progression([1,5,20], target = "t", bounds = (3, 30))
    #comp_performance_multi([b,a], TC_vec = [5,20], first = True, step = 1000)


    #a = QWGraphBuilder.Ring(8) + QWGraphBuilder.Ring(8)
    #a = QWGraphBuilder.SquareCut()
    #plot_performance(a, evt_mode = "TC", TC = 2)

    a = QWGraphBuilder.Ring(3)
    b = QWGraphBuilder.Ring(4)
    c = QWGraphBuilder.SquareCut()

    # c = QWGraph.chain(a,2)
    # print(c.mat)
    # plot_performance(c, mode = "diag")
    # print(c.re_coord)
    # print(c.start,c.target)
    #c.plot()
    plt.show()

    tp = TransportParameters( evt_mode = "first", diag = True)
    #prog = CollectionBuilder().C_progression( bounds = (4,38), step = 2, tp= tp)

    #print ( [gr.code for gr in prog.get_list()])

    #prog.transport_time_lm()

    #a = QWGraph.chain(a, 5)
    #plot_performance(a, mode = "time", evt_mode = "first", sample_step = 500)
    #plot_chain_progression(a, bounds = (1,10), target = "t", fix_phi = np.pi/2)
    #chain_progression(a, bounds = (1,20), target = "t", show = True, L_ref = True)
    #inspect_chain_first_maxima(b, bounds = (1,30), by = .02)
    #comp_t_chain_progression([a,b,c],mode = "best", target = "t", unit_bounds = (1,30), l_ref = True )
    #comp_best_t_chain_progression([a,b,c], target = "t", bounds = (4,20), l_ref = True )

    #t_chain_progression_phases(a, bounds = (1,40), sample_step = 4, mode = "first").legend()

    #multi_2_phases_example(sample = 200)

    #chain_ch_comp( bounds = (3,50))
    #line_speedup_perf_comp( bounds = (6,60), step = 5, su_bounds = (.1,2.5,1000), mode = "first", TC = 1).legend()

    #a.rephase(np.pi*-.5)
    #b = QWGraphBuilder.SquareCut()
    #plot_performance(b, evt_mode = "first")
    #plot_performance(b, evt_mode = "TC")
    #plot_evo_mat_heatmap(b)
    #plot_evo_vs_phase(b ,end = 10)

    #plot_evo_chain_speedup(a, 10, bounds = (0,100), step = .2,su_sample = 300, show = True)
    #plot_evo_line_speedup(10)
    #plot_evo_line_speedup(11)


    #parallel can chase the Bancho optimum
    # a = QWGraphBuilder.Parallel(4,2)
    # a = QWGraph.chain(a,32)
    # plot_performance(a, mode = "diag")

    #b = QWGraphBuilder.Ring(50)
    #b.re_coord[0] = (2,1)
    #plot_performance(b)
    #b.rephase(np.pi*.5)
    #print(b.mat)

    #b.krylov_basis(mode = "basis_plot")
    #b.krylov_basis(mode = "link_plot")
    #plot_evo_mat_heatmap(b)


    #a.krylov_basis(mode = "basis_plot")
    #a.krylov_basis(mode = "link_plot")

    tp = TransportParameters(evt_mode = "first")
    #plot_size_progression_multi( bounds = (4,40), step = 2, loglog = True, target = "p", tp = tp).legend()
    #plot_odd_even_progression( bounds = (3,40), target = "p", tp = tp).legend()
    #
    #plot_chain_progression_multi_loglog(bounds = (5,500), points = 50, target = "p", tp = tp, fast = True)
    #plot_chain_ch_progression(bounds = (5,30), loglog = False, target = "p", tp = tp)

    #plot_chain_progression_multi(bounds = (5,50), target = "p", tp = tp)

#    general_test()

    #odd_even_time_lm(HANDLES = False)
    #a = QWGraphBuilder.Ring(4)
    #a = QWGraphBuilder.SquareCut()
    #plot_performance(a, evt_mode = "TC")
    #plot_chain_progression(a,bounds = (1,20), target = "t")
    #time_chain_progression_lm(a)

    #t_size_progression_phases( bounds = (4,16), step = 2, tp = tp).legend()

    #pickling test

    def is_picklable(obj):
        try:
            dill.dumps(obj)

        except dill.PicklingError:
            return False
        return True
    
    def get_data(gr : QWGraph):
        return performance(gr, target = "t", tp = tp)


    print(is_picklable(get_data))

    if len(sys.argv) > 1:
        plt.savefig(sys.argv[1])
    else :
        file = input("Indicare il nome del file per salvare i risultati (lasciare vuoto per visualizzare)")   

        if file == "":
            plt.show() 
        else :
            plt.savefig(file)

    ######################
    #theta as char

    base = 950 -5

    greek =  [chr(base+i) for i in range(20)]
    print(greek)
    #####################
