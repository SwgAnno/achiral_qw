from Graph import QWGraph as qwg
from Graph import *
from plotter import *
from simulator import *
from trends import *
from article import *
from collection import *
import sys

import matplotlib.pyplot as plt

#global parameters for plotter methods
TC = 5
qut = False

#confront locate_max results with actual evolution profile
def check_locate_max( test_gr = qwg.Ring(6), mode = None, qutip = True, TC = 1):

    fix, ax = plot_evo_vs_derivative(test_gr, show = False)

    tester = Analyzer(test_gr, qutip = False, TC = TC)
    points = np.zeros( (4,2) )

    tester.mode = "first"
    print("############First###############")
    points[0,:] = tester.locate_max()
    tester.mode = "TC"
    print("############TC###############")
    points[1,:] = tester.locate_max()
    
    if qutip:
        tester.set_qutip(True)

        print("############qutip First###############")
        tester.mode = "first"
        points[2,:] = tester.locate_max()
        print("############qutip TC###############")
        tester.mode = "TC"
        points[3,:] = tester.locate_max()

    print(points)

    plt.scatter(points[:,0], points[:,1], color= "green")

    plt.show()


#comparison between eigenvalue method and Qutip solver on computed evolution
def check_evo_vs_qutip( test_gr = qwg.Ring(6), l = 0, start = 0, end = None, by = .1, deriv = True):
    if not end:
        global TC
        end = test_gr.N * TC

    seq = np.arange(start,end,by)
    if not deriv: 
        solver = SESolver(test_gr, qutip = False)
        evo = solver.target_p(seq)

        solver.set_qutip(True)
        evo_q = solver.target_p(seq)

    else :
        solver = SESolver(test_gr, qutip = False)
        evo = solver.target_p_prime(seq)

        solver.set_qutip(True)
        evo_q = solver.target_p_prime(seq)

    print( np.abs(evo-evo_q) < .0001)

    fig, ax = plt.subplots()

    ax.plot(seq, evo)
    ax.plot(seq, evo_q)
    
    ax.set_xlabel('Time')
    ax.set_ylabel('P')

    ax.legend(["P", "P qutip"])
    plt.show()

#see how qutip (or the old method) works with random poits evaluation
def check_evo_vs_qutip_scatter( test_gr = qwg.Ring(6), l = 0, start = 0, end = None, by = .1, deriv = False):
    if not end:
        global TC
        end = test_gr.N * TC

    seq = np.arange(start,end,by)
    seq_q = np.linspace(start,end, 20)
    test = np.random.rand(40)
    test = test * (end-start) + start

    test_qut = test[0:20]
    test_old = test[20:40]
    
    if not deriv: 
        solver = SESolver(test_gr, qutip = False)
        evo = solver.target_p(seq)
        evo_scat = solver.target_p(test_old)

        solver.set_qutip(True)
        evo_q = solver.target_p(test_qut)
        evo_seq = solver.target_p(seq_q)

    else :
        solver = SESolver(test_gr, qutip = False)
        evo = solver.target_p_prime(seq)
        evo_scat = solver.target_p_prime(test_old)

        solver.set_qutip(True)
        evo_q = solver.target_p_prime(test_old)
        evo_seq = solver.target_p_prime(seq_q)

    fig, ax = plt.subplots()

    ax.plot(seq, evo)
    ax.scatter(test_old, evo_scat, color = "green")
    ax.scatter(test_qut, evo_q, color = "red")
    ax.scatter(seq_q, evo_seq, color = "yellow")
    
    ax.set_xlabel('Time')
    ax.set_ylabel('P')

    ax.legend(["P", "Scatter evaluation noQ","Scatter evaluation Q","Scat ordered evaluation Q"])
    plt.show()

#single qutip evaluation check
def check_qutip( test_gr = qwg.Ring(6), t_0 = [2], deriv = False):

    solver = SESolver(test_gr, qutip = True)

    for i in range(2):
        if not deriv:
            print( solver.target_p(t_0))
        else :
            print( solver.target_p_prime(t_0))

#graphic comparison of optimum phase result
def check_optimum_phase( test_gr = qwg.Ring(6), mode = None, an_mode = "first", qutip = True):

    if mode == "diag":
        plot_performance(test_gr, mode = mode, an_mode = an_mode, show = False)
    else :
        plot_performance(test_gr, an_mode = an_mode, show = False)

    tester = Analyzer(test_gr, mode = an_mode, qutip = qutip)

    print("running optimization algorithm")
    if mode == "diag":
        res = tester.optimum_phase_minimize(diag = True)
    elif mode == "yolo":
        res = tester.optimum_phase_yolo()
    elif mode == "smart":
        res = tester.optimum_phase_smart()
    else:
        res = tester.optimum_phase_minimize()

    plt.scatter(res[0], tester.performance(res[0]))

    plt.show()

def check_line_bessel_dist(l = 5, end = 10):

    grid = np.arange(0,end, 1)
    eval = np.zeros( len(grid))

    for t in range(len(grid)):
        for i in range(l):
            eval[t] = eval[t] + line_evo_bessel(l-1,i, grid[t])

    print(eval)

def check_ring_bessel_dist(l = 5, end = 10):

    grid = np.arange(0,end, 1)
    eval = np.zeros( len(grid))

    for t in range(len(grid)):
        for i in range(l):
            eval[t] = eval[t] + ring_evo_bessel(l,i, grid[t])

    print(eval)

def inspect_chain_first_maxima( gr_unit, bounds = (1,10), by = .1):

    sample = 4
    aTC = 4

    fig, axx = plt.subplots(1,4, figsize = (18,6))
    t_chain_progression_phases(gr_unit, bounds = bounds, sample_step = sample, ax = axx[0])
    axx[0].legend()

    chain1 = qwg.chain(gr_unit, int(bounds[1]/4))
    chain2 = qwg.chain(gr_unit, int(bounds[1]/2))
    chain3 = qwg.chain(gr_unit, bounds[1])

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
        phi_seq = np.arange(0, sample) * np.pi/sample
    else :
        phi_seq = np.arange(0, sample) * 2*np.pi/sample
    print(phi_seq)

    phase_vec = np.exp( 1j * phi_seq)

    for phase, ax in zip(phase_vec, axx2):
        chain1.rephase( np.repeat( phase, chain1.get_phase_n() ))
        plot_evo_mat(chain1, TC = 1, filter = "target", ax = ax)
        ax.set_title( "phase = pi*%.2f" % (np.imag(np.log(phase))/np.pi))

    plt.show()

#todo: cleanup
def random():
    #a  = qwg.chain(qwg.Ring(3), 10)
    #a = qwg.Paralle l(3,2)
    #a = qwg.chain(a,2)
    ##a = a+a
    a = qwg.Ring(6)


    plot_performance(a)

    plot_evo_vs_derivative(a)



#######################################


if __name__ == "__main__" :
    #plot_line_vs_bessel(l = 5, trace_conn = False)
    #bessel_progression(bounds = (2,50), target = "p", L_ref = True)

    #plot_line_bessel_evolution(5,end = 10)
    #plot_evo_mat(qwg.Ring(8), end = 10)
    #check_line_bessel_dist(5)

    #plot_speedup_performance_multi(bounds = (4,24), target = "p", show = True)
    #plot_speedup_performance_multi_chain(qwg.Ring(3), bounds = (4,20), target = "p", show = True
    #plot_performance(qwg.chain(qwg.Ring(3), 10, speedup = 1), sample_step= 1000, mode = "diag")
    #plot_evo_mat(qwg.chain(qwg.Ring(3), 5, speedup = 10), end = 50)
    #plot_evo_mat(qwg.Line(6, speedup = 20), end = 20)

    # su_vec = [1,2,4,5]
    # a = qwg.Ring(3)
    # plot_performance(a, an_mode = "TC", TC = 2)
    # chain_performance_multi_speedup(a, su_vec = su_vec,rep = 10, sample_step = 1000)
    #plot_evo_vs_phase( qwg.chain(a, 10, speedup = 5), by = .05, TC = 5, phase_by = .01)

    #a = qwg.chain(qwg.Ring(4), 1, speedup = 20)
    #a.plot()
    #b = qwg.Line(4, speedup = 20)
    #b.plot()

    #a = qwg.Ring(9)
    #a.plot()
    #plot_evo_vs_phase(a, TC = 2, phase_by = .01)
    #plot_performance(a,TC = 2, target = "t", sample_step = 500)

    a = qwg.Ring(7, HANDLES = False)
    b = qwg.Ring(6, HANDLES = True)
    #plot_performance_1_multi(a, TC_vec = [5,20], first = True, step = 1000)
    #plot_performance_odd_even(step = 1000)
    #comp_evo_vs_phase(a, [2,12], phase_by = .02, by = .05)


    #comp_size_progression([1,5,20], target = "t", bounds = (3, 30))
    #comp_performance_multi([b,a], TC_vec = [5,20], first = True, step = 1000)


    #a = qwg.Ring(8) + qwg.Ring(8)
    #a = qwg.SquareCut()
    #plot_performance(a, an_mode = "TC", TC = 2)

    a = qwg.Ring(3)
    b = qwg.Ring(4)
    c = qwg.SquareCut()

    # c = qwg.chain(a,2)
    # print(c.mat)
    # plot_performance(c, mode = "diag")
    # print(c.re_coord)
    # print(c.start,c.target)
    #c.plot()
    plt.show()

    an = Analyzer( mode = "first", diag = True)
    prog = CollectionBuilder().C_progression( bounds = (4,38), step = 2, analyzer= an)

    print ( [gr.code for gr in prog.get_list()])

    #prog.transport_time_lm()

    #a = qwg.chain(a, 5)
    #plot_performance(a, mode = "time", an_mode = "first", sample_step = 500)
    #plot_chain_progression(a, bounds = (1,10), target = "t", fix_phi = np.pi/2)
    #chain_progression(a, bounds = (1,20), target = "t", show = True, L_ref = True)
    #inspect_chain_first_maxima(b, bounds = (1,30), by = .02)
    #comp_t_chain_progression([a,b,c],mode = "best", target = "t", unit_bounds = (1,30), l_ref = True )
    #comp_best_t_chain_progression([a,b,c], target = "t", unit_bounds = (1,30), l_ref = True )

    #t_chain_progression_phases(a, bounds = (1,40), sample_step = 4, mode = "first").legend()

    #multi_2_phases_example(sample = 200)

    #chain_ch_comp( bounds = (3,50))
    #line_speedup_perf_comp( bounds = (6,60), step = 5, su_bounds = (.1,2.5,1000), mode = "first", TC = 1).legend()

    #a.rephase(-1j)
    #b = qwg.SquareCut()
    #plot_performance(b, an_mode = "first")
    #plot_performance(b, an_mode = "TC")
    #plot_evo_mat_heatmap(b)
    #plot_evo_vs_phase(b ,end = 10)

    #plot_evo_chain_speedup(a, 10, bounds = (0,100), step = .2,su_sample = 300, show = True)
    #plot_evo_line_speedup(10)
    #plot_evo_line_speedup(11)


    #parallel can chase the Bancho optimum
    # a = qwg.Parallel(4,2)
    # a = qwg.chain(a,32)
    # plot_performance(a, mode = "diag")

    b = qwg.Ring(50)
    #b.re_coord[0] = (2,1)
    #plot_performance(b)
    #b.rephase(1j)
    #print(b.mat)

    #b.krylov_basis(mode = "basis_plot")
    #b.krylov_basis(mode = "link_plot")
    #plot_evo_mat_heatmap(b)


    #a.krylov_basis(mode = "basis_plot")
    #a.krylov_basis(mode = "link_plot")

    an = Analyzer(mode = "first", TC = 20)
    #plot_size_progression_multi( bounds = (4,40), step = 2, loglog = True, target = "p", analyzer = an).legend()
    #plot_odd_even_progression( bounds = (3,40), target = "p", analyzer = an).legend()
    #
    plot_chain_progression_multi_loglog(bounds = (5,500), points = 50, target = "p", analyzer = an, fast = True)

    #b = qwg.Ring(3)
    #b = qwg.chain(b, 30)
    #b.re_coord[0] = (2,1)
    #plot_performance(b)
    #b.rephase( np.repeat(-1j, len(b.re_coord)))
    #print(b.mat)

    #b.krylov_basis(mode = "basis_plot")
    #b.krylov_basis(mode = "link_plot")
    #plot_evo_mat_heatmap(b)


    #odd_even_time_lm(HANDLES = False)
    #a = qwg.Ring(4)
    #a = qwg.SquareCut()
    #plot_performance(a, an_mode = "TC")
    #plot_chain_progression(a,bounds = (1,20), target = "t")
    #time_chain_progression_lm(a)

    #t_size_progression_phases( bounds = (4,16), step = 2, analyzer = an).legend()


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
