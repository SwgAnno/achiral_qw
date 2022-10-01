from Graph import *
from plotter import *
from simulator import *
from trends import *
from article import *

#global parameters for plotter methods
TC = 5
qut = False

#confront locate_max results with actual evolution profile
def check_locate_max( test_gr = QWGraph.Ring(6), mode = None, qutip = True, TC = 1):

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
def check_evo_vs_qutip( test_gr = QWGraph.Ring(6), l = 0, start = 0, end = None, by = .1, deriv = True):
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
def check_evo_vs_qutip_scatter( test_gr = QWGraph.Ring(6), l = 0, start = 0, end = None, by = .1, deriv = False):
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
def check_qutip( test_gr = QWGraph.Ring(6), t_0 = [2], deriv = False):

    solver = SESolver(test_gr, qutip = True)

    for i in range(2):
        if not deriv:
            print( solver.target_p(t_0))
        else :
            print( solver.target_p_prime(t_0))

#graphic comparison of optimum phase result
def check_optimum_phase( test_gr = QWGraph.Ring(6), mode = None, an_mode = "first", qutip = True):

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

#todo: cleanup
def random():
    #a  = QWGraph.chain(QWGraph.Ring(3), 10)
    #a = QWGraph.Paralle l(3,2)
    #a = QWGraph.chain(a,2)
    ##a = a+a
    a = QWGraph.Ring(6)


    plot_performance(a)

    plot_evo_vs_derivative(a)



#######################################

#plot_line_vs_bessel(l = 5, trace_conn = False)
#bessel_progression(bounds = (2,50), target = "p", L_ref = True)

#plot_line_bessel_evolution(5,end = 10)
#plot_evo_mat(QWGraph.Ring(8), end = 10)
#check_line_bessel_dist(5)

#plot_speedup_performance_multi(bounds = (4,24), target = "p", show = True)

#a = QWGraph.Ring(9)
#plot_evo_vs_phase(a, TC = 2, phase_by = .01)
#plot_performance(a,TC = 2, target = "t", sample_step = 500)

#a = QWGraph.Ring(6, HANDLES = True)
#plot_performance_1_multi(a, TC_vec = [5,20], first = True, step = 1000)
#plot_performance_odd(step = 100)

#a = QWGraph.Ring(8) + QWGraph.Ring(8)
#a = QWGraph.SquareCut()
#plot_performance(a, an_mode = "TC", TC = 2)

#FIX!!!!!
a = QWGraph.SquareCut()
#a = QWGraph.chain(a, 5)
#plot_performance(a, mode = "time", an_mode = "first", sample_step = 500)
#plot_chain_progression(a, bounds = (1,10), target = "t", fix_phi = np.pi/2)
#chain_progression(a, bounds = (1,20), target = "t", show = True, L_ref = True)
t_chain_progression_multi(a, bounds = (1,50), sample_step = 4)



#a.rephase(-1j)
#b = QWGraph.SquareCut()
#plot_performance(b, an_mode = "first")
#plot_performance(b, an_mode = "TC")
#plot_evo_mat_heatmap(b)
#plot_evo_vs_phase(b ,end = 10)

#plot_evo_chain_speedup(a, 10, bounds = (0,100), step = .2,su_sample = 300, show = True)

#b = QWGraph.chain(b,10, speedup = 1)

#b.krylov_basis(mode = "basis_plot")
#b.krylov_basis(mode = "link_plot")
#plot_evo_mat_heatmap(b)


#a.krylov_basis(mode = "basis_plot")
#a.krylov_basis(mode = "link_plot")


#plot_simpler_topo_progression( bounds = (3,30), mode = "TC", target = "t", TC = 10)
#plot_odd_even_progression( bounds = (3,40), mode = "first", target = "t", TC = 1)

#odd_even_time_lm(HANDLES = False)




######################
#theta as char
print(chr(952))
#####################