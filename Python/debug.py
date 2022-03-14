from plotter import *
from simulator import *
from Graph import *


#global parameters for plotter methods
TC = 5
qut = True

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

def check_qutip( test_gr = QWGraph.Ring(6), t_0 = [2], deriv = False):

    solver = SESolver(test_gr, qutip = True)

    for i in range(2):
        if not deriv:
            print( solver.target_p(t_0))
        else :
            print( solver.target_p_prime(t_0))

#todo: cleanup
def random():
    #a  = QWGraph.chain(QWGraph.Ring(3), 10)
    #a = QWGraph.Parallel(3,2)
    #a = QWGraph.chain(a,2)
    ##a = a+a
    a = QWGraph.Ring(6)


    plot_performance(a)

    plot_evo_vs_derivative(a)

    test = Analyzer(a, qutip = True)


    test.mode = "first"
    print(test.locate_max())
    test.mode = "TC"
    print(test.locate_max())


    #plot_evo_vs_derivative(a)
    #plot_evo_vs_qutip(a)
    #plot_evo_mat(a)
    ##plot_performance(a,100, mode = "diag", an_mode = "TC")
    ##a.plot()
    #plot_evo_vs_derivative(a|a+a)

    #print(test.optimum_phase_minimize(diag = True))



#######################################

a = QWGraph.Ring(7)
a = a+a
plot_performance(a, an_mode = "first")
