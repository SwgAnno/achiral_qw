from simulator import *
import matplotlib.pyplot as plt

#todo: decide on the fate of those global variables
TC = 1
qut = False


#Plot probability evolution on each and every site of the graph
#todo: discuss show buffering option
def plot_evo_mat(gr , start = 0, end = None, by = .1, show = True):

    if not end:
        global TC
        end = gr.N * TC
    global qut

    seq = np.arange(start,end,by)
    
    solver = SESolver(gr, qutip = qut)
    evo = solver.evo_p_psi( gr.get_start_state(qut), seq)

    print("Grid-wise maximum: ", max(evo[gr.target][:]))

    fig, ax = plt.subplots()
    for i in range(gr.N) :
        ax.plot(seq, evo[i][:])
    ax.set_xlabel('Time')
    ax.set_ylabel('Expectation values')
    
    ax.legend(list(range(0,gr.N)))

    #display finished plot or pass parameter for firther additions
    if show:
        plt.show()
    else :
        return fig, ax

#comparison between target site rpobability and its derivative
def plot_evo_vs_derivative(gr, l = 0, start = 0, end = None, by = .1, show = True):

    if not end:
        global TC
        end = gr.N * TC
    global qut
    
    seq = np.arange(start,end,by)
    
    solver = SESolver(gr,qutip = qut)
    evo = solver.target_p(seq)
    deriv = solver.target_p_prime(seq)

    fig, ax = plt.subplots()

    ax.plot(seq, evo)
    ax.plot(seq, deriv)
    
    ax.set_xlabel('Time')
    ax.set_ylabel('P')

    ax.legend(["P", "dP/dt"])
    
    #display finished plot or pass parameter for firther additions
    if show:
        plt.show()
    else :
        return fig, ax

# phase-dependent best transport maximua
# actually a router method for more specific performance-like routines
def plot_performance(gr, sample_step = 100, mode = None, an_mode = "TC", show = True):

    global TC
    global qut
    an = Analyzer(gr, TC = TC, qutip = qut, mode = an_mode)

    if mode == "diag":
        plot_performance_diag(sample_step, an, show)
    elif gr.get_phase_n() == 1 :
        plot_performance_1(sample_step, an, show)
    elif gr.get_phase_n() == 2 :
        plot_performance_2(sample_step, an, show)
    else :
        print("")
        
# equal-phases setting best transport maxima
def plot_performance_diag(sample_step, an, show ):
    seq = np.linspace(0, np.pi*2, sample_step)

    perf = []
    perf = an.performance_full_diag(sample_step)

    fig, ax = plt.subplots()
    
    ax.plot(seq, perf)
    
    ax.set_xlabel('phi')
    ax.set_ylabel('max P')

    #display finished plot or pass parameter for firther additions
    if show:
        plt.show()
    else :
        return fig, ax
    
# 1-phased graph performance
def plot_performance_1(sample_step, an, show):
    seq = np.linspace(0, np.pi*2, sample_step)

    perf = an.performance_full(sample_step)

    fig, ax = plt.subplots()
    
    ax.plot(seq, perf)
    
    ax.set_xlabel('phi')
    ax.set_ylabel('max P')

    #display finished plot or pass parameter for firther additions
    if show:
        plt.show()
    else :
        return fig, ax
    
# 2-phased graph performance
def plot_performance_2(sample_step, an, show):
    seq = np.linspace(0, np.pi*2, sample_step)

    perf = an.performance_full(sample_step)

    fig, ax = plt.subplots()
    
    c = ax.pcolormesh(seq, seq, perf)
    
    ax.set_xlabel('phi 1')
    ax.set_ylabel('phi 2')

    fig.colorbar(c, ax=ax)

    #display finished plot or pass parameter for firther additions
    if show:
        plt.show()
    else :
        return fig, ax
    
################



