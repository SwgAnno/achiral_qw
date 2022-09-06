from simulator import *
from bessel import *
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

#todo: decide on the fate of those global variables
TC = 1
qut = False


#Plot probability evolution on each and every site of the graph
#todo: discuss show buffering option
def plot_evo_mat(gr , start = 0, end = None, by = .1, filter = None, show = True):

    if not end:
        global TC
        end = gr.N * TC
    global qut

    seq = np.arange(start,end,by)
    
    solver = SESolver(gr, qutip = qut)
    evo = solver.evo_p_psi( gr.get_start_state(qut), seq)

    print("Grid-wise maximum: ", max(evo[gr.target][:]))

    if filter == None :
        selection = range(gr.N)
    elif filter == "target" :
        selection = [gr.target]
    elif filter == "start" :
        selection = [gr.start]
    else :
        if type(filter) == int :
            selection = [filter]
        else :
            selection = filter

    fig, ax = plt.subplots()
    for i in selection :
        ax.plot(seq, evo[i][:], label = str(i))
    ax.set_xlabel('Time')
    ax.set_ylabel('Expectation values')


    #display finished plot or pass parameter for firther additions
    if show:
        ax.legend()
        plt.show()
    else :
        return fig, ax

#display evolution on a 2D heatmap
def plot_evo_mat_heatmap(gr , start = 0, end = None, by = .1, filter = None, show = True):
    
    if not end:
        global TC
        end = gr.N * TC
    global qut

    seq = np.arange(start,end,by)
    
    solver = SESolver(gr, qutip = qut)
    evo = solver.evo_p_psi( gr.get_start_state(qut), seq)

    print("Grid-wise maximum: ", max(evo[gr.target][:]))

    if filter == None :
        selection = np.arange(0, gr.N)
    elif filter == "target" :
        selection = [gr.target]
    elif filter == "start" :
        selection = [gr.start]
    else :
        if type(filter) == int :
            selection = [filter]
        else :
            selection = filter

    fig, ax = plt.subplots()

    c = ax.pcolormesh(seq, selection, evo, label = gr.code)

    fig.colorbar(c, ax=ax)

    ax.set_xlabel('t')
    ax.set_ylabel('site')
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))


    #display finished plot or pass parameter for firther additions
    if show:
        ax.legend()
        plt.show()
    else :
        return fig, ax
    pass

def plot_evo_vs_phase(gr , start = 0, end = None, by = .1, phase_by = .1, show = True):
    
    if not end:
        global TC
        end = gr.N * TC
    global qut

    seq = np.arange(start,end,by)
    phase_seq = np.arange(0,2*np.pi, phase_by)

    phase_vec = np.exp( 1j * phase_seq)

    data = np.ndarray( (len(phase_seq), len(seq)))
    
    an = Analyzer(gr, qutip = qut)

    for i in range( len(phase_vec)):
        an.rephase_gr( np.repeat( phase_vec[i], \
                                  an.dim() ))

        data[i:] = an.evo_full( bounds = (start,end), step = by)

    fig, ax = plt.subplots()

    c = ax.pcolormesh(seq, phase_seq, data, label = gr.code)
    
    ax.set_xlabel('t')
    ax.set_ylabel('phase')

    fig.colorbar(c, ax=ax)


    #display finished plot or pass parameter for firther additions
    if show:
        ax.legend()
        plt.show()
    else :
        return fig, ax
    pass

def plot_line_bessel_evolution(l = 5, end = 30):

    grid = np.arange(0,end, .1)
    eval = []
    for i in range(l):
        eval.append( [])

    for t in grid:
        for i in range(l):
            eval[i].append( line_evo_bessel(l-1,i,t))

    fig, ax = plt.subplots()

    for i in range(l):
        ax.plot(grid,eval[i], label = str(i))

    ax.legend()

    plt.show()

def plot_ring_bessel_evolution(l = 5, end = 30):

    grid = np.arange(0,end, .1)
    eval = []
    for i in range(l):
        eval.append( [])

    for t in grid:
        for i in range(l):
            eval[i].append( ring_evo_bessel(l,i,t))

    fig, ax = plt.subplots()

    for i in range(l):
        ax.plot(grid,eval[i], label = str(i))

    ax.legend()

    plt.show()

#compare line evolution on target site with analytical dynamics with bessel functions
#todo: make it work :(
def plot_line_vs_bessel(n = 5, l= None, end = 30, trace_conn = False):

    if l == None:
        l = n

    gr = QWGraph.Line(n+1)
    if trace_conn:
        gr.retrace_conn()

    fig,ax = plot_evo_mat(gr, end = end, filter = l, show = False)

    grid = np.arange(0,end, .1)
    eval = []
    for t in grid:
        eval.append( line_evo_bessel(n,l,t))

    ax.plot(grid,eval, label = "Bessel")
    ax.legend()

    plt.show()

#compare ring evolution on target site with analytical dynamics with bessel functions
def plot_ring_vs_bessel(l = 5, end = 30):

    gr = QWGraph.Ring(l)

    fig,ax = plot_evo_mat(gr, end = end,filter = "target", show = False)

    grid = np.arange(0,end, .1)
    eval = []
    for t in grid:
        eval.append( ring_evo_bessel(l,None,t)/2)

    ax.plot(grid,eval, label = "Bessel")
    ax.legend()

    plt.show()

#comparison between target site rpobability and its derivative
def plot_evo_vs_derivative(gr, l = 0, start = 0, end = None, by = .1, show = True):

    if not end:
        global TC
        end = gr.N * TC
    global qut

    print(qut)
    
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

    #display finished plot or pass parameter for further additions
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



