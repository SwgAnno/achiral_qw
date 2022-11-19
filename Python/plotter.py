from concurrent.futures.thread import _shutdown
from datetime import time
from itertools import chain
from simulator import *
from bessel import *
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
from matplotlib.ticker import MaxNLocator

#todo: decide on the fate of those global variables
TC = 1
qut = False



def plot_evo_mat(gr , start = 0, end = None, by = .1, filter = None, TC = None, ax = None):
    """
    Plot probability evolution on each and every site of the graph
    todo: discuss show buffering option   
    """

    assert end or TC , "plot_evo_mat, no time bounds given"

    if not end:
        end = gr.distance()*TC
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

    #no previous plot given, must create from scatch
    if ax == None:
        fig, ax = plt.subplots()
        
        ax.set_xlabel('Time')
        ax.set_ylabel('Expectation values')

    for i in selection :
        ax.plot(seq, evo[i][:], label = str(i))


    #return ax parameter for further additions
    ax.legend()
    return ax

def plot_evo_mat_heatmap(gr , start = 0, end = None, by = .1, filter = None, TC = None, fig = None, ax = None):
    """
    Display probability evolution on each node of a graph in a 2D heatmap
    """
    
    assert end or TC , "plot_evo_mat_heatmap, no time bounds given"

    if not end:
        end = gr.distance()*TC
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

    #no previous plot
    if ax == None :
        fig, ax = plt.subplots()

    ax.set_xlabel('t')
    ax.set_ylabel('site')
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))

    c = ax.pcolormesh(seq, selection, evo,vmin = 0, vmax = 1,label = gr.code)

    fig.colorbar(c, ax=ax)

    #dPass parameter for further additions
    return ax

def plot_evo_vs_phase(gr , start = 0, end = None, by = .1, phase_by = .1, TC = None, ax = None):
    """
    Plot probability evolution on target site as a function of time and phase
        Graph must have 1 free phase to work properly
    """

    assert end or TC , "plot_evo_vs_phase error, no time bounds given"
    if not end:
        end = gr.distance()*TC
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

    if ax == None :
        fig, ax = plt.subplots()

    c = ax.pcolormesh(seq, phase_seq, data, cmap = "inferno",vmin = 0, vmax = 1, label = gr.code)
    
    ax.set_xlabel('t')
    ax.set_ylabel(chr(952))

    fig.colorbar(c, ax=ax)

    #Pass parameter for further additions
    return ax

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

def plot_evo_vs_derivative(gr, l = 0, start = 0, end = None, by = .1, TC = None, ax = None):
    """
    Comparison between target site probability and its derivative
    """
    assert end or TC , "plot_evo_vs_derivative error, no time bounds given"
    if not end:
        end = gr.distance()*TC
    global qut

    print(qut)
    
    seq = np.arange(start,end,by)
    
    solver = SESolver(gr,qutip = qut)
    evo = solver.target_p(seq)
    deriv = solver.target_p_prime(seq)

    if ax == None :
        fig, ax = plt.subplots()
        ax.set_xlabel('Time')
        ax.set_ylabel('P')

    ax.plot(seq, evo)
    ax.plot(seq, deriv)

    ax.legend(["P", "dP/dt"])
    
    #Pass parameter for firther additions
    return ax

########################################
#plot performance methods

def plot_performance(gr, sample_step = 100, target = "p", mode = None, an_mode = "TC", TC = 1, \
                    ax = None, **kwargs):
    """
    Generic wrapper to plot transport performance as a function of phases
    actually a router method for graph-specific routines       
    """
    global qut
    an = Analyzer(gr, TC = TC, qutip = qut, mode = an_mode)

    if mode == "diag":
        return plot_performance_diag(sample_step, target, an, ax, **kwargs)
    elif gr.get_phase_n() == 1 :
        return plot_performance_1(sample_step, target, an, ax, **kwargs)
    elif gr.get_phase_n() == 2 :
        return plot_performance_2(sample_step, target, an, ax, **kwargs)
    elif mode == "time":
        return plot_performance_time(sample_step, an, ax, **kwargs)
    else :
        raise NotImplementedError("Plot mode not found or graph not supported")

def plot_performance_diag(sample_step, target, an, ax = None):
    """
    Transport performance in equal phases setting as a function of one parameter
    """
    seq = np.linspace(0, np.pi*2, sample_step)

    perf = []
    perf = an.performance_full_diag(sample_step = sample_step, target = target)

    if ax == None :
        fig, ax = plt.subplots()

        set_performance_plot(ax, target)

    ax.plot(seq, perf, label = an.get_label())
    
    #Pass parameter for further additions
    return ax

def plot_performance_time( sample_step, an, ax = None):
    """
    Plot diagonal trsnsport probability performance
    + time of arrival information as color
    """
    seq = np.linspace(0, np.pi*2, sample_step)

    perf = []
    perf = an.performance_full_diag(sample_step = sample_step, target = "p")

    time = []
    time = an.performance_full_diag(sample_step = sample_step, target = "t")

    if ax == None :
        fig, ax = plt.subplots()

        set_performance_plot(ax, target = "p")

    ax.scatter(seq, perf, s = 10, c = time , cmap = "viridis")

    norm = colors.Normalize(vmin= min(time), vmax=max(time))
    map = cm.ScalarMappable(norm, cmap="viridis")
    plt.colorbar( map)
    
    #Pass parameter for further additions
    return ax

def plot_performance_1(sample_step, target, an, ax = None):
    """
    Transport performance for 1-phase graphs
    """

    seq = np.linspace(0, np.pi*2, sample_step)

    perf = an.performance_full(sample_step = sample_step, target = target)

    if ax == None :
        fig, ax = plt.subplots( figsize = (6,5))

        set_performance_plot(ax,target)

    ax.plot(seq, perf, label = an.get_label())

    #Pass parameter for further additions
    return ax
    
# 2-phased graph performance
def plot_performance_2(sample_step, target, an, ax = None, verbose = False):
    """
    Transport performance for 2-phase graphs
    """

    if verbose :
        print("plot_performance_2 : ", an.get_label())
    
    seq = np.linspace(0, np.pi*2, sample_step)

    perf = an.performance_full(sample_step = sample_step, target = target)

    if ax == None:
        fig, ax = plt.subplots()

    set_performance_plot(ax, target = "p", dim = 2)
    
    c = ax.pcolormesh(seq, seq, perf, cmap = "inferno", vmin = 0, vmax = 1,label = an.solver.gr.code)

    #Pass parameter for further additions
    return ax
    
############################################
#utils

def set_performance_plot(ax,target = "p", dim = 1):
        
    if dim == 1 :
        ax.set_xlabel( "$\tetha$")
        ax.set_xlim(0,6.28)

        if target == "p":
            ax.set_ylabel('$P_{max}$')
            ax.set_ylim(bottom = 0, top = 1)
        else :
            ax.set_ylabel("t")
            ax.set_ylim(bottom = 0, auto = True)

    elif dim == 2:
        ax.set_xlabel("$\\theta_1$")
        ax.set_ylabel("$\\theta_2$")

def set_progression_plot(ax,x_mode = "dist", target = "p"):

    ax.set_xlabel(x_mode)

    if target == "p":
        ax.set_ylabel('$P_{max}$')
        ax.set_ylim(bottom = 0, top = 1)
    elif target == "t":
        ax.set_ylabel("t")
        ax.set_ylim(bottom = 0, auto = True)
    else :
        raise NotImplementedError("target mode not supported")

#probably usless helper 
def probability_colorbar_map():
    norm = colors.Normalize(vmin= 0, vmax=1)
    map = cm.ScalarMappable(norm, cmap="inferno")
    return map


