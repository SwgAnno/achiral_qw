from achiralqw.simulator import Analyzer, SESolver, EigenSESolver, QutipSESolver
from achiralqw.graph import QWGraph, QWGraphBuilder
import achiralqw.bessel as bessel
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import networkx as nx
import matplotlib.colors as colors
import matplotlib.ticker as ticker

##############################################
# Graph info plotting

def plot_qwgraph(gr: QWGraph, ax = None):
    """
    General graph visualization tool
    """

    ref = gr.to_networkx()
    nx.draw_spring(ref, with_labels = True, ax  = ax)

    #todo: more accurate plotting

        # names = []
        # for i in range(self.N):
        #     names.append(str(i))

        # cols = []
        # for e in ref.es:
        #     if e.tuple in self.re_coord or e.tuple[::-1] in self.re_coord:
        #         cols.append("red")
        #     else :
        #         cols.append("black")

        # v_cols = ["yellow"]* self.N
        # v_cols[self.start] = "green"
        # v_cols[self.target] = "red"

        # ref.vs["label"] = names
        # ref.vs["color"] = v_cols
        # ref.es["color"] = cols
        # ig.plot( ref , layout = ig.Graph.layout_fruchterman_reingold(ref))

def plot_graph_eigenbasis(gr: QWGraph):
    """
    Represent eigenvector basis in the site basis with
    modulus of projection(left plot)
    argument of projection(right plot)
    """

    x_range = np.arange(0,gr.N)
    y_range = np.arange(0,len(gr.eig_vec))

    modulus = np.zeros( ( gr.N, len(gr.eig_vec)) )
    phase = np.zeros( ( gr.N, len(gr.eig_vec)) )

    for j in range(len(gr.eig_val)) :
        modulus[:,j] = np.abs( gr.eig_vec[j])
        phase[:,j] = np.angle( gr.eig_vec[j])

    fig, ax = plt.subplots( nrows = 1, ncols = 2, sharex = True, sharey = True)

    c1 = ax[0].pcolormesh(x_range, y_range, modulus, label = "mod")
    c2 = ax[1].pcolormesh(x_range, y_range, phase, label = "phase")

    fig.colorbar(c1, ax = ax[0])
    fig.colorbar(c2, ax = ax[1])

    ax[0].set_title("Projection modulus")
    ax[1].set_title("Relative phase")
    
    for stuff in ax :    
        stuff.set_xlabel('site')
        stuff.set_ylabel('eig_n')

    return fig, ax

def plot_krylov_basis(gr : QWGraph, add_phase = False, **kwargs):
    """
    Represent krylov basis vector in the site basis with
    modulus of projection(left plot)
    argument of projection(right plot)
    """

    k_basis, k_E, k_A = gr.krylov_basis(mode = "",**kwargs)

    modulus = np.zeros( (len(k_basis),gr.N))
    phase = np.zeros( (len(k_basis),gr.N))

    for j in range(len(k_basis)) :
        modulus[j,:] = np.abs( k_basis[j][:,0])
        phase[j,:] = np.angle( k_basis[j][:,0])

    x_range = np.arange(0,gr.N)
    y_range = np.arange(0,len(k_basis))

    if add_phase:
        fig, axx = plt.subplots( nrows = 1, ncols = 2, sharex = True, sharey = True, figsize=(8, 4))

        c1 = axx[0].pcolormesh(x_range, y_range, modulus, label = "mod")
        c2 = axx[1].pcolormesh(x_range, y_range, phase, label = "phase")

        fig.colorbar(c1, ax = axx[0])
        fig.colorbar(c2, ax = axx[1])

        axx[0].set_title("Projection modulus")
        axx[1].set_title("Relative phase")
    
        for ax in axx :    
            ax.set_xlabel('$|n\\rangle$')
            ax.set_ylabel('$|k_n \\rangle$')

        return fig, axx

    else :
        fig, ax = plt.subplots( 1,1, figsize=(6, 5))

        c1 = ax.pcolormesh(x_range, y_range, modulus, label = "mod")

        fig.colorbar(c1, ax = ax)
        ax.set_xlabel('$|n\\rangle$')
        ax.set_ylabel('$|k_n \\rangle$')

        return fig, ax

def plot_krylov_couplings(gr: QWGraph, ax = None, **kwargs):
    """
    Plot coupling between krylov basis states in a scatter plot
    """

    k_basis, k_E, k_A = gr.krylov_basis(mode = "",**kwargs)
    
    k_A = k_A[1:]
    
    if ax == None:
        fig, ax = plt.subplots(  figsize=(6, 5))

    x = np.arange(1,len(k_A)+1)
    
    ax.scatter(x, k_A)
    
    ax.set_xlabel('n')
    ax.set_ylabel('$\\langle k_n | H | k_{n+1}\\rangle$')

    ax.legend()
    
    return ax

 ####################################
 # Simple evolution plotting

def plot_evo_mat(gr : QWGraph, start = 0, end = None, by = .1, filter = None, TC = None, solver = EigenSESolver(), ax = None):
    """
    Plot probability evolution on each and every site of the graph
    todo: discuss show buffering option   
    """

    assert end or TC , "plot_evo_mat, no time bounds given"

    if not end:
        end = gr.distance()*TC

    seq = np.arange(start,end,by)
    
    evo = solver.evolve_state_p(gr, gr.get_start_state(), seq)

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

def plot_evo_mat_heatmap(gr , start = 0, end = None, by = .1, filter = None, TC = None, fig = None, solver = EigenSESolver(), ax = None):
    """
    Display probability evolution on each node of a graph in a 2D heatmap
    """
    
    assert end or TC , "plot_evo_mat_heatmap, no time bounds given"

    if not end:
        end = gr.distance()*TC

    seq = np.arange(start,end,by)
    
    evo = solver.evolve_state_p(gr,  gr.get_start_state(), seq)

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
    ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))

    c = ax.pcolormesh(seq, selection, evo,vmin = 0, vmax = 1,label = gr.code)

    fig.colorbar(c, ax=ax)

    #dPass parameter for further additions
    return ax

def plot_evo_vs_phase(gr , start = 0, end = None, by = .1, phase_by = .1, TC = None, plot_max = False, ax = None):
    """
    Plot probability evolution on target site as a function of time and phase
        Graph must have 1 free phase to work properly
    """

    assert end or TC , "plot_evo_vs_phase error, no time bounds given"
    if not end:
        end = gr.distance()*TC

    seq = np.arange(start,end,by)
    phase_seq = np.arange(0,2*np.pi, phase_by)

    data = np.ndarray( (len(phase_seq), len(seq)))
    max_data = np.empty( (2,len(phase_seq)))
    
    an = Analyzer(gr, mode = "first")

    for i in range( len(phase_seq)):
        data[i:] = an.evolution_grid( bounds = (start,end), step = by, phase_vec = np.repeat( phase_seq[i],an.dim() ))
        max_data[:,i] = [ an.locate_max()[0], phase_seq[i]]

    if ax == None :
        fig, ax = plt.subplots()

    c = ax.pcolormesh(seq, phase_seq, data, cmap = "inferno",vmin = 0, vmax = 1, label = gr.code)

    if plot_max :
        ax.scatter(max_data[0], max_data[1], c = "green", s= 3)
    
    ax.set_xlim(0, seq[-1])
    ax.set_xlabel('t')
    ax.set_ylabel(chr(952))

    #Pass parameter for further additions
    return ax

def plot_line_bessel_evolution(l = 5, end = 30):

    grid = np.arange(0,end, .1)
    eval = []
    for i in range(l):
        eval.append( [])

    for t in grid:
        for i in range(l):
            eval[i].append( bessel.line_evo_bessel(l-1,i,t))

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
            eval[i].append( bessel.ring_evo_bessel(l,i,t))

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

    gr = QWGraphBuilder.Line(n+1)
    if trace_conn:
        gr.retrace_conn()

    fig,ax = plot_evo_mat(gr, end = end, filter = l, show = False)

    grid = np.arange(0,end, .1)
    eval = []
    for t in grid:
        eval.append( bessel.line_evo_bessel(n,l,t))

    ax.plot(grid,eval, label = "Bessel")
    ax.legend()

    plt.show()

#compare ring evolution on target site with analytical dynamics with bessel functions
def plot_ring_vs_bessel(l = 5, end = 30):

    gr = QWGraphBuilder.Ring(l)

    fig,ax = plot_evo_mat(gr, end = end,filter = "target", show = False)

    grid = np.arange(0,end, .1)
    eval = []
    for t in grid:
        eval.append( bessel.ring_evo_bessel(l,None,t)/2)

    ax.plot(grid,eval, label = "Bessel")
    ax.legend()

    plt.show()

def plot_evo_vs_derivative(gr, l = 0, start = 0, end = None, by = .1, TC = None, solver = EigenSESolver(), ax = None):
    """
    Comparison between target site probability and its derivative
    """
    assert end or TC , "plot_evo_vs_derivative error, no time bounds given"
    if not end:
        end = gr.distance()*TC
    
    seq = np.arange(start,end,by)
    
    evo = solver.evolve_default_p(gr, seq)
    deriv = solver.evolve_default_p_deriv(gr, seq)

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

    an = Analyzer(gr, TC = TC, mode = an_mode)

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
    perf = an.performance_grid_diag(sample_step = sample_step, target = target)

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
    perf = an.performance_grid_diag(sample_step = sample_step, target = "p")

    time = []
    time = an.performance_grid_diag(sample_step = sample_step, target = "t")

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

    perf = an.performance_grid(sample_step = sample_step, target = target)

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

    perf = an.performance_grid(sample_step = sample_step, target = target)

    if ax == None:
        fig, ax = plt.subplots()

    set_performance_plot(ax, target = "p", dim = 2)
    
    c = ax.pcolormesh(seq, seq, perf, cmap = "inferno", vmin = 0, vmax = 1,label = an.get_gr().code)

    #Pass parameter for further additions
    return ax
    
############################################
#utils

def set_performance_plot(ax,target = "p", dim = 1):
        
    if dim == 1 :
        ax.set_xlabel( "$\\theta$")
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

def set_progression_plot(ax,x_mode = "dist", target = "p", loglog = False, bounds = None,  **kwargs):


    if x_mode == "size":
        ax.set_xlabel("N")
    else:
        ax.set_xlabel(x_mode)

    if loglog and bounds is not None:
        ax.set_xticks( np.geomspace(*bounds, num = 10, dtype=int))
        ax.get_xaxis().set_major_formatter(ticker.ScalarFormatter())


    if target == "p":
        ax.set_ylabel('$P_{max}$')
        if loglog :
            ax.set_ylim(top = 1, bottom = None, auto = True)
        else:
            ax.set_ylim(bottom = 0, top = 1)
    elif target == "t":
        ax.set_ylabel("$t^*$")
        ax.set_ylim(bottom = 0, auto = True)
    else :
        raise NotImplementedError("target mode not supported")

#probably usless helper 
def probability_colorbar_map():
    norm = colors.Normalize(vmin= 0, vmax=1)
    map = cm.ScalarMappable(norm, cmap="inferno")
    return map

#trick to sneak in a global colorbar into a figure
def fit_colorbar(fig ):
    fig.subplots_adjust(.05, hspace = .25)
    cbar_ax = fig.add_axes([0.92, 0.1, 0.025, 0.8])
    map = probability_colorbar_map()
    fig.colorbar(map, cax=cbar_ax)


