from concurrent.futures.thread import _shutdown
from datetime import time
from itertools import chain
from simulator import *
from bessel import *
from trends import *
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
from matplotlib.ticker import MaxNLocator

#todo: decide on the fate of those global variables
TC = 1
qut = False


#Plot probability evolution on each and every site of the graph
#todo: discuss show buffering option
def plot_evo_mat(gr , start = 0, end = None, by = .1, filter = None, TC = None, show = True):

    if not end and not TC:
        print( "plot_evo_mat error, no time bounds given")
    elif not end:
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
def plot_evo_mat_heatmap(gr , start = 0, end = None, by = .1, filter = None, TC = None, show = True):
    
    if not end and not TC:
        print( "plot_evo_mat_heatmap error, no time bounds given")
    elif not end:
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

    fig, ax = plt.subplots()

    c = ax.pcolormesh(seq, selection, evo,label = gr.code)

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

def plot_evo_vs_phase(gr , start = 0, end = None, by = .1, phase_by = .1, TC = None, show = True):
    
    if not end and not TC:
        print( "plot_evo_vs_phase error, no time bounds given")
    elif not end:
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

    fig, ax = plt.subplots()

    c = ax.pcolormesh(seq, phase_seq, data, cmap = "inferno", label = gr.code)
    
    ax.set_xlabel('t')
    ax.set_ylabel(chr(952))

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
def plot_evo_vs_derivative(gr, l = 0, start = 0, end = None, by = .1, TC = None, show = True):

    if not end and not TC:
        print( "plot_evo_vs_derivative error, no time bounds given")
    elif not end:
        end = gr.distance()*TC
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

# phase-dependent best transport maxima
# can select either its value or the time of arrival
# actually a router method for more specific performance-like routines
def plot_performance(gr, sample_step = 100, target = "p", mode = None, an_mode = "TC", TC = 1, \
                     show = True, fig = None, ax = None):

    global qut
    an = Analyzer(gr, TC = TC, qutip = qut, mode = an_mode)

    if mode == "diag":
        return plot_performance_diag(sample_step, target, an, show, fig, ax)
    elif mode == "time":
        return plot_performance_time(sample_step, an, show, fig, ax)
    elif gr.get_phase_n() == 1 :
        return plot_performance_1(sample_step, target, an, show,  fig, ax)
    elif gr.get_phase_n() == 2 :
        return plot_performance_2(sample_step, target, an, show)
    else :
        print("")

# wrapper for plot performance for a list imput
def plot_performance_list(gr_list, sample_step = 100, target = "p", mode = None, an_mode = "TC", TC = 1, \
                     show = True):

    fig = None
    ax = None
    for i in range( len(gr_list) -1):
        fig, ax = plot_performance(gr_list[i], sample_step, target, mode, an_mode, TC,show = False, fig = fig, ax = ax)

    if show :
        plot_performance(gr_list[-1], sample_step, target, mode, an_mode, TC, \
                        show = True, fig = fig, ax = ax)
    else :
        return plot_performance(gr_list[-1], sample_step, target, mode, an_mode, TC, \
                                show = False, fig = fig, ax = ax)

# equal-phases setting best transport maxima
def plot_performance_diag(sample_step, target, an, show , fig = None , ax = None):
    seq = np.linspace(0, np.pi*2, sample_step)

    perf = []
    perf = an.performance_full_diag(sample_step = sample_step, target = target)

    if fig == None or ax == None :
        fig, ax = plt.subplots()

        ax.set_xlabel( chr(952))
        ax.set_xlim(0,6.28)

        if target == "p":
            ax.set_ylabel('max P')
            ax.set_ylim(bottom = 0, top = 1)
        else :
            ax.set_ylabel("t")
            ax.set_ylim(bottom = 0)

    ax.plot(seq, perf, label = an.get_label())
    
    #display finished plot or pass parameter for further additions
    if show:
        ax.legend()
        plt.show()
    else :
        return fig, ax

# diagonal performance + time of arrival as color

def plot_performance_time( sample_step, an, show , fig = None, ax = None):
    seq = np.linspace(0, np.pi*2, sample_step)

    perf = []
    perf = an.performance_full_diag(sample_step = sample_step, target = "p")

    time = []
    time = an.performance_full_diag(sample_step = sample_step, target = "t")

    if fig == None or ax == None :
        fig, ax = plt.subplots()

        ax.set_xlabel( chr(952))
        ax.set_xlim(0,6.28)

        ax.set_ylabel('max P')
        ax.set_ylim(bottom = 0, top = 1)

    ax.scatter(seq, perf, s = 10, c = time , cmap = "viridis")

    norm = colors.Normalize(vmin= min(time), vmax=max(time))
    map = cm.ScalarMappable(norm, cmap="viridis")
    plt.colorbar( map)
    
    #display finished plot or pass parameter for further additions
    if show:
        ax.legend()
        plt.show()
    else :
        return fig, ax


# 1-phased graph performance
def plot_performance_1(sample_step, target, an, show, fig = None, ax = None):
    seq = np.linspace(0, np.pi*2, sample_step)

    perf = an.performance_full(sample_step = sample_step, target = target)

    if fig == None or ax == None :
        fig, ax = plt.subplots()

        ax.set_xlabel( chr(952))
        ax.set_xlim(0,6.28)

        if target == "p":
            ax.set_ylabel('max P')
            ax.set_ylim(bottom = 0, top = 1)
        else :
            ax.set_ylabel("t")
            ax.set_ylim(bottom = 0)
    
    ax.plot(seq, perf, label = an.get_label())


    #display finished plot or pass parameter for further additions
    if show:
        ax.legend()
        plt.show()
    else :
        return fig, ax
    
# 2-phased graph performance
def plot_performance_2(sample_step, target, an, show):
    seq = np.linspace(0, np.pi*2, sample_step)

    perf = an.performance_full(sample_step = sample_step, target = target)

    fig, ax = plt.subplots()
    
    c = ax.pcolormesh(seq, seq, perf, cmap = "inferno", label = an.solver.gr.code)
    
    ax.set_xlabel(chr(952) + '1')
    ax.set_ylabel(chr(952) + '2')

    fig.colorbar(c, ax=ax)

    #display finished plot or pass parameter for further additions
    if show:
        plt.show()
    else :
        return fig, ax

# chain_progression() data plotting function, adapted for multiple draw call
def plot_chain_progression(gr_unit, bounds = (1,10), target = "p", \
                            x_mode = "dist", HANDLES = True, fix_phi = None, \
                            show = True, fig = None, ax = None):
    
    x, data = chain_progression( gr_unit = gr_unit, bounds = bounds, target = target, \
                                x_mode = x_mode, HANDLES = HANDLES, fix_phi = fix_phi,\
                                show = False)

    if fig == None or ax == None :
        fig, ax = plt.subplots()

        ax.set_xlabel( "#unit")
        if target == "p":
            ax.set_ylabel('max P')
            ax.set_ylim(bottom = 0, top = 1)

        else :
            ax.set_ylabel("t")
            #ax.set_ylim(bottom = 0, top = 'auto')

    phi_label = ""
    if fix_phi != None:
        phi_label = str(int( fix_phi // (np.pi/6))) + "/6" + chr(960)

    
    ax.plot(x, data, label = gr_unit.code + " " + phi_label)
    
    #display finished plot or pass parameter for further additions
    if show:
        ax.legend()
        plt.show()
    else :
        return fig, ax
    
################



