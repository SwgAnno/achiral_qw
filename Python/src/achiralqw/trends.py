from achiralqw.simulator import Analyzer
from achiralqw.graph import QWGraphBuilder
from achiralqw.bessel import *
from achiralqw.collection import CollectionBuilder, QWGraphCollection, get_line_data
from scipy.optimize import minimize_scalar
from achiralqw.plotter import set_progression_plot
import matplotlib.pyplot as plt


# plot tick API
from matplotlib.ticker import MaxNLocator


#todo: this shouldnt be able to plot just give the progression
def bessel_progression( bounds = (3,12), target = "p", L_ref = True, approx_ref = True):
    """
    Progression of first maxima of squared bessel function
    """

    def f(x,l):
        return -4* np.power( np.abs( jv(l,2*x)), 2)

    data = np.empty((2, bounds[1]-bounds[0]+1))

    for i in range(bounds[1]-bounds[0] +1):

        order = bounds[0]+i

        res = minimize_scalar( f, bounds = (order/2, order), method="bounded", args = (order) )
        data[:,i] = [ order, res.x]

    print(data)

    if target == "p":
        for i in range( data.shape[1]):
            data[1,i] = -1* f(data[1,i], data[0,i])

    ax = plot_standard_progression(data, target = target, x_mode = "order")

    if L_ref:
        line_data = get_line_data( bounds, target = target, x_mode = "size")
        ax.plot(line_data[0]+1, line_data[1], color = "green")

    if approx_ref:

        def approx(x):
            k1 = .6748851
            k2 = .16172
            k3 = .02918
            max_trend =  k1*np.power(x, -1/3) * (1 -k2*np.power(x,-2/3) + k3*np.power(x,-4/3))
            return 4*max_trend*max_trend
        
        grid = np.arange(bounds[0], bounds[1], .1)
        eval = []
        for t in grid:
            eval.append( approx(t))

        ax.plot(grid, eval)

    plt.show()

def line_progression_vs_bessel_expansion(bounds = (3,10)):
    """
    Compare line behaviour with Abramovich(???) handbook expansion
    """

    x = np.arange(bounds[0], bounds[1]+1)

    def approx(x):
        k1 = .6748851
        k2 = .16172
        k3 = .02918
        max_trend =  k1*np.power(x, -1/3) * (1 -k2*np.power(x,-2/3) + k3*np.power(x,-4/3))
        return 4* max_trend*max_trend


    line_data = get_line_data(  bounds, target = "p", x_mode = "dist")

    ax = plot_standard_progression(line_data, target = "p", x_mode = "dist")
    
    grid = np.arange(bounds[0], bounds[1], .1)
    eval = []
    for t in grid:
        eval.append( approx(t))

    ax.plot(grid, eval)

    plt.show()
    
def plot_speedup_performance(N, target = "p", ax = None) :
    """
    L(N) best transport time/probability as a function of internal speedup
    """

    sample = np.linspace(.1,10, 1000)
    data = np.empty( len(sample))

    cur = QWGraphBuilder.Line(N)
    an = Analyzer(cur, mode = "first")

    for i in range(len(sample)):
        cur = QWGraphBuilder.Line(N, speedup = sample[i])
        an.set_gr(cur)

        data[i] = an.locate_max()[1] if target == "p" else an.locate_max()[0]

    if ax == None:
        fig, ax = plt.subplots()
        set_progression_plot(ax, x_mode = "$\mu$ speedup", target=target)
    
    ax.plot(sample, data, label = "P" + str(N))

    print(max(data), sample[np.argmax(data)])

    #Pass parameter for further additions
    return ax

def plot_evo_line_speedup(N, bounds = (0,50,.5),su_bounds = (.01, 10, 1000), fig = None, ax = None):
    """
    L(N) evolution as a function of time and speedup
    """
    
    gr = QWGraphBuilder.Line(N)
    sample = np.geomspace(*su_bounds)
    t_sample = np.arange(*bounds)
    data = np.empty( ( len(sample), len(t_sample)))

    an = Analyzer(gr, mode = "first")

    for m in range( len(sample)):
        cur = QWGraphBuilder.Line(N, speedup = sample[m])
        an.set_gr(cur)

        data[m,:] = an.evolution_grid(bounds = bounds[0:2], step = bounds[2] )

    if ax == None :
        fig, ax = plt.subplots()
        ax.set_yscale("log")
    
    c = ax.pcolormesh(t_sample, sample, data, label = "LN")
    
    ax.set_xlabel('t')
    ax.set_ylabel('$\mu$')

    fig.colorbar(c, ax=ax)

    #Pass parameter for further additions
    return ax

def plot_evo_chain_speedup(gr, rep, bounds = None, step = .1, su_bounds = (.1, 10), su_sample = 1000, fig = None, ax = None):
    """
    No phase chain graph best transport time/probability as a function of time and speedup
    """
    
    if bounds == None :
        bounds = (0, 50)
    sample = np.linspace(su_bounds[0], su_bounds[1], su_sample)
    t_sample = np.arange(bounds[0], bounds[1], step)
    data = np.empty( ( len(sample), len(t_sample)))

    an = Analyzer(gr.chain( rep), mode = "first")

    for m in range( len(sample)):

        print (m, "of", len(sample))
        cur = gr.chain( rep, speedup = sample[m])
        an.set_gr(cur)
        
        data[m,:] = an.evolution_grid(bounds = bounds, step = step)

    if ax == None :
        fig, ax = plt.subplots()
    
    c = ax.pcolormesh(t_sample, sample, data, label = "LN")
    
    ax.set_xlabel('t')
    ax.set_ylabel('speedup')

    fig.colorbar(c, ax=ax)

    #Pass parameter for further additions
    return ax

def plot_speedup_performance_multi(bounds = (4,20), target = "p",fig = None, ax = None) :

    sample = np.linspace(.1,10, 1000)
    y_sample = np.arange(bounds[0],bounds[1]+1)
    data = np.empty( ( len(y_sample), len(sample)))

    cur = QWGraphBuilder.Line(4)
    an = Analyzer(cur, mode = "first")

    for m in range( len(y_sample)):
        for i in range(len(sample)):
            cur = QWGraphBuilder.Line(y_sample[m], speedup = sample[i])
            an.set_gr(cur)

            data[m,i] = an.locate_max()[1] if target == "p" else an.locate_max()[0]

    if ax == None :
        fig, ax = plt.subplots()
    
    c = ax.pcolormesh(sample, y_sample, data, label = "LN")

    ax.vlines(np.sqrt(2), bounds[0], bounds[1], color = "red")
    
    ax.set_xlabel('speedup')
    ax.set_ylabel('N')

    fig.colorbar(c, ax=ax)

    #Pass parameter for further additions
    return ax

def plot_speedup_performance_multi_chain(gr_unit, bounds = (4,20), target = "p", ax = None) :

    sample = np.linspace(.1,10, 1000)
    y_sample = np.arange(bounds[0],bounds[1]+1)
    data = np.empty( ( len(y_sample), len(sample)))

    cur = QWGraphBuilder.Line(4)
    an = Analyzer(cur, mode = "first", diag = True)

    for m in range( len(y_sample)):
        print("Graph", m, "out of", len(y_sample) )

        #get hypothetical best phase for the given chain
        cur = gr_unit.chain(rep = y_sample[m], speedup = 1)
        an.set_gr(cur)

        best_phi = an.optimum_phase_minimize()[0]
        #print(best_phi)
        
        for i in range(len(sample)):
            cur = gr_unit.chain( rep = y_sample[m], speedup = sample[i])
            an.set_gr(cur)

            data[m,i] = an.performance_diag( phi = best_phi, target = target)

    if ax == None :
        fig, ax = plt.subplots()
    
    c = ax.pcolormesh(sample, y_sample, data, label = "LN")

    ax.vlines(np.sqrt(2), bounds[0], bounds[1], color = "red")
    
    ax.set_xlabel('speedup')
    ax.set_ylabel('N')

    fig.colorbar(c, ax=ax)

    #Pass parameter for further additions
    return ax

#############################################
# progression plots

def plot_standard_progression(prog : QWGraphCollection, target = "p", x_mode = "dist", label = "", ax = None, **kwargs):
    """
    Plot helper for a standard progression ouput
        prog is a QWGraphCollection object
    """
    if ax == None:
        fig, ax = plt.subplots(1, 1, figsize = (6,5))
        set_progression_plot(ax, x_mode=x_mode, target=target)
    
    x, data = prog.evaluate( target = target, x_mode = x_mode, **kwargs)
    ax.plot( x, data, marker = ".", label = label)

    #force integer ticks (discrete families of graphs)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))

    #Pass parameter for further additions
    return ax

def plot_base_progression(g_type, target = "p", \
                            x_mode = "dist", \
                            label = "",ax = None, **kwargs):
    """
    Optimized Transport time/probability for a topology from the 3 basic ones
    """
    if ax == None :
        fig, ax = plt.subplots()

        set_progression_plot(ax, x_mode = x_mode, target = target)

    prog = CollectionBuilder().base_progression(g_type= g_type, **kwargs)
    plot_standard_progression(prog, ax = ax, label = label, target = target)

    #Pass parameter for further additions
    return ax


# chain_progression() data plotting function, adapted for multiple draw call
def plot_chain_progression(gr_unit, target = "p", \
                            x_mode = "dist", \
                            label = "", ax = None, **kwargs):
    """
    Optimized Transport time/probability for a family of chain graphs
    built from a given unit
    """

    if ax == None :
        fig, ax = plt.subplots()

        set_progression_plot(ax, x_mode = x_mode, target = target)

    prog = CollectionBuilder().chain_progression(gr_unit = gr_unit, **kwargs)
    plot_standard_progression(prog, ax = ax, label = label, target = target)
    
    #Pass parameter for further additions
    return ax

def plot_line_data(ax = None, bounds = None, **kwargs):

    if ax == None:
        fig, ax = plt.subplots(1,1,figsize = (6,5))

    if bounds == None:
        lim = ax.get_xlim()
        bounds = (int(max(2, lim[0])), int(lim[1])+1)

    #print(bounds)

    x, data = get_line_data(bounds=bounds,**kwargs)
    ax.plot(x,data, color = "black", marker = ".", label = "P", alpha = .5)