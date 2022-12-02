from simulator import *
from Graph import *
from bessel import *
from scipy.optimize import minimize_scalar
from plotter import set_progression_plot
import scipy.stats as stats

# plot tick API
from matplotlib.ticker import MaxNLocator

def unit_list(bounds, unit):

    b_0 = max(0, bounds[0]// unit.distance())

    return np.arange(b_0, bounds[1]// unit.distance())


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


def size_progression(g_type = "C", bounds = (3,12), target = "p", x_mode = "dist", speedup = None, L_ref = False, show = False,  **kwargs):
    """
    Transport Time/Probability of a family of graph from a selected list of topologies
        suppoerted topologies: C, Ch, L
    """
    
    gr_list = []
    for i in range( bounds[1]- bounds[0] +1):
        if      g_type == "C":
            gr_list.append( QWGraph.Ring( bounds[0]+ i*2) )
        elif    g_type == "Ch":
            gr_list.append( QWGraph.Ring( bounds[0]+ i*2, HANDLES = True) )
        elif    g_type == "L":
            gr_list.append( QWGraph.Line( bounds[0]+ i, speedup = speedup))

    out = optimized_progression(gr_list, target = target, **kwargs)
    x = get_list_x(gr_list, x_mode = x_mode)

    #todo: show option on data methods is bullshit
    if show:
        if L_ref:
            line_data = get_line_data( (min(x), max(x)), target = target, x_mode = x_mode)

            ax = plot_standard_progression([x, out[1]], target = target, x_mode = x_mode)
            ax.plot(line_data[0], line_data[1], color = "green")

        else :
            plot_standard_progression([x, out[1]], target = target, x_mode = x_mode)

        plt.show()
    else:
        return x, out[1]

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
    pass



def get_line_data(bounds = (3,10), target = "p", x_mode = "dist", **kwargs):
    """
    Simple wrapper for L graphs progression references   
    """
    return size_progression("L", bounds = bounds, target = target, x_mode = x_mode)

def plot_line_data(ax = None, bounds = None, **kwargs):

    if ax == None:
        fig, ax = plt.subplots(1,1,figsize = (6,5))

    if bounds == None:
        lim = ax.get_xlim()
        bounds = (int(max(2, lim[0])), int(lim[1])+1)

    print(bounds)

    x, data = get_line_data(bounds=bounds,**kwargs)
    ax.plot(x,data, color = "green", label = "L")
    
def chain_progression( gr_unit = QWGraph.Ring(4), bounds = (1,10), target = "p", \
                        x_mode = "dist", HANDLES = True, fix_phi = None):
    """
    Transport Time/Probability of a family of chain graph from a given unit
    """

    gr_list = []
    sizes = unit_list(bounds,gr_unit)

    for size in sizes:
        gr_list.append( QWGraph.chain(gr_unit,size, HANDLES = HANDLES))

    if fix_phi != None:
        out = optimized_progression(gr_list, target = target, opt_mode = "fix", opt_phi = fix_phi)
    else:
        out = optimized_progression(gr_list, target = target)
        
    x = get_list_x(gr_list, x_mode = x_mode)

    return x, out[1]

def time_progression_lm( gr_list, x_mode = "dist"):
    """
    Construct a linear model of the best transport time from a generic list of graphs
    """
    
    data = optimized_progression(gr_list, target = "t")

    #print(data)
    x = get_list_x(gr_list, x_mode = x_mode)
    out = stats.linregress(x,data[1])

    print("m: ", out.slope, " +- ", out.stderr)
    print("q: ", out.intercept, " +- ", out.intercept_stderr)
    print("r: ", out.rvalue)

    return out.slope, out.intercept

def time_size_progression_lm(g_type = "C", x_mode = "dist"):
    """
    Construct a linear model of the best transport time for a selected famoly of graphs
        supported topologies: C,Ch, L (relies on size_progression)
    """
    #todo: should it rely on time_progression?
    gr_list = []

    #carefully chosen(?)
    bounds = [5,20]
    
    data = size_progression(g_type, bounds, target = "t")

    #print(data)
    x = get_list_x(gr_list, x_mode = x_mode)
    out = stats.linregress(x,data[1])

    print("m: ", out.slope, " +- ", out.stderr)
    print("q: ", out.intercept, " +- ", out.intercept_stderr)
    print("r: ", out.rvalue)

    return out.slope, out.intercept

def time_chain_progression_lm(gr_unit = QWGraph.Ring(4), x_mode = "dist", HANDLES = True, L_ref = True):
    """
    Construct a linear model of the best transport time from a family of chain graphs
    built from a given unit
    """
    
    gr_list = []

    #carefully chosen(?)
    bounds = [1,10]
    
    data = chain_progression(gr_unit, bounds, HANDLES = HANDLES, target = "t")

    out = stats.linregress(data[0],data[1])

    print("m: ", out.slope, " +- ", out.stderr)
    print("q: ", out.intercept, " +- ", out.intercept_stderr)
    print("r: ", out.rvalue)

    if not L_ref :
        return out.slope, out.intercept

    print("Line lm:")
    gr_list = []
    for i in range(21):
        gr_list.append( QWGraph.Line( 10+ i))

    l_coeff = time_progression_lm(gr_list, x_mode = "dist")
    #print( get_line_data( bounds = (4,20), target = "t"))

    print( gr_unit.code + "chain/line m : \t" , out.slope/l_coeff[0])

    return out.slope, out.intercept

    

def optimized_progression( g_list, target = "p", mode = "first",TC = None, diag = True, opt_mode = None, opt_phi = None):
    """
    Return best trasnport time/probability from a given list of graphs
    """
    
    perf = np.empty( len(g_list))

    print("########\nOptimized progression", target, mode)
    for i in range(len(g_list)):
        
        print( i,"of", len(g_list),  "//",g_list[i].code)
        tester = Analyzer(g_list[i], mode = mode, TC = TC, qutip = False)
        best_phi = 0

        if opt_mode == "smart" :
            best_phi = tester.optimum_phase_smart()[0]
        elif opt_mode == "fix" :
            best_phi = opt_phi
        else :
            best_phi = tester.optimum_phase_minimize(diag = diag)[0]

        #print(best_phi, opt_mode)
        target_t = (target != "p")

        if diag:
            perf[i] = tester.performance_diag(best_phi, t = target_t)
        else:
            perf[i] = tester.performance(best_phi, t = target_t)
            
        #print( i, perf[i])

    x = np.arange(1, len(g_list)+1)

    return x , perf
    
def plot_speedup_performance(N, target = "p", ax = None) :
    """
    L(N) best transport time/probability as a function of internal speedup
    """

    sample = np.linspace(.1,10, 1000)
    data = np.empty( len(sample))

    cur = QWGraph.Line(N)
    an = Analyzer(cur, mode = "first")

    for i in range(len(sample)):
        cur = QWGraph.Line(N, speedup = sample[i])
        an.set_gr(cur)

        data[i] = an.locate_max()[1] if target == "p" else an.locate_max()[0]

    if ax == None:
        fig, ax = plt.subplots()
        set_progression_plot(ax, x_mode = "$\mu$ speedup", target=target)
    
    ax.plot(sample, data, label = "L" + str(N))

    print(max(data), sample[np.argmax(data)])

    #Pass parameter for further additions
    return ax

def plot_evo_line_speedup(N, bounds = (0,50,.5),su_bounds = (.01, 10, 1000), fig = None, ax = None):
    """
    L(N) evolution as a function of time and speedup
    """
    
    gr = QWGraph.Line(N)
    sample = np.geomspace(*su_bounds)
    t_sample = np.arange(*bounds)
    data = np.empty( ( len(sample), len(t_sample)))

    an = Analyzer(gr, mode = "first")

    for m in range( len(sample)):
        cur = QWGraph.Line(N, speedup = sample[m])
        an.set_gr(cur)

        data[m,:] = an.evo_full(bounds = bounds[0:2], step = bounds[2] )

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

    an = Analyzer(QWGraph.chain(gr, rep), mode = "first")

    for m in range( len(sample)):

        print (m, "of", len(sample))
        cur = QWGraph.chain(gr, rep, speedup = sample[m])
        an.set_gr(cur)
        
        data[m,:] = an.evo_full(bounds = bounds, step = step)

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

    cur = QWGraph.Line(4)
    an = Analyzer(cur, mode = "first")

    for m in range( len(y_sample)):
        for i in range(len(sample)):
            cur = QWGraph.Line(y_sample[m], speedup = sample[i])
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

    cur = QWGraph.Line(4)
    an = Analyzer(cur, mode = "first")

    for m in range( len(y_sample)):
        print("Graph", m, "out of", len(y_sample) )

        #get hypothetical best phase for the given chain
        cur = QWGraph.chain(gr_unit, rep = y_sample[m], speedup = 1)
        an.set_gr(cur)

        best_phi = an.optimum_phase_minimize( diag= True)[0]
        #print(best_phi)
        
        for i in range(len(sample)):
            cur = QWGraph.chain(gr_unit, rep = y_sample[m], speedup = sample[i])
            an.set_gr(cur)

            data[m,i] = an.performance_diag( phi = best_phi, t = (target !="p"))

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

def plot_standard_progression(prog, target = "p", x_mode = "dist", label = "", ax = None):
    """
    Plot helper for a standard progression ouput
        prog[0] is a lsist containing x_mode indices
        prog[1] is a list containing the relative performances
    """
    if ax == None:
        fig, ax = plt.subplots(1, 1, figsize = (6,5))
        set_progression_plot(ax, x_mode=x_mode, target=target)
    
    ax.plot(prog[0], prog[1], marker = ".", label = label)

    #force integer ticks (discrete families of graphs)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))

    #Pass parameter for further additions
    return ax

# chain_progression() data plotting function, adapted for multiple draw call
def plot_chain_progression(gr_unit, target = "p", \
                            x_mode = "dist", fix_phi = None, \
                            ax = None, **kwargs):
    """
    Optimized Transport time/probability for a family of chain graphs
    built from a given unit
    """
    
    x, data = chain_progression( gr_unit = gr_unit, target = target, \
                                x_mode = x_mode, fix_phi = fix_phi,\
                                **kwargs)

    if ax == None :
        fig, ax = plt.subplots()

        set_progression_plot(ax, x_mode = x_mode, target = target)

    if fix_phi != None:
        phi_label =  ("%.2f" % (fix_phi / np.pi)) + chr(960)
    else:
        phi_label = ""
    
    ax.plot(x, data, ".-", label = gr_unit.code + " " + phi_label)
    
    #Pass parameter for further additions
    return ax

def plot_size_progression(g_type = "C", target = "p", x_mode = "dist", speedup = None, \
                         ax = None, **kwargs) :

    """
    Optimized Transport time/probability for a given family of graphs
    Supported topologies: C, Ch, L (relies on size prorgession)
    """
    prog_data = size_progression( g_type = g_type, target = target, x_mode = x_mode, **kwargs)

    plot_standard_progression( prog_data, target = target, x_mode = x_mode, label = g_type, ax = ax)

    #Pass parameter for further additions
    return ax