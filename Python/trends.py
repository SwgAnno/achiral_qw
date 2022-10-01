from plotter import *
from simulator import *
from Graph import *
from scipy.optimize import minimize_scalar
import scipy.stats as stats

#progression of first maxima of squared bessel function
def bessel_progression( bounds = (3,12), target = "p", L_ref = True, approx_ref = True,show = True):

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

    fig, ax = plot_standard_progression(data, target = target, x_mode = "order", show = False)

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

        plt.plot(grid, eval)

    plt.show()


def size_progression(g_type = "C", bounds = (3,12), target = "p", x_mode = "dist", speedup = None, L_ref = False, show = False):
    
    gr_list = []
    for i in range( bounds[1]- bounds[0] +1):
        if      g_type == "C":
            gr_list.append( QWGraph.Ring( bounds[0]+ i) )
        elif    g_type == "Ch":
            gr_list.append( QWGraph.Ring( bounds[0]+ i, HANDLES = True) )
        elif    g_type == "L":
            gr_list.append( QWGraph.Line( bounds[0]+ i, speedup = speedup))

    out = optimized_progression(gr_list, target = target)
    x = get_list_x(gr_list, x_mode = x_mode)

    if show:
        if L_ref:
            line_data = get_line_data( (min(x), max(x)), target = target, x_mode = x_mode)

            fig, ax = plot_standard_progression([x, out[1]], target = target, x_mode = x_mode, show = False)
            ax.plot(line_data[0], line_data[1], color = "green")

            plt.show()
        else :
            fig, ax = plot_standard_progression([x, out[1]], target = target, x_mode = x_mode, show = True)
    else:
        return out

#Compare line behaviour with Abramovich(???) handbook expansion
def line_progression_vs_bessel_expansion(bounds = (3,10)):

    x = np.arange(bounds[0], bounds[1]+1)

    def approx(x):
        k1 = .6748851
        k2 = .16172
        k3 = .02918
        max_trend =  k1*np.power(x, -1/3) * (1 -k2*np.power(x,-2/3) + k3*np.power(x,-4/3))
        return 4* max_trend*max_trend


    line_data = get_line_data(  bounds, target = "p", x_mode = "dist")

    fig, ax = plot_standard_progression(line_data, target = "p", x_mode = "dist", show = False)
    
    grid = np.arange(bounds[0], bounds[1], .1)
    eval = []
    for t in grid:
        eval.append( approx(t))

    plt.plot(grid, eval)

    plt.show()
    pass


#simple wrapper for L class references
def get_line_data(bounds = (3,10), target = "p", x_mode = "dist"):
    return size_progression("L", bounds = bounds, target = target, x_mode = x_mode)
    
def chain_progression( gr_unit = QWGraph.Ring(4), bounds = (1,10), target = "p", \
                        x_mode = "dist", HANDLES = True, L_ref = True, fix_phi = None, show = False):
    gr_list = []

    for i in range( bounds[1]- bounds[0] +1):
        gr_list.append( QWGraph.chain(gr_unit, bounds[0]+ i, HANDLES = HANDLES))

    if fix_phi != None:
        out = optimized_progression(gr_list, target = target, opt_mode = "fix", opt_phi = fix_phi)
    else:
        out = optimized_progression(gr_list, target = target)
        
    x = get_list_x(gr_list, x_mode = x_mode)

    if show:
        if L_ref:
            fig, ax = plot_standard_progression([x, out[1]], target = target, x_mode = x_mode, show = False)

            line_data = get_line_data( (min(x), max(x)), target = target, x_mode = x_mode)
            ax.plot(line_data[0], line_data[1], color = "green")

        else :
            fig, ax = plot_standard_progression([x, out[1]], target = target, x_mode = x_mode, show = True)
            
        plt.show()
    else:
        return x, out[1]

def time_progression_lm( gr_list, x_mode = "dist"):
    
    data = optimized_progression(gr_list, target = "t")

    #print(data)
    x = get_list_x(gr_list, x_mode = x_mode)
    out = stats.linregress(x,data[1])

    print("m: ", out.slope, " +- ", out.stderr)
    print("q: ", out.intercept, " +- ", out.intercept_stderr)
    print("r: ", out.rvalue)

    return out.slope, out.intercept

def time_size_progression_lm(g_type = "C", x_mode = "dist"):
    
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

def time_chain_progression_lm(gr_unit = QWGraph.Ring(4), x_mode = "dist", HANDLES = True):
    
    gr_list = []

    #carefully chosen(?)
    bounds = [1,10]
    
    data = chain_progression(gr_unit, bounds, HANDLES = HANDLES, target = "t")

    #print(data)
    x = get_list_x(gr_list, x_mode = x_mode)
    out = stats.linregress(x,data[1])

    print("m: ", out.slope, " +- ", out.stderr)
    print("q: ", out.intercept, " +- ", out.intercept_stderr)
    print("r: ", out.rvalue)

    return out.slope, out.intercept

    

def optimized_progression( g_list, target = "p", mode = "first",TC = None, diag = True, opt_mode = None, opt_phi = None):
    perf = np.empty( len(g_list))
    for i in range(len(g_list)):
        
        print( g_list[i].code)
        tester = Analyzer(g_list[i], mode = mode, TC = TC, qutip = False)
        best_phi = 0

        if opt_mode == "smart" :
            best_phi = tester.optimum_phase_smart()[0]
        elif opt_mode == "fix" :
            best_phi = opt_phi
        else :
            best_phi = tester.optimum_phase_minimize(diag = diag)[0]

        target_t = (target != "p")

        if diag:
            perf[i] = tester.performance_diag(best_phi, t = target_t)
        else:
            perf[i] = tester.performance(best_phi, t = target_t)
            
        #print( i, perf[i])

    x = np.arange(1, len(g_list)+1)

    return x , perf

#plot best transport performance for a general collection of graph
def plot_standard_progression(prog, target = "p", x_mode = "dist", label = "",show = False, fig = None, ax = None):

    if fig == None:
        fig, ax = plt.subplots()
    
    ax.plot(prog[0], prog[1], marker = ".", label = label)

    if(target == "p"):
        plt.ylim(0,1.1)
            
        ax.set_xlabel(x_mode)
        ax.set_ylabel('max P')
    if(target == "t"):
        ax.set_xlabel(x_mode)
        ax.set_ylabel('t')

    #display finished plot or pass parameter for further additions
    if show:
        ax.legend()
        plt.show()
    else :
        return fig, ax
    
def plot_speedup_performance(N, target = "p", show = False) :

    sample = np.linspace(.1,10, 1000)
    data = np.empty( len(sample))

    cur = QWGraph.Line(N)
    an = Analyzer(cur, mode = "first")

    for i in range(len(sample)):
        cur = QWGraph.Line(N, speedup = sample[i])
        an.set_gr(cur)

        data[i] = an.locate_max()[1] if target == "p" else an.locate_max()[0]

    fig, ax = plt.subplots()
    
    ax.plot(sample, data, label = "L" + str(N))

    if(target == "p"):
        plt.ylim(0,1.1)
            
        ax.set_xlabel("speedup")
        ax.set_ylabel('max P')
    if(target == "t"):
        ax.set_xlabel("speedup")
        ax.set_ylabel('t')

    print(max(data), sample[np.argmax(data)])

    #display finished plot or pass parameter for further additions
    if show:
        plt.show()
        ax.legend()
    else :
        return fig, ax

def plot_evo_line_speedup(N, bounds = (0,50), step = .5, show = False):
    
    gr = QWGraph.Line(N)
    sample = np.linspace(.1,10, 1000)
    t_sample = np.arange(bounds[0], bounds[1], step)
    data = np.empty( ( len(sample), len(t_sample)))

    an = Analyzer(gr, mode = "first")

    for m in range( len(sample)):
        cur = QWGraph.Line(N, speedup = sample[m])
        an.set_gr(cur)

        data[m,:] = an.evo_full(bounds = bounds, step = step)

    fig, ax = plt.subplots()
    
    c = ax.pcolormesh(t_sample, sample, data, label = "LN")
    
    ax.set_xlabel('t')
    ax.set_ylabel('speedup')

    fig.colorbar(c, ax=ax)

    #display finished plot or pass parameter for further additions
    if show:
        plt.show()
        ax.legend()
    else :
        return fig, ax

def plot_evo_chain_speedup(gr, rep, bounds = None, step = .1, su_bounds = (.1, 10), su_sample = 1000, show = False):
    
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

    fig, ax = plt.subplots()
    
    c = ax.pcolormesh(t_sample, sample, data, label = "LN")
    
    ax.set_xlabel('t')
    ax.set_ylabel('speedup')

    fig.colorbar(c, ax=ax)

    #display finished plot or pass parameter for further additions
    if show:
        plt.show()
        ax.legend()
    else :
        return fig, ax

def plot_speedup_performance_multi(bounds = (4,20), target = "p", show = False) :

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

    fig, ax = plt.subplots()
    
    c = ax.pcolormesh(sample, y_sample, data, label = "LN")

    ax.vlines(np.sqrt(2), bounds[0], bounds[1], color = "red")
    
    ax.set_xlabel('speedup')
    ax.set_ylabel('N')

    fig.colorbar(c, ax=ax)

    #display finished plot or pass parameter for further additions
    if show:
        plt.show()
        ax.legend()
    else :
        return fig, ax