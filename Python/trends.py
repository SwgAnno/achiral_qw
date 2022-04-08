from plotter import *
from simulator import *
from Graph import *
from scipy.optimize import minimize_scalar
import scipy.stats as stats

#progression of first maxima of squared bessel function
def bessel_progression( bounds = (3,12), target = "p", L_ref = True,show = True):

    def f(x,l):
        return -4* np.power( np.abs( jv(l,2*x)), 2)

    data = np.empty((2, bounds[1]-bounds[0]+1))

    for i in range(bounds[1]-bounds[0] +1):

        order = bounds[0]+i

        res = minimize_scalar( f, bounds = (order/2, 3/4*order), method="bounded", args = (order) )
        data[:,i] = [ order, res.x]

    print(data)

    if target == "p":
        for i in range( data.shape[1]):
            data[1,i] = -1* f(data[1,i], data[0,i])

    if L_ref:
        line_data = get_line_data( bounds, target = target, x_mode = "size")

        fig, ax = plot_standard_progression(data, target = target, x_mode = "order", show = False)
        ax.plot(line_data[0]+1, line_data[1], color = "green")

        plt.show()
    else:
        plot_standard_progression(data, target = target, x_mode = "order", show = True)


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

#simple wrapper for L class references
def get_line_data(bounds = (3,10), target = "p", x_mode = "dist"):
    return size_progression("L", bounds = bounds, target = target, x_mode = x_mode)
    
def chain_progression( gr_unit = QWGraph.Ring(4), bounds = (1,10), target = "p", \
                        x_mode = "dist", HANDLES = True, L_ref = True, show = False):
    gr_list = []

    for i in range( bounds[1]- bounds[0] +1):
        gr_list.append( QWGraph.chain(gr_unit, bounds[0]+ i, HANDLES = HANDLES))

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

    

def optimized_progression( g_list, target = "p", mode = "first", diag = True):
    perf = np.empty( len(g_list))
    for i in range(len(g_list)):
        
        print( g_list[i].code)
        tester = Analyzer(g_list[i], mode = mode, qutip = False)
        best_phi = tester.optimum_phase_minimize(diag = diag)[0]

        target_t = (target != "p")

        if diag:
            perf[i] = tester.performance_diag(best_phi, t = target_t)
        else:
            perf[i] = tester.performance(best_phi, t = target_t)
            
        print( i, perf[i])

    x = np.arange(1, len(g_list)+1)

    return x , perf

#plot best transport performance for a general collection of graph
def plot_standard_progression(prog, target = "p", x_mode = "dist", show = False):

    fig, ax = plt.subplots()
    
    ax.plot(prog[0], prog[1])

    if(target == "p"):
        plt.ylim(0,1.1)
            
        ax.set_xlabel(x_mode)
        ax.set_ylabel('max P')
    if(target == "t"):
        ax.set_xlabel(x_mode)
        ax.set_ylabel('t')

    #display finished plot or pass parameter for further additions
    if show:
        plt.show()
    else :
        return fig, ax
    
