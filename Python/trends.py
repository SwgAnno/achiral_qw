from plotter import *
from simulator import *
from Graph import *
import scipy.stats as stats

def size_progression(g_type = "C", bounds = (3,12), target = "p", show = False):
    
    gr_list = []
    for i in range( bounds[1]- bounds[0] +1):
        if      g_type == "C":
            gr_list.append( QWGraph.Ring( bounds[0]+ i) )
        elif    g_type == "Ch":
            gr_list.append( QWGraph.Ring( bounds[0]+ i, HANDLES = True) )
        elif    g_type == "L":
            gr_list.append( QWGraph.Line( bounds[0]+ i))

    out = optimized_progression(gr_list, target = target)
    x = np.arange( bounds[0], bounds[1]+1)
    if show:
        plot_standard_progression([x, out[1]], target = target, show = True)
    else:
        return out
    
def chain_progression( gr_unit = QWGraph.Ring(4), bounds = (1,10),target = "p", HANDLES = True, show = False):
    gr_list = []

    for i in range( bounds[1]- bounds[0] +1):
        gr_list.append( QWGraph.chain(gr_unit, bounds[0]+ i))

    out = optimized_progression(gr_list, target = target)
    x = np.arange( bounds[0], bounds[1]+1)
    if show:
        plot_standard_progression([x, out[1]], target = target, show = True)
    else:
        return out
    
def time_size_progression_lm(g_type = "C"):
    
    gr_list = []

    #carefully chosen(?)
    bounds = [5,20]
    
    for i in range( bounds[1]- bounds[0] +1):
        if      g_type == "C":
            gr_list.append( QWGraph.Ring( bounds[0]+ i) )
        elif    g_type == "Ch":
            gr_list.append( QWGraph.Ring( bounds[0]+ i, HANDLES = True) )
        elif    g_type == "L":
            gr_list.append( QWGraph.Line( bounds[0]+ i))

    data = optimized_progression(gr_list, target = "t")

    #print(data)
    x = np.arange( bounds[0], bounds[1]+1)
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
def plot_standard_progression(prog, target = "p", show = False):

    fig, ax = plt.subplots()
    
    ax.plot(prog[0], prog[1])

    if(target == "p"):
        plt.ylim(0,1)
            
        ax.set_xlabel('n')
        ax.set_ylabel('max P')
    if(target == "t"):
        ax.set_xlabel('n')
        ax.set_ylabel('t')

    #display finished plot or pass parameter for further additions
    if show:
        plt.show()
    else :
        return fig, ax
    
