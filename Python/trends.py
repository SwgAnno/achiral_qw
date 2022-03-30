from plotter import *
from simulator import *
from Graph import *

def size_progression(g_type = "C", bounds = (3,12)):
    
    gr_list = []
    for i in range( bounds[1]- bounds[0] +1):
        if      g_type == "C":
            gr_list.append( QWGraph.Ring( bounds[0]+ i) )
        elif    g_type == "Ch":
            gr_list.append( QWGraph.Ring( bounds[0]+ i, HANDLES = True) )
        elif    g_type == "L":
            gr_list.append( QWGraph.Line( bounds[0]+ i))

    plot_optimized_list(gr_list, x = np.arange( bounds[0], bounds[1]+1), show = True)

def chain_progression( gr_unit = QWGraph.Ring(4), bounds = (1,10), HANDLES = True):
    gr_list = []

    for i in range( bounds[1]- bounds[0] +1):
        gr_list.append( QWGraph.chain(gr_unit, bounds[0]+ i))

    plot_optimized_list(gr_list, x = np.arange( bounds[0], bounds[1]+1), show = True)

#plot best transport performance for a general collection of graph
def plot_optimized_list(g_list, x = None, mode = "first", show = False):

    perf = np.empty( len(g_list))
    for i in range(len(g_list)):
        
        print( g_list[i].code)
        tester = Analyzer(g_list[i], mode = mode, qutip = False)
        perf[i] = tester.optimum_phase_minimize(diag = True)[1]
        print( i, perf[i])

    fig, ax = plt.subplots()

    plt.ylim(0,1)

    if not x.any():
        x = np.arange(1, len(g_list)+1)
    
    ax.plot(x, perf)
    
    ax.set_xlabel('phi')
    ax.set_ylabel('max P')

    #display finished plot or pass parameter for firther additions
    if show:
        plt.show()
    else :
        return fig, ax
    
