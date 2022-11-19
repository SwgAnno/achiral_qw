#list of methonds to draw specific plots in the article

from trends import *
from plotter import *
from matplotlib.pyplot import figure

def odd_even_time_lm( L_ref = True, HANDLES = False):
    """
    Compute the two linear model from the transport time trends of,
    possibly handled, odd and even C graphs
    """

    coeff = []

    print("Even lm:")
    gr_list = []
    for i in range(16):
        gr_list.append( QWGraph.Ring( 4+ i*2, HANDLES = HANDLES))

    coeff.append (time_progression_lm(gr_list, x_mode = "dist") )

    print("Odd lm:")
    gr_list = []
    for i in range(16):
        gr_list.append( QWGraph.Ring( 5+ i*2, HANDLES = HANDLES))

    coeff.append( time_progression_lm(gr_list, x_mode = "dist") )

    if not L_ref :
        return

    print("Line lm:")
    gr_list = []
    for i in range(21):
        gr_list.append( QWGraph.Line( 10+ i))

    coeff.append( time_progression_lm(gr_list, x_mode = "dist") )
    #print( get_line_data( bounds = (4,20), target = "t"))

    print( "even/line m : \t" , coeff[0][0]/coeff[2][0])
    print( "odd/line m : \t" , coeff[1][0]/coeff[2][0])

def plot_performance_odd( step = 100):
    """
    Example plot of transport probability for ODD C graphs
    """

    gr_list = []

    gr_list.append( QWGraph.Ring(5, HANDLES = True))
    gr_list.append( QWGraph.Ring(7))
    gr_list.append( QWGraph.Ring(9, HANDLES = True))
    gr_list.append( QWGraph.Ring(11))

    fig, ax = plt.subplots(1,1, figsize = (6,5))

    for gr in gr_list:
        plot_performance(gr, sample_step = step , an_mode="first", ax = ax)

    ax.legend()
    plt.show()

def plot_performance_even( step = 100):
    """
    Example plot of transport probability for EVEN C graphs
    """
    gr_list = []

    gr_list.append( QWGraph.Ring(4))
    gr_list.append( QWGraph.Ring(6, HANDLES = True))
    gr_list.append( QWGraph.Ring(8))
    gr_list.append( QWGraph.Ring(10, HANDLES = True))

    fig, ax = plt.subplots()

    for gr in gr_list:
        plot_performance(gr, sample_step = step , an_mode="first", ax = ax)

    ax.legend()
    plt.show()

def plot_performance_1_multi(gr,TC_vec, first = False, step = 100, ax = None):
    """
    Ad hoc plot performance wrapper for multiple  TC
    """
    if ax == None :
        fig, ax = plt.subplots(1,1, figsize = (6,5))

    if first:
        plot_performance(gr, an_mode = "first", sample_step = step, ax = ax)

    for TC in TC_vec:
        plot_performance(gr, an_mode= "TC", TC = TC, sample_step = step, \
                         ax = ax)

    return ax 

def t_chain_progression_multi(gr_unit = QWGraph.Ring(3), bounds = (1,10), sample_step = 5,l_ref = True,\
                             ax = None):
    """
    Transport time performance for a given family of chain plotted for a range of fixed phase values
    """
    if gr_unit.code[0] == "C" and gr_unit.N % 2 == 0:
        phi_vec = np.arange(0, sample_step) * np.pi/sample_step
    else :
        phi_vec = np.arange(0, sample_step) * 2*np.pi/sample_step
    print(phi_vec)

    if ax == None :
        fig, ax = plt.subplots(1,1, figsize = (6,5))

    for phase in phi_vec:
        plot_chain_progression(gr_unit = gr_unit, bounds = bounds, target = "t", fix_phi = phase, \
                               ax = ax)

    if l_ref :
        a = ax.get_xlim()
        l_bounds = (2 ,int(a[1]))
        L_x, L_data = get_line_data( l_bounds, target = "t")
        ax.plot(L_x, L_data, label = "L", color = "green")

    return ax

def t_size_progression_multi(g_type = "C", bounds = (3,15), sample_step = 5, l_ref = False, ax = None, **kwargs):
    """
    Transport time performance for a given family from a selected topologies,
    plotted for a range of fixed phase values
    """

    if g_type == "C" and bounds[0] % 2 == 0:
        phi_vec = np.arange(0, sample_step) * np.pi/sample_step
    else :
        phi_vec = np.arange(0, sample_step) * 2*np.pi/sample_step
    print(phi_vec)

    gr_list = []
    for i in range( bounds[1]- bounds[0] +1):
        gr_list.append( QWGraph.Ring( bounds[0]+ i*2))

    if ax == None:
        fig,ax = plt.subplots(1,1, figsize = (6,5) )

    for phase in phi_vec :

        data = optimized_progression(gr_list,target = "t", opt_mode = "fix", opt_phi = phase, **kwargs)
        plot_standard_progression(data ,target = "t",ax = ax)

    if l_ref:    
        a = ax.get_xlim()
        l_bounds = (2 ,int(a[1]))
        L_x, L_data = get_line_data( l_bounds, target = "t")
        ax.plot(L_x, L_data, label = "L", color = "green")

    return ax

def plot_size_progression_multi( bounds = (3,12), target = "p") :
    """
    Example plot of progression for the 3 base topologies together
    """

    fig, ax = plt.subplots(1,1, figsize = (6,5))

    plot_size_progression("L", bounds = bounds, target = target, ax = ax)
    plot_size_progression("C", bounds = bounds, target = target, ax = ax)
    plot_size_progression("Ch", bounds = bounds, target = target, ax = ax)

    ax.legend()
    plt.show()

def plot_chain_progression_multi( bounds = (3,20), target = "p"):
    """
    Example plot of progression for the 3 best chain graph families
    """

    #bounds[0] has to be >=3
    chain3_bounds = (bounds[0]-2, bounds[1]-3)
    chain4_bounds = ( (bounds[0]-2)//2, (bounds[1]-3)//2)

    fig, ax = plt.subplots(1,1, figsize = (6,5))

    plot_chain_progression( QWGraph.Ring(3), bounds = chain3_bounds, target = target, ax = ax)
    plot_chain_progression( QWGraph.Ring(4), bounds = chain4_bounds, target = target, ax = ax)
    plot_chain_progression( QWGraph.SquareCut(), bounds = chain3_bounds, target = target, ax = ax)
    plot_size_progression( "L", bounds = bounds, target = target, ax = ax)

    ax.legend()
    plt.show()

def plot_odd_even_progression( bounds = (3,12), target = "p", x_mode = "dist",mode = "TC", TC = 1) :
    """
    Show differences in behaviour between odd and even C graph
    You can plot either best transport time or probability
    """

    fig, ax = plt.subplots(1,1, figsize = (6,5))

    gr_list = []
    for i in range( bounds[1]- bounds[0] +1):
        gr_list.append( QWGraph.Line( bounds[0]+ i))
    
    out = optimized_progression(gr_list, target = target, mode = mode, TC = TC, opt_mode = "smart")
    x = get_list_x(gr_list, x_mode = x_mode)

    plot_standard_progression([x, out[1]], target = target, x_mode = x_mode, label = "L", ax = ax)

    label = []
    if bounds[0]%2 == 0 :
        label.append("C even")
        label.append("C odd")
    else:
        label.append("C odd")
        label.append("C even")

    gr_list = []
    for i in range( bounds[1]- bounds[0] +1):
        gr_list.append( QWGraph.Ring( bounds[0]+ i*2))
    
    out = optimized_progression(gr_list, target = target, mode = mode, TC = TC, opt_mode = "smart")
    x = get_list_x(gr_list, x_mode = x_mode)

    plot_standard_progression([x, out[1]], target = target, x_mode = x_mode, label = label[0], ax = ax)

    gr_list = []
    for i in range( bounds[1]- bounds[0] +1):
        gr_list.append( QWGraph.Ring( bounds[0]+1+ i*2))
    
    out = optimized_progression(gr_list, target = target, mode = mode, TC = TC, opt_mode = "smart")
    x = get_list_x(gr_list, x_mode = x_mode)

    plot_standard_progression([x, out[1]], target = target, x_mode = x_mode, label = label[1], ax = ax)

    ax.legend()
    plt.show()

def chain_ch_comp( bounds = (3,20), target = "p"):
    """
    Example plot with progression comparison between C3 chain and Ch topology
    """

    chain3_bounds = (bounds[0]-2, bounds[1]-3)

    fix, ax = plt.subplots(1,1, figsize = (6,5))

    plot_chain_progression( QWGraph.Ring(3), bounds = chain3_bounds, target = target, ax = ax)
    plot_size_progression( "Ch", bounds = bounds, target = target, ax = ax)
    plot_size_progression( "L", bounds = bounds, target = target, ax = ax)

    ax.legend()
    plt.show()
    
def chain_performance_multi_speedup( gr_unit,su_vec, rep = 10, sample_step = 100, target = "p"):
    """
    Plot performance wrapper to compare performance plots for a chain built with different speedups
    """

    fig, ax = plt.subplots(1,1, figsize = (6,5))

    for i in range(len(su_vec)):
        print("Graph", i+1, "out of", len(su_vec) )
        chain = QWGraph.chain(gr_unit, rep, speedup = su_vec[i])
        chain.code = "C3chain su" + str(su_vec[i])
        plot_performance( chain, sample_step= sample_step, mode = "diag",ax = ax)

    ax.legend()
    plt.show()

def multi_2_phases_example( sample = 100):
    """
    Example plot with 2phase performance comparison between first maxima and TC10 search for a sample of graphs
    """

    gr_num = 3
    gr_list = []

    gr_list.append( QWGraph.Ring(9)| QWGraph.Ring(4))
    gr_list.append( QWGraph.Ring(8)+ QWGraph.Ring(8))
    gr_list.append( QWGraph.SquareCut())

    fig, axx = plt.subplots(2,gr_num, figsize = (9,6))

    for i in range(len(gr_list)):
        plot_performance( gr_list[i],sample_step=sample, an_mode = "first", ax = axx[0][i], verbose = True)
        plot_performance( gr_list[i],sample_step=sample, an_mode = "TC", TC = 10, ax = axx[1][i],  verbose = True)


    # add single colorbar trick
    # https://stackoverflow.com/questions/13784201/how-to-have-one-colorbar-for-all-subplots

    fig.subplots_adjust(.05, hspace = .25)
    cbar_ax = fig.add_axes([0.92, 0.1, 0.025, 0.8])
    map = probability_colorbar_map()
    fig.colorbar(map, cax=cbar_ax)

    plt.show()




