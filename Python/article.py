#list of methonds to draw specific plots in the article

from trends import *
from plotter import *
from collection import *
from matplotlib.pyplot import figure

def odd_even_time_lm( L_ref = True, HANDLES = False):
    """
    Compute the two linear model from the transport time trends of,
    possibly handled, odd and even C graphs
    """

    coeff = []
    an = Analyzer( mode = "first")

    print("Even lm:")

    prog = CollectionBuilder().C_progression( bounds = (4,38), step = 2, analyzer= an)
    coeff.append (prog.transport_time_lm() )

    print("Odd lm:")
    prog = CollectionBuilder().C_progression( bounds = (5,39), step = 2, analyzer= an)
    coeff.append (prog.transport_time_lm() )

    if not L_ref :
        return

    print("Line lm:")
    prog = CollectionBuilder().C_progression( bounds = (10,31), step = 1, analyzer= an)
    coeff.append (prog.transport_time_lm() )

    #print( get_line_data( bounds = (4,20), target = "t"))

    print( "even/line m : \t" , coeff[0][0]/coeff[2][0])
    print( "odd/line m : \t" , coeff[1][0]/coeff[2][0])

def plot_performance_odd_even( step = 100):
    """
    Example plot of transport probability for ODD vs EVEN C graphs
    """

    fig, axx = plt.subplots(1,2, figsize = (10,4))

    gr_list_odd = []

    gr_list_odd.append( QWGraph.Ring(5, HANDLES = True))
    gr_list_odd.append( QWGraph.Ring(7))
    gr_list_odd.append( QWGraph.Ring(9, HANDLES = True))
    gr_list_odd.append( QWGraph.Ring(11))

    gr_list_even = []

    gr_list_even.append( QWGraph.Ring(4))
    gr_list_even.append( QWGraph.Ring(6, HANDLES = True))
    gr_list_even.append( QWGraph.Ring(8))
    gr_list_even.append( QWGraph.Ring(10, HANDLES = True))

    for gr_o, gr_e in zip(gr_list_odd, gr_list_even):
        plot_performance(gr_o, sample_step = step , an_mode="first", ax = axx[0])
        plot_performance(gr_e, sample_step = step , an_mode="first", ax = axx[1])

    axx[0].legend()
    axx[1].legend()
    plt.show()

def comp_evo_vs_phase(gr, TC_vec, **kwargs):

    n = len(TC_vec)

    fig, axx = plt.subplots(1,n, figsize = (5*n,4)) #5*n definitely is not scalable but what did you expect

    for ax, TC in zip(axx,TC_vec) :
        plot_evo_vs_phase(gr, TC = TC, ax = ax, **kwargs )
        ax.set_title("$\\nu$ ={}".format(TC))

    fit_colorbar(fig)
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

def comp_performance_multi( gr_list, **kwargs):

    fig, axx = plt.subplots(1, len(gr_list), figsize = (11,4))

    for ax , gr in zip(axx, gr_list):
        print("Performance_multi for {}".format(gr.code))
        plot_performance_1_multi(gr,ax = ax, **kwargs)
        ax.legend()

    plt.show()

def t_chain_progression_phases(gr_unit = QWGraph.Ring(3), bounds = (1,10), sample_step = 5,l_ref = True,\
                             ax = None, **kwargs):
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

    an = Analyzer(opt_mode = "fix", **kwargs)
    prog = CollectionBuilder().chain_progression(gr_unit=gr_unit, bounds = bounds, analyzer=an)
    label = gr_unit.code + " {:1.2f}$\pi$"

    for phase in phi_vec:
        an.set_fix_phi(phase)
        plot_standard_progression(prog, target = "t", \
                                    label = label.format(phase/np.pi), ax =ax)

    if l_ref :
        a = ax.get_xlim()
        l_bounds = (2 ,int(a[1]))
        plot_line_data( bounds = l_bounds , target = "t", ax = ax)

    return ax

def comp_t_chain_progression(gr_list, mode = "best", l_ref = True, **kwargs):
    
    fig, axx = plt.subplots(1, len(gr_list), figsize = (10,4))

    for ax , gr in zip(axx, gr_list):
        print("Performance_multi for {}".format(gr.code))
        if mode == "best":
            plot_chain_progression(gr, ax = ax, **kwargs)
        else:
            t_chain_progression_phases(gr,ax = ax, **kwargs)
        ax.set_title(gr.code)

    max_t = max( [ax.get_ylim()[1] for ax in axx])
    max_s = max( [ax.get_xlim()[1] for ax in axx])
    for ax in axx:
        ax.set_ylim(0,max_t)
        ax.set_xlim(right = max_s)

        if l_ref:
            plot_line_data(ax,**kwargs)

        ax.legend()

    plt.show()

def comp_best_t_chain_progression(gr_list, l_ref = True, ax = None, **kwargs):

    if ax == None:
        fig , ax = plt.subplots(1,1, figsize = (6,5))

    for gr in gr_list:
        print("Performance_best for {}".format(gr.code))
        plot_chain_progression(gr, ax = ax, **kwargs)

    if l_ref:
        plot_line_data(ax,**kwargs)

    ax.legend()

    return ax

def t_size_progression_phases(g_type = "C", bounds = (3,15), sample_step = 5, l_ref = False, ax = None, **kwargs):
    """
    Transport time performance for a given family from a selected topologies,
    plotted for a range of fixed phase values
    """

    if g_type == "C" and bounds[0] % 2 == 0:
        phi_vec = np.arange(0, sample_step) * np.pi/sample_step
    else :
        phi_vec = np.arange(0, sample_step) * 2*np.pi/sample_step
    print(phi_vec)

    if ax == None:
        fig,ax = plt.subplots(1,1, figsize = (6,5) )
    prog = CollectionBuilder().base_progression(g_type, bounds = bounds, **kwargs)
    label = g_type + " {:1.2f}$\pi$"

    for phase in phi_vec:
        plot_standard_progression(prog, target = "t", opt_mode = "fix", opt_phi = phase, \
                                    label = label.format(phase/np.pi), ax =ax)

    if l_ref:    
        a = ax.get_xlim()
        l_bounds = (2 ,int(a[1]))
        L_x, L_data = get_line_data( l_bounds, target = "t")
        ax.plot(L_x, L_data, label = "P", color = "green")

    return ax

def plot_size_progression_multi( ax = None, loglog = False, **kwargs) :
    """
    Example plot of progression for the 3 base topologies together
    """

    if ax == None:
        fig, ax = plt.subplots(1,1, figsize = (6,5))
    
    if loglog:
        ax.set_xscale("log")
        ax.set_yscale("log")

    plot_base_progression("C", ax = ax, **kwargs, label = "C")
    plot_base_progression("Ch", ax = ax, **kwargs, label = "Ch")
    plot_line_data(ax = ax)

    return ax

def comp_size_progression(TC_vec = [1,5,20], **kwargs):

    fig, axx = plt.subplots(1, len(TC_vec), figsize = (10,4))

    for ax , TC in zip(axx, TC_vec):
        plot_size_progression_multi(ax = ax, mode = "TC", TC = TC, **kwargs)
    
    axx[-1].legend()

    plt.show()



def plot_chain_progression_multi( bounds = (3,20), target = "p", analyzer = None):
    """
    Example plot of progression for the 3 best chain graph families
    """

    fig, ax = plt.subplots(1,1, figsize = (6,5))

    set_progression_plot(ax, x_mode = "dist", target="p")
    ax.set_ylim(0.1,1)

    plot_chain_progression( QWGraph.Ring(3), bounds = bounds, target = target, ax = ax, analyzer = analyzer)
    plot_chain_progression( QWGraph.Ring(4),bounds = bounds, target = target, ax = ax, analyzer = analyzer)
    plot_chain_progression( QWGraph.SquareCut(), bounds = bounds, target = target, ax = ax, analyzer = analyzer)
    plot_line_data(target = target, ax = ax)

    ax.legend()
    
    return fig, ax

def plot_chain_progression_multi_loglog( bounds = (3,20), points = 50, target = "p", analyzer = None):


    fig, ax = plt.subplots(1,1, figsize = (6,5))
    ax.set_xscale("log")
    ax.set_yscale("log")

    select = np.geomspace(*bounds, num = points, dtype=int)
    select = set(select)
    select = [x for x in select]
    select.sort()
    select = np.array(select)

    print(select)
    
    plot_chain_progression( QWGraph.Ring(3)    , select = select, target = target, ax = ax, analyzer = analyzer, label = "C3")
    plot_chain_progression( QWGraph.Ring(4)    , select = select, target = target, ax = ax, analyzer = analyzer, label = "C4")
    plot_chain_progression( QWGraph.SquareCut(), select = select, target = target, ax = ax, analyzer = analyzer, label = "DiC4")
    plot_base_progression(  "P", select = select, target = target, label = "P", ax = ax)

    ax.legend()
    
    return fig, ax

def plot_odd_even_progression( bounds = (3,12),target = "p", ax = None, **kwargs) :
    """
    Show differences in behaviour between odd and even C graph
    You can plot either best transport time or probability
    """

    if ax == None:
        fig, ax = plt.subplots(1,1, figsize = (6,5))

    label = []
    if bounds[0]%2 == 0 :
        label.append("C even")
        label.append("C odd")
    else:
        label.append("C odd")
        label.append("C even")

    prog = CollectionBuilder().C_progression(bounds = bounds, step = 2, **kwargs)

    plot_standard_progression(prog, target = target,label = label[0], ax = ax)
    
    prog = CollectionBuilder().C_progression(bounds = (bounds[0]+1, bounds[1]+1), step = 2, **kwargs)

    plot_standard_progression(prog, target = target,label = label[0], ax = ax)

    plot_line_data(target = target, ax = ax)

    return ax

def chain_ch_comp( bounds = (3,20), target = "p"):
    """
    Example plot with progression comparison between C3 chain and Ch topology
    """

    chain3_bounds = (bounds[0]-2, bounds[1]-3)

    fix, ax = plt.subplots(1,1, figsize = (6,5))

    plot_chain_progression( QWGraph.Ring(3), bounds = chain3_bounds, target = target, ax = ax)
    plot_base_progression( "Ch", bounds = bounds, target = target, ax = ax)
    plot_base_progression( "P", bounds = bounds, target = target, ax = ax)

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
    fit_colorbar(fig)

    plt.show()


def line_speedup_perf_comp( bounds = (4,20),step = 3, su_bounds = (.1,2, 1000), target = "p", ax = None , **kwargs):
    """
    Best transport performance of a sample L(N) graphs as a function of speedup
    **kwargs are passed to the analyzer
    """
    
    sample = np.linspace(*su_bounds)
    y_sample = np.arange(start = bounds[0], stop = bounds[1]+1, step = step)
    print(y_sample)
    data = np.empty( ( len(y_sample), len(sample)))
    labels = []

    cur = QWGraph.Line(4)
    an = Analyzer(cur, **kwargs)

    for m in range( len(y_sample)):
        print("Speedup iteration {:.2f}%".format((m+1)/len(y_sample)*100), end='\r')
        for i in range(len(sample)):
            cur = QWGraph.Line(y_sample[m], speedup = sample[i])
            an.set_gr(cur)

            data[m,i] = an.locate_max()[1] if target == "p" else an.locate_max()[0]
        labels.append(an.get_label( mode = False))

    if ax == None :
        fig, ax = plt.subplots()

        ax.set_ylim(0,1)

    for i in range( len(y_sample)):
        ax.plot(sample, data[i], label = labels[i])

    ax.vlines(np.sqrt(2), 0,1,linestyle='dashed', color = "red")
    
    ax.set_xlabel('$\mu$')
    ax.set_ylabel('$P_{max}$')

    #Pass parameter for further additions
    return ax




