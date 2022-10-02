#list of methonds to draw specific plots in the article

from ctypes.wintypes import tagRECT
from pickletools import TAKEN_FROM_ARGUMENT4U
from turtle import st
from trends import *
from plotter import *
from matplotlib.pyplot import figure

def plot_simpler_topo_progression( bounds = (3,12), target = "p", x_mode = "dist",mode = "TC", TC = 1):

    gr_list = []
    for i in range( bounds[1]- bounds[0] +1):
        gr_list.append( QWGraph.Line( bounds[0]+ i))
    
    out = optimized_progression(gr_list, target = target, mode = mode, TC = TC)
    x = get_list_x(gr_list, x_mode = x_mode)

    fig, ax = plot_standard_progression([x, out[1]], target = target, x_mode = x_mode, label = "L", show = False)

    gr_list = []
    for i in range( bounds[1]- bounds[0] +1):
        gr_list.append( QWGraph.Ring( bounds[0]+ i*2))
    
    out = optimized_progression(gr_list, target = target, mode = mode, TC = TC)
    x = get_list_x(gr_list, x_mode = x_mode)

    fig, ax = plot_standard_progression([x, out[1]], target = target, x_mode = x_mode, label = "C", show = False, fig = fig, ax = ax)

    gr_list = []
    for i in range( bounds[1]- bounds[0] +1):
        gr_list.append( QWGraph.Ring( bounds[0]+ i*2, HANDLES = True))
    
    out = optimized_progression(gr_list, target = target, mode = mode, TC = TC)
    x = get_list_x(gr_list, x_mode = x_mode)

    plot_standard_progression([x, out[1]], target = target, x_mode = x_mode, label = "Ch", show = True, fig = fig, ax = ax)

def plot_odd_even_progression( bounds = (3,12), target = "p", x_mode = "dist",mode = "TC", TC = 1) :

    gr_list = []
    for i in range( bounds[1]- bounds[0] +1):
        gr_list.append( QWGraph.Line( bounds[0]+ i))
    
    out = optimized_progression(gr_list, target = target, mode = mode, TC = TC, opt_mode = "smart")
    x = get_list_x(gr_list, x_mode = x_mode)

    fig, ax = plot_standard_progression([x, out[1]], target = target, x_mode = x_mode, label = "L",show = False)

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

    fig, ax = plot_standard_progression([x, out[1]], target = target, x_mode = x_mode, label = label[0], show = False, fig = fig, ax = ax)

    gr_list = []
    for i in range( bounds[1]- bounds[0] +1):
        gr_list.append( QWGraph.Ring( bounds[0]+1+ i*2))
    
    out = optimized_progression(gr_list, target = target, mode = mode, TC = TC, opt_mode = "smart")
    x = get_list_x(gr_list, x_mode = x_mode)

    plot_standard_progression([x, out[1]], target = target, x_mode = x_mode, label = label[1] , show = True, fig = fig, ax = ax)

def odd_even_time_lm( L_ref = True, HANDLES = False):

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

def plot_performance_1_multi(gr,TC_vec, first = False, step = 100):
    
    if first:
        fig, ax = plot_performance(gr, an_mode = "first", sample_step = step, show = False)
    else :
        fig = None
        ax = None

    for i in range(len(TC_vec) -1):
        fig, ax = plot_performance(gr, an_mode= "TC", TC = TC_vec[i], sample_step = step, \
                                    show = False, fig = fig, ax = ax)

    plot_performance(gr, an_mode= "TC", TC = TC_vec[-1], sample_step = step, \
                    show = True, fig = fig, ax = ax)

def plot_performance_odd( step = 100):

    gr_list = []

    gr_list.append( QWGraph.Ring(5, HANDLES = True))
    gr_list.append( QWGraph.Ring(7))
    gr_list.append( QWGraph.Ring(9, HANDLES = True))
    gr_list.append( QWGraph.Ring(11))

    plot_performance_list(gr_list, sample_step = step, an_mode = "first")

def plot_performance_even( step = 100):

    gr_list = []

    gr_list.append( QWGraph.Ring(4))
    gr_list.append( QWGraph.Ring(6, HANDLES = True))
    gr_list.append( QWGraph.Ring(8))
    gr_list.append( QWGraph.Ring(10, HANDLES = True))

    plot_performance_list(gr_list, sample_step = step, an_mode = "first")

def t_chain_progression_multi(gr_unit = QWGraph.Ring(3), bounds = (1,10), sample_step = 5):

    if gr_unit.code[0] == "C" and gr_unit.N % 2 == 0:
        phi_vec = np.arange(0, sample_step) * np.pi/sample_step
    else :
        phi_vec = np.arange(0, sample_step) * 2*np.pi/sample_step
    print(phi_vec)

    fig = None
    ax = None
    for i in range(len(phi_vec)):
        fig, ax = plot_chain_progression(gr_unit = gr_unit, bounds = bounds, target = "t", fix_phi = phi_vec[i], \
                                         show = False, fig = fig, ax = ax)

    a = ax.get_xlim()
    l_bounds = (2 ,int(a[1]))
    L_x, L_data = get_line_data( l_bounds, target = "t")
    ax.plot(L_x, L_data, label = "L", color = "green")

    ax.legend()
    plt.show()

def plot_size_progression_multi( bounds = (3,12), target = "p") :

    fig, ax = plot_size_progression("L", bounds = bounds, target = target, show = False)
    fig, ax = plot_size_progression("C", bounds = bounds, target = target, show = False, fig = fig, ax = ax)
    plot_size_progression("Ch", bounds = bounds, target = target, show = True, fig = fig, ax = ax)

def plot_chain_progression_multi( bounds = (3,20), target = "p"):

    #bounds[0] has to be >=3
    chain3_bounds = (bounds[0]-2, bounds[1]-3)
    chain4_bounds = ( (bounds[0]-2)//2, (bounds[1]-3)//2)

    fig, ax = plot_chain_progression( QWGraph.Ring(3), bounds = chain3_bounds, target = target, show = False)
    fig, ax = plot_chain_progression( QWGraph.Ring(4), bounds = chain4_bounds, target = target, show = False, fig = fig, ax = ax)
    fig, ax = plot_chain_progression( QWGraph.SquareCut(), bounds = chain3_bounds, target = target, show = False, fig = fig, ax = ax)
    plot_size_progression( "L", bounds = bounds, target = target, show = True, fig = fig, ax = ax)

def chain_ch_comp( bounds = (3,20), target = "p"):

    chain3_bounds = (bounds[0]-2, bounds[1]-3)

    fig, ax = plot_chain_progression( QWGraph.Ring(3), bounds = chain3_bounds, target = target, show = False)
    fig, ax = plot_size_progression( "Ch", bounds = bounds, target = target, show = False, fig = fig, ax = ax)
    plot_size_progression( "L", bounds = bounds, target = target, show = True, fig = fig, ax = ax)
