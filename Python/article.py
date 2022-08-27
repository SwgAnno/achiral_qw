#list of methonds to draw specific plots in the article

from trends import *

def plot_simpler_topo_progression( bounds = (3,12), target = "p", x_mode = "dist",mode = "TC", TC = 1):

    gr_list = []
    for i in range( bounds[1]- bounds[0] +1):
        gr_list.append( QWGraph.Line( bounds[0]+ i))
    
    out = optimized_progression(gr_list, target = target, mode = mode, TC = TC)
    x = get_list_x(gr_list, x_mode = x_mode)

    fig, ax = plot_standard_progression([x, out[1]], target = target, x_mode = x_mode, show = False)

    gr_list = []
    for i in range( bounds[1]- bounds[0] +1):
        gr_list.append( QWGraph.Ring( bounds[0]+ i*2))
    
    out = optimized_progression(gr_list, target = target, mode = mode, TC = TC)
    x = get_list_x(gr_list, x_mode = x_mode)

    fig, ax = plot_standard_progression([x, out[1]], target = target, x_mode = x_mode, show = False, fig = fig, ax = ax)

    gr_list = []
    for i in range( bounds[1]- bounds[0] +1):
        gr_list.append( QWGraph.Ring( bounds[0]+ i*2, HANDLES = True))
    
    out = optimized_progression(gr_list, target = target, mode = mode, TC = TC)
    x = get_list_x(gr_list, x_mode = x_mode)

    plot_standard_progression([x, out[1]], target = target, x_mode = x_mode, show = True, fig = fig, ax = ax)

def plot_odd_even_progression( bounds = (3,12), target = "p", x_mode = "dist",mode = "TC", TC = 1) :

    gr_list = []
    for i in range( bounds[1]- bounds[0] +1):
        gr_list.append( QWGraph.Line( bounds[0]+ i))
    
    out = optimized_progression(gr_list, target = target, mode = mode, TC = TC)
    x = get_list_x(gr_list, x_mode = x_mode)

    fig, ax = plot_standard_progression([x, out[1]], target = target, x_mode = x_mode, show = False)

    gr_list = []
    for i in range( bounds[1]- bounds[0] +1):
        gr_list.append( QWGraph.Ring( bounds[0]+ i*2))
    
    out = optimized_progression(gr_list, target = target, mode = mode, TC = TC)
    x = get_list_x(gr_list, x_mode = x_mode)

    fig, ax = plot_standard_progression([x, out[1]], target = target, x_mode = x_mode, show = False, fig = fig, ax = ax)

    gr_list = []
    for i in range( bounds[1]- bounds[0] +1):
        gr_list.append( QWGraph.Ring( bounds[0]+1+ i*2))
    
    out = optimized_progression(gr_list, target = target, mode = mode, TC = TC)
    x = get_list_x(gr_list, x_mode = x_mode)

    plot_standard_progression([x, out[1]], target = target, x_mode = x_mode, show = True, fig = fig, ax = ax)

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
    print( get_line_data( bounds = (4,20), target = "t"))

    print( "even/line m : \t" , coeff[0][0]/coeff[2][0])
    print( "odd/line m : \t" , coeff[1][0]/coeff[2][0])


