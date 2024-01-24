import numpy as np
from scipy.special import jv

def line_evo_bessel(n, l = None, x = .0): 

    if l == None:
        l = n
    diff = 2*n-2*l
    x = 2*x

    #always do at least an iteration
    out = np.power(-1j,l)*jv(l,x) + np.power(-1j,l+diff)*jv(l+diff,x)
    l = l+2*n

    while l < 2*x:
        out = out + np.power(-1j,l)*jv(l,x) + np.power(-1j,l+diff)*jv(l+diff,x)
        l = l+2*n

    #todo: is this legal?
    if l%n != 0 :
        out = out * np.sqrt(2)

    return np.power(abs(out), 2)

def ring_evo_bessel(n, l = None, x = .0): 

    if l == None:
        l = int(n/2)
    diff = n-2*l
    x = 2*x

    #always do at least an iteration
    out = np.power(-1j,l)*jv(l,x) + np.power(-1j,l+diff)*jv(l+diff,x)
    l = l+n

    while l < 2*x:
        out = out + np.power(-1j,l)*jv(l,x) + np.power(-1j,l+diff)*jv(l+diff,x)
        l = l+n

    return np.power(abs(out), 2)

#simulate first maximum transport performance with besse functions
#def ring_pefrormance_bessel(l,  )


#########################################
## plotting helpers
def check_line_bessel_dist(l = 5, end = 10):

    grid = np.arange(0,end, 1)
    eval = np.zeros( len(grid))

    for t in range(len(grid)):
        for i in range(l):
            eval[t] = eval[t] + line_evo_bessel(l-1,i, grid[t])

    print(eval)

def check_ring_bessel_dist(l = 5, end = 10):

    grid = np.arange(0,end, 1)
    eval = np.zeros( len(grid))

    for t in range(len(grid)):
        for i in range(l):
            eval[t] = eval[t] + ring_evo_bessel(l,i, grid[t])

    print(eval)

if __name__ == "__main__" :
    #plot_line_vs_bessel(l = 5, trace_conn = False)
    #bessel_progression(bounds = (2,50), target = "p", L_ref = True)

    #plot_line_bessel_evolution(5,end = 10)
    check_line_bessel_dist(5)

    #plt.show()