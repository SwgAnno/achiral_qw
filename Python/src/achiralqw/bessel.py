import numpy as np
from scipy.special import *

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