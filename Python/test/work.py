import sys
from achiralqw.article import *
from achiralqw.plotter import *
from achiralqw.graph import * 
from achiralqw.collection import *
from achiralqw.trends import *
from achiralqw.simulator import *

an = Analyzer(mode = "first")
plot_chain_progression_multi_loglog(bounds = (50,100), points = 20, target = "p", analyzer = an, fast = True)

if len(sys.argv) > 1:
    plt.savefig(sys.argv[1])
else :
    file = input("Indicare il nome del file per salvare i risultati (lasciare vuoto per visualizzare)")   

    if file == "":
        plt.show() 
    else :
        plt.savefig(file)
