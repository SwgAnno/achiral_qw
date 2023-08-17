import sys
from article import *
from plotter import *
from graph import * 
from collection import *
from trends import *
from simulator import *

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
