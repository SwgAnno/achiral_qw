from plotter import *
from simulator import *
from Graph import *

#a  = QWGraph.chain(QWGraph.Ring(3), 10)
#a = QWGraph.Parallel(3,2)
#a = QWGraph.chain(a,2)
##a = a+a
a = QWGraph.Ring(6)
TC = 1
qut = True

plot_performance(a)

plot_evo_vs_derivative(a)

test = Analyzer(a, qutip = True)

test.mode = "first"
print(test.locate_max())
test.mode = "TC"
print(test.locate_max())


#plot_evo_vs_derivative(a)
#plot_evo_vs_qutip(a)
#plot_evo_mat(a)
##plot_performance(a,100, mode = "diag", an_mode = "TC")
##a.plot()
#plot_evo_vs_derivative(a|a+a)

#print(test.optimum_phase_minimize(diag = True))
