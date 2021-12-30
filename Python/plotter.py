from simulator import *
import matplotlib.pyplot as plt

TC = 1

def plot_evo_mat(gr , start = 0, end = None, by = .1):

    if not end:
        global TC
        end = gr.N * TC

    seq = np.arange(start,end,by)
    
    solver = SESolver(gr)
    evo = solver.evo_p_l( gr.get_start_state(), seq)

    print("Grid-wise maximum: ", max(evo[gr.target,:]))

    fig, ax = plt.subplots()
    for i in range(gr.N) :
        ax.plot(seq, evo[i,:])
    ax.set_xlabel('Time')
    ax.set_ylabel('Expectation values')
    
    ax.legend(list(range(0,gr.N)))
    plt.show()

def plot_evo_vs_derivative(gr, l = 0, start = 0, end = None, by = .1):

    if not end:
        global TC
        end = gr.N * TC
    
    seq = np.arange(start,end,by)
    
    solver = SESolver(gr)
    evo = solver.evo_p_l( gr.get_start_state(), seq)
    deriv = solver.deriv_p_l( gr.get_start_state(), seq)

    fig, ax = plt.subplots()

    ax.plot(seq, evo[l,:])
    ax.plot(seq, deriv[l,:])
    
    ax.set_xlabel('Time')
    ax.set_ylabel('P')

    ax.legend(["P", "dP/dt"])
    plt.show()

def plot_performance(gr, sample_step = 100):

    global TC
    an = Analyzer(gr, TC = TC)

    seq = np.linspace(0, np.pi*2, sample_step)

    perf = an.performance(sample_step)

    fig, ax = plt.subplots()
    
    ax.plot(seq, perf)
    
    ax.set_xlabel('phi')
    ax.set_ylabel('max P')

    plt.show()

################

a = QWGraph.Ring(7)
TC = 5


plot_performance(a)
#plot_evo_mat(a)
#plot_evo_vs_derivative(a|a+a)

