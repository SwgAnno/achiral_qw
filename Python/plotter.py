from simulator import *
import matplotlib.pyplot as plt

TC = 1
qut = False

def plot_evo_mat(gr , start = 0, end = None, by = .1):

    if not end:
        global TC
        end = gr.N * TC
    global qut

    seq = np.arange(start,end,by)
    
    solver = SESolver(gr, qutip = qut)
    evo = solver.evo_p_psi( gr.get_start_state(qut), seq)

    print("Grid-wise maximum: ", max(evo[gr.target][:]))

    fig, ax = plt.subplots()
    for i in range(gr.N) :
        ax.plot(seq, evo[i][:])
    ax.set_xlabel('Time')
    ax.set_ylabel('Expectation values')
    
    ax.legend(list(range(0,gr.N)))
    plt.show()

def plot_evo_vs_derivative(gr, l = 0, start = 0, end = None, by = .1):

    if not end:
        global TC
        end = gr.N * TC
    global qut
    
    seq = np.arange(start,end,by)
    
    solver = SESolver(gr,qutip = qut)
    evo = solver.target_p(seq)
    deriv = solver.target_p_prime(seq)

    fig, ax = plt.subplots()

    ax.plot(seq, evo)
    ax.plot(seq, deriv)
    
    ax.set_xlabel('Time')
    ax.set_ylabel('P')

    ax.legend(["P", "dP/dt"])
    plt.show()

def plot_evo_vs_qutip(gr, l = 0, start = 0, end = None, by = .1):
    if not end:
        global TC
        end = gr.N * TC

    seq = np.arange(start,end,by)

    solver = SESolver(gr)
    evo = solver.target_p(seq)

    solver.qut = True
    evo_q = solver.target_p(seq)

    print( np.abs(evo-evo_q) < .0001)

    fig, ax = plt.subplots()

    ax.plot(seq, evo)
    ax.plot(seq, evo_q)
    
    ax.set_xlabel('Time')
    ax.set_ylabel('P')

    ax.legend(["P", "P qutip"])
    plt.show()

def plot_performance(gr, sample_step = 100, diag = False, mode = "TC"):

    global TC
    global qut
    an = Analyzer(gr, TC = TC, qutip = qut, mode = mode)

    if diag:
        plot_performance_diag(sample_step, an)
    elif gr.get_phase_n() == 1 :
        plot_performance_1(sample_step, an)
    elif gr.get_phase_n() == 2 :
        plot_performance_2(sample_step, an)
    else :
        print("")
        

def plot_performance_diag(sample_step, an):
    seq = np.linspace(0, np.pi*2, sample_step)

    perf = []
    perf = an.performance_diag(sample_step)

    fig, ax = plt.subplots()
    
    ax.plot(seq, perf)
    
    ax.set_xlabel('phi')
    ax.set_ylabel('max P')

    plt.show()    

def plot_performance_1(sample_step, an):
    seq = np.linspace(0, np.pi*2, sample_step)

    perf = an.performance(sample_step)

    fig, ax = plt.subplots()
    
    ax.plot(seq, perf)
    
    ax.set_xlabel('phi')
    ax.set_ylabel('max P')

    plt.show()

def plot_performance_2(sample_step, an):
    seq = np.linspace(0, np.pi*2, sample_step)

    perf = an.performance(sample_step)

    fig, ax = plt.subplots()
    
    ax.pcolormesh(seq, seq, perf)
    
    ax.set_xlabel('phi 1')
    ax.set_ylabel('phi 2')

    plt.show()

    pass
    
################

#a  = QWGraph.chain(QWGraph.Ring(3), 10)
a = QWGraph.Ring(6)
a = a+a
TC = 1
qut = False

test = Analyzer(a, qutip = False)

test.mode = "first"
print(test.locate_max())
test.mode = "TC"
print(test.locate_max())


#plot_evo_vs_derivative(a)
#plot_evo_vs_qutip(a)
#plot_evo_mat(a)
plot_performance(a,100, diag = False, mode = "TC")
#plot_evo_vs_derivative(a|a+a)

