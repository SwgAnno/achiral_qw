from simulator import *
import matplotlib.pyplot as plt

def plot_evo_mat(gr , start = 0, end = 5, by = .1):

    seq = np.arange(start,end,by)
    
    solver = SESolver(gr)
    evo = solver.evo_p_l_t( gr.get_start_state(), seq)

    fig, ax = plt.subplots()
    for i in range(gr.N) :
        ax.plot(seq, evo[i,:])
    ax.set_xlabel('Time')
    ax.set_ylabel('Expectation values')

    ax.legend(list(range(0,gr.N)))
    plt.show()





################
