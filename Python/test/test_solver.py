import matplotlib.pyplot as plt


from achiralqw.graph import QWGraph, QWGraphBuilder
from achiralqw.simulator import EigenSESolver,QutipSESolver
from achiralqw.plotter import *


def test_sesolver_eigen():

    gb = QWGraphBuilder()
    t_vec = [1.0, 2.0, 3.0, 4.0, 5.0]

    test = gb.Ring(6, COMPUTE_EIGEN= True)

    solver = EigenSESolver()

    res = solver.evolve_default_p(test, t_vec)
    exp = [0.06650288, 0.7369635,  0.03505077, 0.69609166, 0.20971131]
    np.testing.assert_allclose(res, exp, rtol = 1e-5)

    res = solver.evolve_default_p_deriv(test, t_vec)
    exp = [ 0.3288674,   0.27184382, -0.48680875,  0.56527321, -0.68552995]
    np.testing.assert_allclose(res, exp, rtol = 1e-5)

    test.rephase(np.pi/2)
    #new eigen are not computed automatically
    test.update_eigen()
    res = solver.evolve_default_p(test, t_vec)
    exp = [0.03325401, 0.37171176, 0.03695258, 0.05451863, 0.23905042]
    np.testing.assert_allclose(res, exp, rtol = 1e-5)

    res = solver.evolve_default_p_deriv(test, t_vec)
    exp = [ 0.16446268,  0.14973757, -0.29957085, -0.10177815,  0.79405903]
    np.testing.assert_allclose(res, exp, rtol = 1e-5)

    test = gb.Ring(5, COMPUTE_EIGEN= True)

    res = solver.evolve_default_p(test, t_vec)
    exp = [0.14114771, 0.32800068, 0.14961465, 0.12044615, 0.07631292]
    print(res)
    np.testing.assert_allclose(res, exp, rtol = 1e-5)

    res = solver.evolve_default_p_deriv(test, t_vec)
    exp = [ 0.39848401, -0.25359153,  0.09991786, -0.20888723,  0.05222153]
    print(res)
    np.testing.assert_allclose(res, exp, rtol = 1e-5)

    test.rephase(np.pi/2)
    #new eigen are not computed automatically
    test.update_eigen()
    res = solver.evolve_default_p(test, t_vec)
    exp = [0.23196217, 0.61333875, 0.04019297, 0.23850944, 0.25618063]
    print(res)
    np.testing.assert_allclose(res, exp, rtol = 1e-5)

    res = solver.evolve_default_p_deriv(test, t_vec)
    exp = [ 0.73747511, -0.69543845,  0.42919989, -0.61452427,  0.87168303]
    print(res)
    np.testing.assert_allclose(res, exp, rtol = 1e-5)

def test_sesolver_qutip():

    gb = QWGraphBuilder()
    t_vec = [1,2,3,4,5]

    test = gb.Ring(6, COMPUTE_EIGEN= True)

    solver = QutipSESolver()

    res = solver.evolve_default_p(test, t_vec)
    exp = [0.06650288, 0.7369635,  0.03505077, 0.69609166, 0.20971131]
    np.testing.assert_allclose(res, exp, rtol = 1e-5)

    res = solver.evolve_default_p_deriv(test, t_vec)
    exp = [ 0.3288674,   0.27184382, -0.48680875,  0.56527321, -0.68552995]
    np.testing.assert_allclose(res, exp, rtol = 1e-5)

    t_vec = [0,1,2,3,4,5]

    res = solver.evolve_default_p(test, t_vec)
    exp = [0, 0.06650288, 0.7369635,  0.03505077, 0.69609166, 0.20971131]
    np.testing.assert_allclose(res, exp, rtol = 1e-5)

    res = solver.evolve_default_p_deriv(test, t_vec)
    exp = [0, 0.3288674,   0.27184382, -0.48680875,  0.56527321, -0.68552995]
    np.testing.assert_allclose(res, exp, rtol = 1e-5)

    t_vec = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]

    res = solver.evolve_default_p(test, t_vec)
    exp = [0, 0.06650288, 0.7369635,  0.03505077, 0.69609166, 0.20971131]
    np.testing.assert_allclose(res, exp, rtol = 1e-5)

    res = solver.evolve_default_p_deriv(test, t_vec)
    exp = [0, 0.3288674,   0.27184382, -0.48680875,  0.56527321, -0.68552995]
    np.testing.assert_allclose(res, exp, rtol = 1e-5)

    t_vec = [1,2,3,4,5]
    test.rephase(np.pi/2)
    res = solver.evolve_default_p(test, t_vec)
    exp = [0.03325401, 0.37171176, 0.03695258, 0.05451863, 0.23905042]
    np.testing.assert_allclose(res, exp, rtol = 1e-5)

    res = solver.evolve_default_p_deriv(test, t_vec)
    exp = [ 0.16446268,  0.14973757, -0.29957085, -0.10177815,  0.79405903]
    np.testing.assert_allclose(res, exp, rtol = 1e-5)

    test = gb.Ring(5, COMPUTE_EIGEN= True)

    res = solver.evolve_default_p(test, t_vec)
    exp = [0.14114771, 0.32800068, 0.14961465, 0.12044615, 0.07631292]
    print(res)
    np.testing.assert_allclose(res, exp, rtol = 1e-5)

    res = solver.evolve_default_p_deriv(test, t_vec)
    exp = [ 0.39848401, -0.25359153,  0.09991786, -0.20888723,  0.05222153]
    print(res)
    np.testing.assert_allclose(res, exp, rtol = 1e-5)

    test.rephase(np.pi/2)
    res = solver.evolve_default_p(test, t_vec)
    exp = [0.23196217, 0.61333875, 0.04019297, 0.23850944, 0.25618063]
    print(res)
    np.testing.assert_allclose(res, exp, rtol = 1e-5)

    res = solver.evolve_default_p_deriv(test, t_vec)
    exp = [ 0.73747511, -0.69543845,  0.42919989, -0.61452427,  0.87168303]
    print(res)
    np.testing.assert_allclose(res, exp, rtol = 1e-5)

def test_specific_solver():
    #compare the non default methods for evolution
    
    gr = QWGraphBuilder.Ring(5)

    esolver = EigenSESolver()
    qsolver = QutipSESolver()

    t = 5

    #state vector
    exp = [-0.66088495-0.47713899j, \
           -0.21553422+0.21217166j, \
           -0.21553422+0.21217166j, \
            0.12644093-0.24561273j, \
            0.12644093-0.24561273j ]
    eres = esolver.evolve_state(gr, psi = gr.get_start_state(), t = t)
    eres = [eres[i][0] for i in range(len(exp))]
    np.testing.assert_allclose(eres, exp, rtol = 1e-5)
    qres = qsolver.evolve_state(gr, psi = gr.get_start_state(), t = t)
    qres = [qres[i][0] for i in range(len(exp))]
    np.testing.assert_allclose(qres, exp, rtol = 1e-5)

    #state vector derivative
    exp = [-0.42434333-0.43106844j, \
            0.72275171-0.53444402j, \
            0.72275171-0.53444402j, \
            0.03344106-0.08909329j, \
            0.03344106-0.08909329j ]
    eres = esolver.evolve_state_deriv(gr, psi = gr.get_start_state(), t = t)
    eres = [eres[i][0] for i in range(len(exp))]
    np.testing.assert_allclose(eres, exp, rtol = 1e-5)
    qres = qsolver.evolve_state_deriv(gr, psi = gr.get_start_state(), t = t)
    qres = [qres[i][0] for i in range(len(exp))]
    np.testing.assert_allclose(qres, exp, rtol = 1e-5)

    # site prob
    exp = [ 0.66443053, \
            0.09147181, \
            0.09147181, \
            0.07631292, \
            0.07631292]
    eres = esolver.evolve_state_p(gr, psi = gr.get_start_state(), t = t)
    eres = [eres[i][0] for i in range(len(exp))]
    np.testing.assert_allclose(eres, exp, rtol = 1e-5)
    qres = qsolver.evolve_state_p(gr, psi = gr.get_start_state(), t = t)
    qres = [qres[i][0] for i in range(len(exp))]
    np.testing.assert_allclose(qres, exp, rtol = 1e-5)

    # site prob derivative
    exp = [ 0.97224335, \
           -0.53834321, \
           -0.53834321, \
            0.05222153, \
            0.05222153]
    eres = esolver.evolve_state_p_deriv(gr, psi = gr.get_start_state(), t = t)
    eres = [eres[i][0] for i in range(len(exp))]
    np.testing.assert_allclose(eres, exp, rtol = 1e-5)
    qres = qsolver.evolve_state_p_deriv(gr, psi = gr.get_start_state(), t = t)
    qres = [qres[i][0] for i in range(len(exp))]
    np.testing.assert_allclose(qres, exp, rtol = 1e-5)


#plotting helper
#comparison between eigenvalue method and Qutip solver on computed evolution
def check_evo_vs_qutip( test_gr = QWGraphBuilder.Ring(6, COMPUTE_EIGEN= True), l = 0, start = 0, end = None, by = .1, deriv = True):
    if not end:
        end = test_gr.N * 5

    seq = np.arange(start,end,by)
    if not deriv: 
        solver = EigenSESolver()
        evo = solver.evolve_default_p( test_gr, seq)

        solver = QutipSESolver()
        evo_q = solver.evolve_default_p( test_gr, seq)

    else :
        solver = EigenSESolver()
        evo = solver.evolve_default_p_deriv( test_gr, seq)

        solver = QutipSESolver()
        evo_q = solver.evolve_default_p_deriv( test_gr, seq)

    print( np.abs(evo-evo_q) < .0001)

    fig, ax = plt.subplots()

    ax.plot(seq, evo)
    ax.plot(seq, evo_q)
    
    ax.set_xlabel('Time')
    ax.set_ylabel('P')

    ax.legend(["P", "P qutip"])
    plt.show()

#plotting helper
#see how qutip (or the old method) works with random poits evaluation
def check_evo_vs_qutip_scatter( test_gr = QWGraphBuilder.Ring(6, COMPUTE_EIGEN=True), l = 0, start = 0, end = None, by = .1, deriv = False):
    if not end:
        global TC
        end = test_gr.N * TC

    seq = np.arange(start,end,by)
    seq_q = np.linspace(start,end, 20)
    test = np.random.rand(40)
    test = test * (end-start) + start

    test_qut = test[0:20]
    test_old = test[20:40]
    
    if not deriv: 
        solver = EigenSESolver()
        evo = solver.evolve_default_p(test_gr, seq)
        evo_scat = solver.evolve_default_p(test_gr, test_old)

        solver = QutipSESolver()
        evo_q = solver.evolve_default_p(test_gr, test_qut)
        evo_seq = solver.evolve_default_p(test_gr, seq_q)

    else :
        solver = EigenSESolver()
        evo = solver.evolve_default_p_deriv(test_gr, seq)
        evo_scat = solver.evolve_default_p_deriv(test_gr, test_old)

        solver = QutipSESolver()
        evo_q = solver.evolve_default_p_deriv(test_gr, test_old)
        evo_seq = solver.evolve_default_p_deriv(test_gr, seq_q)

    fig, ax = plt.subplots()

    ax.plot(seq, evo)
    ax.scatter(test_old, evo_scat, color = "green")
    ax.scatter(test_qut, evo_q, color = "red")
    ax.scatter(seq_q, evo_seq, color = "yellow")
    
    ax.set_xlabel('Time')
    ax.set_ylabel('P')

    ax.legend(["P", "Scatter evaluation noQ","Scatter evaluation Q","Scat ordered evaluation Q"])
    plt.show()

#dumb helper
#single qutip evaluation check
def check_qutip( test_gr = QWGraphBuilder.Ring(6, COMPUTE_EIGEN=True), t_0 = [2]):

    solver = QutipSESolver()

    for _ in range(2):
        print( solver.evolve_default_p(test_gr, t_0))
        print( solver.evolve_default_p_deriv(test_gr, t_0))

if __name__ == "__main__":
    
    qwgb = QWGraphBuilder()

    test1 = qwgb.Ring(5, COMPUTE_EIGEN= True)

    fig, axx = plt.subplots(3,1, figsize = (5,10))

    my_solver = EigenSESolver()
    plot_evo_vs_derivative( test1, TC = 4, solver = my_solver, ax = axx[0])
    plot_evo_mat(test1, TC = 4, solver = my_solver, ax = axx[1])
    plot_evo_mat_heatmap(test1, TC = 4, solver = my_solver,fig = fig, ax = axx[2])

    check_evo_vs_qutip(deriv = True)
    check_evo_vs_qutip(deriv = False)
    check_qutip()
    check_evo_vs_qutip_scatter()
    plt.show()


    #my_solver = QutipSESolver()
    #plot_evo_vs_derivative( test1, TC = 4, solver = my_solver, ax = axx[1])

    plt.show()