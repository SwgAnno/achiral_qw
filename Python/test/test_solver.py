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

if __name__ == "__main__":
    
    qwgb = QWGraphBuilder()

    test1 = qwgb.Ring(5, COMPUTE_EIGEN= True)

    fig, axx = plt.subplots(3,1, figsize = (5,10))

    my_solver = EigenSESolver()
    plot_evo_vs_derivative( test1, TC = 4, solver = my_solver, ax = axx[0])
    plot_evo_mat(test1, TC = 4, solver = my_solver, ax = axx[1])
    plot_evo_mat_heatmap(test1, TC = 4, solver = my_solver,fig = fig, ax = axx[2])

    #my_solver = QutipSESolver()
    #plot_evo_vs_derivative( test1, TC = 4, solver = my_solver, ax = axx[1])

    plt.show()