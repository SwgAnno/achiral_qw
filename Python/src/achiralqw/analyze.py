from achiralqw.graph import QWGraphBuilder, QWGraph
from achiralqw.simulator import QutipSESolver, EigenSESolver, SESolver
import numpy as np
from numpy.typing import NDArray
from scipy import optimize as opt

from typing import Tuple, List

def dim(gr : QWGraph):
    """
    shorthand for the number of free phases of the analyzed graph instance
    """
    return gr.get_phase_n()

def phase_sample(step = 100):
    """
    helper function to get a discrete sample of phase values between 0 and 2*pi
    """
    return np.linspace(0, np.pi*2, step)

#DIY optimization method
#it divides the sample bounds into n smaller sections
# according to a given lenght or a set number
#and there it tries to run scipy.minimize
# results try to replicate scipy output
def section_minimize(f, bounds, f_prime = None, n_sec = None, sec_size = None):

    if n_sec != None:
        b_vec = np.linspace(bounds[0], bounds[1], n_sec+1)
    elif sec_size != None :
        b_vec = np.linspace(bounds[0], bounds[1], \
                            (bounds[1]-bounds[0])//sec_size + 2)
    else :
         b_vec = np.linspace(bounds[0], bounds[1], 11)
         #1 order of magintude finer

    #print(n_sec)
         
    sol_vec    = np.empty( len(b_vec)-1)
    f_sol_vec  = np.empty( len(b_vec)-1)

    mini_opts = dict()
    mini_opts["disp"] = False

    for i in range( len(b_vec)-1):

        #print(i, b_vec[i])
        sol = opt.minimize(f, \
                           x0 = (b_vec[i+1]+ b_vec[i])/2, \
                           bounds = [(b_vec[i],b_vec[i+1])], \
                           jac = f_prime, method="L-BFGS-B", \
                            options = mini_opts)
        
        sol_vec[i] = sol.x[0]
        f_sol_vec[i] =  f( sol.x[0])

    out = dict()
    out["x"] =sol_vec[np.argmin(f_sol_vec)]
    out["f"] = min(f_sol_vec)

    return out

class TransportParameters(dict):
    """
    This object holds the collection of relevant parameters used to analyze transport performance
    Having a single object instead of many parameters allows for easy and readable consistency of the transport analysyz
    while considering a big collection of graphs or families

    Parameter list:
        # event_s: time extent of transport events (used as a reference for brute force optimization)
        # evt_mode: Criterion of selection of best transport events 
            -TC : best maximum for a given choice of phases is the result of a optimization routine on a linearly scaling window
            -first : best maximum for a given phase choice is the first maximum encountered
        # TIME_CONSTANT: proportionality constant of the TC method time window
        # opt_mode: ALgorithm for selecting the best phase dependent transport event
            -none : no phase optimization
            -min : standard(multivariate) optimization on all the phases space
            -smart : choose the best among multiple of pi/2
            -fix : no phase optimization (but arbitraty choice of a fixed phase stored in fix_phi)
            -yolo : brute force grid evaluation op the phases space of transport performance in search of the global maximum (spoiler: it's inefficient)
        # diag: boolean, flatten the dimensionality of the graph considering only the equal phases subspace of the phase space
        # fix_phi: default parameter for the fix optimtization mode
        # solver_mode: Solver backend
            -eigen :  SESolver is a EigenSESolver
            -qutip : SESolver is QutipSESolver
    """

    _evt_modes = ["TC", "first"]
    _opt_modes = ["none","min","yolo", "smart", "fix"]
    _solver_modes = ["eigen", "qutip"]

    def __init__(self, event_s = 1, TC = 1, evt_mode = "TC", opt_mode = "none",solver_mode = "eigen", diag = True): 

        if not solver_mode in TransportParameters._solver_modes:
            raise ValueError("Solver mode not supported/recognized")
        self.solver_mode = solver_mode

        if not evt_mode in TransportParameters._evt_modes:
            raise ValueError("Mode not supported/recognized")
        else :
            self.evt_mode = evt_mode

        if not opt_mode in TransportParameters._opt_modes:
            raise ValueError("Opt mode not supported/recognized")
        else :
            self.opt_mode = opt_mode

        self.event_size = event_s
        self.TIME_CONSTANT = TC
        
        self.fix_phi = None
        self.diag = diag


    def get_solver(self) -> SESolver:
        if self.solver_mode == "qutip" :
            return QutipSESolver()
        else :
            return EigenSESolver()
        
    def max_search_time(self, gr : QWGraph) -> int:
        """
        Return the time window size of the TC criterion
        """
        return gr.distance() * self.TIME_CONSTANT

    def get_label(self, gr : QWGraph, mode = True) -> str :
        """
        Generate an automatic plot label that represent the internal state of the transprot parameters
        """

        if not mode:
            return gr.code

        mode_label = self.evt_mode

        if self.evt_mode == "TC":
            mode_label = chr(957) + "=" + str( self.TIME_CONSTANT )

        return gr.code + " " + mode_label
    
    def fix_phase(self, gr : QWGraph) :
        """
        Given the old transport parameter get the optimum phase, set it as a local variable and switch to fix mode
        """

        if self.opt_mode != "fix":
            best_phi = optimum_phase(gr, tp = self)
            self.fix_phi = best_phi
            self.opt_mode = "fix"

        return self.fix_phi

###########################################
## Trasport analysis methods

def locate_max_TC(gr : QWGraph, tp : TransportParameters = TransportParameters()) -> Tuple[float,float]:
    """
    Localize the Best Maximum according to TC criterium:
    -Search on a linarly scaling widonws defined as TC*graph.distance()
    -Devide the time window in smaller intervals and perform a standard optimization routine on each one
    -Best Maximum is the highest maximum found

    #Returns
    (float,float) -> time at which betst maximum occurs, best probability

    """
    # minimize works fine for TC near 1 (aka few maxima actually there)
    # Not reliable for higher TC
    # therefore we split the search interval into
    # smaller pieces according to a reference "event size"
        
    start = 0
    end = tp.max_search_time(gr)
    solver = tp.get_solver()

    def f(t):
        return -1*solver.evolve_default_p( gr, t)

    def f_prime(t):
        return -1*solver.evolve_default_p_deriv( gr, t)

    b_vec = np.linspace(start,end, (end-start)//tp.event_size + 2)
    sol_vec = np.empty( len(b_vec)-1)

    for i in range( len(b_vec)-1):
        
        sol = opt.minimize(f, \
                            x0 = (b_vec[i]+b_vec[i+1])/2, \
                            bounds = [(b_vec[i],b_vec[i+1])], \
                            jac = f_prime)
        
        sol_vec[i] = sol.x[0]

    probs = solver.evolve_default_p( gr, sol_vec)
    return ( sol_vec[np.argmax(probs)], max(probs))

def locate_max_first(gr : QWGraph, tp : TransportParameters = TransportParameters)  -> Tuple[float,float]:
    """
    Localize best maximum according to first maximum criterion
    Sample the evolution probability on subsequent interval,
    which increas exponentially as no sign invertion is found

    When the first sign inversion is found perform an optimization routine ( root_scalar) on that interval

    #Returns
    (float,float) -> time at which betst maximum occurs, best probability
    """
    #todo :digits are definitely note carefully researched, a check is needed

    start = .001
    end = gr.distance()/4
    exp_scale = 1.5
    evt_sample_scale = .1
    sign_change_safe = -1e-24
    maxiter = 10

    solver = tp.get_solver()

    res = None
    for i in range(maxiter) :
        #print("looking for first max in "+ str(start) +" - "+ str(end) )

        sample = np.linspace(start, end, int((end-start)//(tp.event_size* evt_sample_scale))+2)

        #print(sample)
        deriv_evo = solver.evolve_default_p_deriv( gr, sample)

        for i in range(len(sample)-1):

            if deriv_evo[i]*deriv_evo[i+1] < sign_change_safe :
##                        print("Found deriv sign inversion at  " \
##                              + str(sample[i]) +" - " \
##                              + str(sample[i+1]) )
##                        print("   -Derivative values at extrama:  " \
##                              + str(deriv_evo[i]) +" - " \
##                              + str(deriv_evo[i+1]) )
                if deriv_evo[i]*deriv_evo[i+1] >= 0:
                    pass
                else :

                    def obj(t):
                        return solver.evolve_default_p_deriv(gr, t)
                    res = opt.root_scalar( obj, \
                                            bracket = [sample[i], \
                                                        sample[i+1]])
                break
                #find solution
        if res:
            break
        else :
            start, end = sample[-1], (end + (end-start)*exp_scale)

    if not res:
        root = end
    else:
        root = res.root

    return (root, solver.evolve_default_p( gr, root)[0])

def locate_max(gr : QWGraph, tp = TransportParameters()) -> Tuple[float,float]:
    """
    Localize the "Best Maximum" according to the chosen criteria and phases

    wrapper for _locate_max_TC and _locate_max_first
    """

    if tp.evt_mode == "TC" :
        return locate_max_TC(gr=gr, tp= tp)

    elif tp.evt_mode == "first" :
        return locate_max_first(gr=gr, tp= tp)
    
def performance(gr : QWGraph, phi_vec = None, target = "p", tp = TransportParameters()) -> float:
    """
    Evaluate the transport performance for a given phases configuration

    #Returns

    float -> loc probability for best max (target == "t") of transport time of best max (target == "p")
    """

    if not np.any(phi_vec):
        phi_vec = np.repeat(0, dim(gr))

    gr.rephase(phi_vec, UPDATE_EIGEN=True)

    if target == "p":
        return locate_max(gr = gr, tp = tp)[1]
    elif target == "t":
        return locate_max(gr = gr, tp = tp)[0]

def performance_diag(gr : QWGraph, phi, target = "p", tp = TransportParameters()) -> float:
    """
    Evaluate the diag transport performance for a given diag phase

    #Returns

    float -> loc probability for best max (target == "t") of transport time of best max (target == "p")
    """

    #print(phi)

    gr.rephase( np.repeat( phi, dim(gr) ),  UPDATE_EIGEN=True)
    
    if target == "p":
        return locate_max(gr = gr, tp = tp)[1]
    elif target == "t":
        return locate_max(gr = gr, tp = tp)[0]
    
def performance_grid(gr : QWGraph, sample_step = 100, target = "p", tp : TransportParameters = TransportParameters()) -> NDArray:
    """
    Evaluate the performance function on a fine grid of points
    Each phase can take #step values in [0,2*pi) 

    #Returns

    ndarray[float] : sample_step^(# phases) array with performance evaluations
    """
    
    sample = phase_sample(sample_step)
    if target == "p":
        target = 1
    else:
        target = 0
        
    n_sample = [sample] * dim(gr)
    grid = np.meshgrid( *n_sample)
        
    out_shape = [sample_step] * dim(gr)
    out = np.empty(sample_step** dim(gr))
    out = np.reshape(out, out_shape)

    it = np.nditer(out, flags=['multi_index'])
    for val in it:

        i = it.multi_index
        
        phi_vec = np.empty(dim(gr))
        for j in range(dim(gr)):
            phi_vec[j] = grid[j][i]

        gr.rephase(phi_vec,  UPDATE_EIGEN=True)
        out[i] = locate_max(gr = gr, tp=tp)[target]

    return out

def performance_grid_diag(gr : QWGraph, sample_step = 100, target = "p", tp : TransportParameters = TransportParameters()) -> NDArray:
    """
    Evaluate the diag performance function on a fine grid of points
    The diag phase can take #step values in [0,2*pi) 

    #Returns

    ndarray[float] : sample_step array with diag performance evaluations
    """

    sample = phase_sample(sample_step)
    if target == "p":
        target = 1
    else:
        target = 0

    out = np.empty(sample_step)
    for i in range(len(sample)):
        gr.rephase( np.repeat( sample[i], dim(gr) ),  UPDATE_EIGEN=True)

        out[i] = locate_max(gr = gr, tp = tp)[target]

    return out
    
def optimum_phase_yolo(gr : QWGraph, step = 100, tp :TransportParameters= TransportParameters()) -> Tuple[float,float] | Tuple[List[float], float]:
    """
    Compute the best transport phases according o the chosen time criterion.
    Transport performance is eavluated on a fine grid where each phase can take
    #step value in [0,2*pi)
    Return the best configuration of the grid
    WARNING: scales as step^(# of phases)

    #Returns

    (list[float], float ) -> vector of best phases, best performance

    or

    (float, float) -> best diag phase, relative best diag performance
    """

    if tp.diag:
        perf = performance_grid_diag(gr,step, tp = tp)
    else :
        perf = performance_grid(gr,step, tp = tp)

    pos_max = np.argmax(perf)

    return pos_max/step*2*np.pi, perf[pos_max]

def optimum_phase_minimize(gr : QWGraph, tp :TransportParameters= TransportParameters()) -> Tuple[float,float] | Tuple[List[float], float]:
    """
    Compute the best transport phases according o the chosen time criterion.
    Best phase configuration os obtained as a result of a n-variables optimization routine
    on the performance function

    Pherformance function with TC criterion has often discontinuities on its derivative,
    results may not be factual

    #Returns

    (list[float], float ) -> vector of best phases, best performance

    or

    (float, float) -> best diag phase, relative best diag performance
    """

    if dim(gr) == 0 :
        return (0, performance(gr, 0, tp = tp))

    #minimize the inverse of perf function
    if tp.diag:
        def perf(x):
            return -1*performance_diag(gr,x, tp=tp)

        #todo: add a smart way to figure out the number of sections

        sol = section_minimize(perf, bounds = [0, 2*np.pi])

    else:
        def perf(x):
            return -1* performance(gr,x, tp = tp)

        sol = opt.minimize(perf, \
                    x0 = np.repeat(.5, dim(gr))   , \
                    bounds = [(0, 2*np.pi)]* dim(gr) )

    return sol["x"], -1*perf(sol["x"])

#todo: does it work???
def optimum_phase_smart(gr : QWGraph, tp : TransportParameters = TransportParameters())  -> Tuple[float,float] | Tuple[List[float], float]:
    """
    Compute the best transport phases according o the chosen time criterion.
    Transport performance is eavluated on a coarse grid where each phase is chosen as a
    multiple of pi/2 (0, pi/2, pi, -pi/2).
    Return the best configuration of the grid
    WARNING: scales as step^(# of phases)

    #Returns

    (list[float], float ) -> vector of best phases, best performance

    or

    (float, float) -> best diag phase, relative best diag performance
    """

    if dim(gr) == 0 :
        return 0 , performance(gr, tp = tp)

    if tp.diag:
        sample = phase_sample(step = 5)

        out = 0
        best = 0

        for phase in sample:

            diag_phase = np.repeat(phase, dim(gr)) 
            gr.rephase(diag_phase,  UPDATE_EIGEN=True)
            cur = locate_max(gr = gr, tp = tp)[1]

            if cur > best:
                best = cur
                out = phase

        return out, performance_diag(gr, out, tp = tp)

    sample = phase_sample(step = 5)
    n_sample = [sample] * dim(gr)
    grid = np.meshgrid( *n_sample)

    best = 0
    out = np.repeat(0, dim(gr))

    it = np.nditer(grid[0], flags=['multi_index'])
    for val in it:

        i = it.multi_index
        
        phi_vec = np.empty(dim(gr))
        for j in range(dim(gr)):
            phi_vec[j] = grid[j][i]

        gr.rephase(phi_vec,  UPDATE_EIGEN=True)
        cur = locate_max(gr = gr, tp = tp)[1]

        if cur > best:
            best = cur
            out = phi_vec

    return out, performance(gr, out, tp = tp)
    
def optimum_phase(gr : QWGraph, tp : TransportParameters = TransportParameters()) -> float | List[float]:
    """
    Localize the "Best transport phase" according to the chosen transport criteria

    wrapper for all the optimum phases methods
    """
    best_phi = 0

    if tp.opt_mode == "none":
        best_phi = 0
    elif tp.opt_mode == "smart" :
        best_phi = optimum_phase_smart(gr= gr, tp = tp)[0]
    elif tp.opt_mode == "fix" :
        assert tp.fix_phi is not None
        best_phi = tp.fix_phi
    elif tp.opt_mode == "yolo" :
        best_phi = optimum_phase_yolo(gr= gr, tp = tp)[0]
    elif tp.opt_mode == "min" :
        best_phi = optimum_phase_minimize(gr= gr, tp = tp)[0]

    return best_phi

    
def performance_best(gr : QWGraph, target = "p", tp : TransportParameters = TransportParameters()):
    """
    Compute the best transport performance according to current time and phase criteria

    #Returns

    float -> loc probability for best max (target == "t") of transport time of best max (target == "p")
    """
    best_phi = 0

    if tp.opt_mode == "none":
        best_phi = 0
    elif tp.opt_mode == "smart" :
        best_phi = optimum_phase_smart(gr, tp = tp)[0]
    elif tp.opt_mode == "fix" :
        assert tp.fix_phi is not None
        best_phi = tp.fix_phi
    elif tp.opt_mode == "yolo" :
        best_phi = optimum_phase_yolo(gr, tp = tp)[0]
    elif tp.opt_mode == "min" :
        best_phi = optimum_phase_minimize(gr, tp = tp)[0]

    #print(best_phi, tp.opt_mode)

    if tp.diag:
        return performance_diag(gr, phi=best_phi, target = target, tp = tp)
    else:
        return performance(gr, phi_vec=best_phi, target = target, tp = tp)
    
#################################################