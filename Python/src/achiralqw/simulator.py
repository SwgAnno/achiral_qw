from abc import abstractmethod
from achiralqw.graph import QWGraphBuilder, QWGraph
import numpy as np
from scipy import optimize as opt
import qutip as qt

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

    
    

#qutip has problems with non list input
#and with evaulation times not starting with 0
def format_qutip_time( t):

    if type(t) == float:
        t = [t]

    if t[0] != 0:
        if isinstance(t, (np.ndarray, np.generic)):
            t = np.insert(t,0,0.)
        else:
            t.insert(0,0.)
        
        return t,True
    else:
        return t,False


class SESolver(object):
    """
    Absrtact class that define the interface for a Evolution solver for QWGraphs
    """

    def __init__():
        pass

    @abstractmethod
    def evolve_state( self, gr : QWGraph, psi, t):
        """
        Solve the Schrodinger equation for a given pure state in the site basis
        Return the projection of the resulting state on each state
        """
        pass

    #todo check correctness
    @abstractmethod
    def evolve_state_deriv(self, gr : QWGraph, psi,t):
        """
        Compute the derivative of an evolving state at time t (Schrodinger equation)
        """
        pass

    @abstractmethod
    def evolve_state_p(self, gr, psi, t) :
        """
        Solve the Schrodinger equation for a given pure state in the site basis
        Return the squared modulus (on site probability) of the resulting state on each state
        """
        pass

    @abstractmethod
    def evolve_state_p_deriv(self, gr, psi, t):
        """
        Compute the derivative of the localization probability on each site at time t ( evolving psi according to Schrodinger equation)
        """
        pass

    @abstractmethod
    def evolve_default(self, gr : QWGraph, t):
        """
        wrapper for evolve state :
        Solves Schrodinger equation for a state localized on the graph start site
        Only returns resulting projection on the target site
        """
        pass

    @abstractmethod
    def evolve_default_deriv(self, gr : QWGraph, t):
        """
        wrapper for evolve state_deriv :
        Compute derivative according to Schrodinger equation for a state localized on the graph start site
        Only returns resulting projection on the target site
        """
        pass
    @abstractmethod
    def evolve_default_p(self, gr : QWGraph, t):
        """
        wrapper for evolve state_p :
        Solves Schrodinger equation for a state localized on the graph start site
        Only returns localization probability of the target site
        """
        pass

    @abstractmethod
    def evolve_default_p_deriv(self, gr : QWGraph, t):
        """
        wrapper for evolve state__p_deriv :
        Compute derivative of localization probability according to Schrodinger equation
        Start state is localized on the grahp start site
        Only returns resulting projection on the target site
        """
        pass


    

class EigenSESolver(SESolver) :
    """
    Schrodinger equation solver for QWGraph that uses eigenvector basis to decompose and evolve the starting state
    """

    def __init__(self):
        pass

    @staticmethod
    def _to_eiegn_basis( eig_vec , local_A):
        """
        Represent the state in the eigenvector basis
        """

        return  eig_vec.conj().T.dot(local_A)

    @staticmethod
    def _to_state_basis(eig_vec, eigen_A):
        """
        Go from eigenvector basis to state (vertices) basis
        """

        return eig_vec.dot(eigen_A)
    
    @staticmethod
    def _exp_map(eig_val, t):
        """
        return evolution map in the eigenvector basis
        """
        return np.exp(-1j * np.outer( eig_val, t) )

    def evolve_state( self, gr : QWGraph, psi, t):
        """
        Solve the Schrodinger equation for a given pure state in the site basis
        Return the projection of the resulting state on each state
        """

        #carry on calculation in  H eigenvector basis
        A_t = EigenSESolver._to_eiegn_basis(gr.eig_vec, psi) * EigenSESolver._exp_map(gr.eig_val, t)

        return EigenSESolver._to_state_basis(gr.eig_vec, A_t)

    #todo check correctness
    def evolve_state_deriv(self, gr : QWGraph, psi,t):
        """
        Compute the derivative of an evolving state at time t (Schrodinger equation)
        """
        psi_t = self.evolve_state(gr, psi,t)
        return -1j * gr.mat.dot(psi_t)

    def evolve_state_p(self, gr, psi, t) :
        """
        Solve the Schrodinger equation for a given pure state in the site basis
        Return the squared modulus (on site probability) of the resulting state on each state
        """
        return np.power(np.abs(self.evolve_state(gr, psi, t)), 2)

    def evolve_state_p_deriv(self, gr, psi, t):
        """
        Compute the derivative of the localization probability on each site at time t ( evolving psi according to Schrodinger equation)
        """
           
        A_t = EigenSESolver._to_eiegn_basis(gr.eig_vec, psi) * EigenSESolver._exp_map(gr.eig_val, t)

        A_t_prime = gr.eig_val* A_t

        #get probability derivative as 2Re[(a)(a*)']
        out = EigenSESolver._to_state_basis(gr.eig_vec, A_t) * np.conjugate(EigenSESolver._to_state_basis( gr.eig_vec, A_t_prime))

        return -2* np.imag(out)

    def evolve_default(self, gr : QWGraph, t):
        """
        wrapper for evolve state :
        Solves Schrodinger equation for a state localized on the graph start site
        Only returns resulting projection on the target site
        """
        return   self.evolve_state(gr, gr.get_start_state(), t)[gr.target, :]
    
    def evolve_default_deriv(self, gr : QWGraph, t):
        """
        wrapper for evolve state_deriv :
        Compute derivative according to Schrodinger equation for a state localized on the graph start site
        Only returns resulting projection on the target site
        """
        return self.evolve_state_deriv(gr, gr.get_start_state(), t)[gr.target, :]

    def evolve_default_p(self, gr : QWGraph, t):
        """
        wrapper for evolve state_p :
        Solves Schrodinger equation for a state localized on the graph start site
        Only returns localization probability of the target site
        """
        return self.evolve_state_p(gr, gr.get_start_state(), t)[gr.target, :]

    def evolve_default_p_deriv(self, gr : QWGraph, t):
        """
        wrapper for evolve state__p_deriv :
        Compute derivative of localization probability according to Schrodinger equation
        Start state is localized on the grahp start site
        Only returns resulting projection on the target site
        """
        return self.evolve_state_p_deriv(gr, gr.get_start_state(), t)[gr.target, :]

class QutipSESolver(SESolver):
    """
    Schrodinger Equation solver on a QWGraph that uses Qutip library
    """
    
    def __init__(delf):
        pass

    def evolve_state( self, gr : QWGraph, psi, t):
        H = gr.get_h()
        psi = qt.Qobj(psi)

        res = qt.sesolve(H, psi, t, [])     

        return res.states

    #todo check correctness
    def evolve_state_deriv(self, gr : QWGraph, psi,t):
        
        psi_t = self.evolve_state(gr, psi, t)
        psi_t = qt.Qobj(psi_t)
        return -1j * gr.get_h() * psi_t 

    def evolve_state_p(self, gr, psi, t) :
        
        H = gr.get_h()
        psi = qt.Qobj(psi)

        E_list = []
        for i in range(gr.N):
            E_list.append(gr.get_projector(i))

        res = qt.sesolve(H, psi, t, E_list)

        return res.expect

    def evolve_state_p_deriv(self, gr, psi, t):
        
        #qutip has problems with non list input
        #and with evaulation times not starting with 0
        t, strip = format_qutip_time(t)

        H = self.gr.get_h()
        psi = qt.Qobj(psi)

        E_prime = -1j * qt.commutator( self.gr.get_projector(), H)

        res = qt.sesolve(H, psi, t, [E_prime])

        if strip:
            return res.expect[0][1:]
        else:
            return res.expect[0]

    def evolve_default(self, gr : QWGraph, t):
        return   self.evolve_state(gr, gr.get_start_state(), t)[gr.target, :]
    
    def evolve_default_deriv(self, gr : QWGraph, t):
        return self.evolve_state_deriv(gr, gr.get_start_state(), t)[gr.target, :]

    def evolve_default_p(self, gr : QWGraph, t):
        return self.evolve_state_p(gr, gr.get_start_state(), t)[gr.target, :]

    def evolve_default_p_deriv(self, gr : QWGraph, t):
        return self.evolve_state_p_deriv(gr, gr.get_start_state(), t)[gr.target, :]



#############################################################



class Analyzer(object):
    """
    Class that wraps around a QWGraph and extracts/computes information on its evolution
    The criteria on which those information are extracted depends on its state

    In particular one can vary:

    Best maximum criterium:
    -TC : best maximum for a given choice of phases is the result of a optimization routine on a linearly scaling window
    -first : best maximum for a given phase choice is the first maximum encountered

    Best phase choice algorithm:
    -none : no phase optimization
    -smart : choose the best among multiple of pi/2
    -fix : no phase optimization (but arbitraty choice of a fixed phase stored in fix_phi)

    Solver backend:
    -eigen :  SESolver is a EigenSESolver
    -qutip : SESolver is QutipSESolver

    """

    _modes = ["TC", "first"]
    _opt_modes = ["none", "smart", "fix"]
    _solver_modes = ["eigen", "qutip"]

    def __init__(self, gr = QWGraphBuilder.Line(2), event_s = 1, TC = 1, mode = "TC", opt_mode = "none",solver_mode = "eigen", diag = True):

        
        if not solver_mode in Analyzer._solver_modes:
            raise ValueError("Solver mode not supported/recognized")
        else :
            if solver_mode == "qutip" :
                self.solver = QutipSESolver()
            else :
                self.solver = EigenSESolver()

        if not mode in Analyzer._modes:
            raise ValueError("Mode not supported/recognized")
        else :
            self.mode = mode

        if not opt_mode in Analyzer._opt_modes:
            raise ValueError("Opt mode not supported/recognized")
        else :
            self.opt_mode = opt_mode

        self.event_size = event_s
        self.TIME_CONSTANT = TC
        
        self.fix_phi = None
        self.diag = diag

        self._graph = gr 

    def evolution_grid(self, phase_vec = None, bounds =(0,10), step = .1):
        """
        Return the evolution of localization probability on the target state
        Sampled on a mesh of points defined by bounds and step

        #Arguments
        phase_vec : list(float) -> phase choice for the rephasing links
        bounds : (float,float) ->
        step : float -> interval that separates sample points

        #Returns
        list[float] -> the evolution sampled on the mesh
        """

        if phase_vec == None:
            phase_vec = np.repeat(0,self.dim())
        
        self.rephase_graph(phase_vec)

        sample =np.arange(bounds[0], bounds[1], step)

        return self.solver.evolve_default_p(self._graph, sample)

    def locate_max(self):
        """
        Localize the "Best Maximum" according to the chosen criteria and phases

        wrapper for _locate_max_TC and _locate_max_first
        """

        if self.mode == "TC" :
            return self._locate_max_TC()

        elif self.mode == "first" :
            return self._locate_max_first()

    def _locate_max_TC(self):
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
        end = self.max_search_time()
        def f(t):
            return -1*self.solver.evolve_default_p( self._graph, t)

        def f_prime(t):
            return -1*self.solver.evolve_default_p_deriv( self._graph, t)

        b_vec = np.linspace(start,end, (end-start)//self.event_size + 2)
        sol_vec = np.empty( len(b_vec)-1)
        
        for i in range( len(b_vec)-1):
            
            sol = opt.minimize(f, \
                                x0 = (b_vec[i]+b_vec[i+1])/2, \
                                bounds = [(b_vec[i],b_vec[i+1])], \
                                jac = f_prime)
            
            sol_vec[i] = sol.x[0]
        
        probs = self.solver.evolve_default_p( self._graph, sol_vec)
        return ( sol_vec[np.argmax(probs)], max(probs))

    def _locate_max_first(self):
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
        end = self.get_gr().distance()/4
        exp_scale = 1.5
        evt_sample_scale = .1
        sign_change_safe = -1e-24
        maxiter = 10

        res = None
        for i in range(maxiter) :
            #print("looking for first max in "+ str(start) +" - "+ str(end) )

            sample = np.linspace(start, end, int((end-start)//(self.event_size* evt_sample_scale))+2)

            #print(sample)
            deriv_evo = self.solver.evolve_default_p_deriv( self._graph, sample)

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
                            return self.solver.evolve_default_p_deriv(self._graph, t)
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

        return (root, self.solver.evolve_default_p( self._graph, root)[0])

    def deriv_roots(self):
        """
        Compute a list of the position of "all" the stationary point of the probability evolution

        #todo: implement it??
        """
        
        pass
    
    def performance(self, phi_vec = None, target = "p"):
        """
        Evaluate the transport performance for a given phases configuration

        #Returns

        float -> loc probability for best max (target == "t") of transport time of best max (target == "p")
        """


        if not np.any(phi_vec):
            phi_vec = np.repeat(0,self.dim())

        self.rephase_graph(phi_vec)

        if target == "p":
            return self.locate_max()[1]
        elif target == "t":
            return self.locate_max()[0]

    def performance_diag(self, phi, target = "p"):
        """
        Evaluate the diag transport performance for a given diag phase

        #Returns

        float -> loc probability for best max (target == "t") of transport time of best max (target == "p")
        """

        #print(phi)

        self.rephase_graph( np.repeat( phi, self.dim() ))
        
        if target == "p":
            return self.locate_max()[1]
        elif target == "t":
            return self.locate_max()[0]
                
            
    def performance_grid(self, sample_step = 100, target = "p"):
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
         
        n_sample = [sample] * self.dim()
        grid = np.meshgrid( *n_sample)
            
        out_shape = [sample_step] * self.dim()
        out = np.empty(sample_step**self.dim())
        out = np.reshape(out, out_shape)

        it = np.nditer(out, flags=['multi_index'])
        for val in it:

            i = it.multi_index
            
            phi_vec = np.empty(self.dim())
            for j in range(self.dim()):
                phi_vec[j] = grid[j][i]

            self.rephase_graph(phi_vec)
            out[i] = self.locate_max()[target]

        return out

    def performance_grid_diag(self, sample_step = 100, target = "p"):
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
            self.rephase_graph( np.repeat( sample[i], \
                                               self.dim() ))

            out[i] = self.locate_max()[target]

        return out

    def performance_best(self, target = "p"):
        """
        Compute the best transport performance according to current time and phase criteria

        #Returns

        float -> loc probability for best max (target == "t") of transport time of best max (target == "p")
        """
        best_phi = 0

        if self.opt_mode == "none":
            best_phi = 0
        elif self.opt_mode == "smart" :
            best_phi = self.optimum_phase_smart()[0]
        elif self.opt_mode == "fix" :
            assert np.any(self.fix_phi )
            best_phi = self.fix_phi
        else :
            best_phi = self.optimum_phase_minimize()[0]

        #print(best_phi, self.opt_mode)

        if self.diag:
            return self.performance_diag(best_phi, target)
        else:
            return self.performance(best_phi, target)

    def optimum_phase_yolo(self, step = 100):
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

        sample = phase_sample(step)

        if self.diag:
            perf = self.performance_grid_diag(step)
        else :
            perf = self.performance_grid(step)

        pos_max = np.argmax(perf)

        return pos_max/step*2*np.pi, perf[pos_max]
    
    def optimum_phase_minimize(self):
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

        #minimize the inverse of perf function
        if self.diag:
            def perf(x):
                return -1*self.performance_diag(x)

            #todo: add a smart way to figure out the number of sections

            sol = section_minimize(perf, bounds = [0, 2*np.pi])

        else:
            def perf(x):
                return -1* self.performance(x)

            sol = opt.minimize(perf, \
                       x0 = np.repeat(.5, self.dim())   , \
                       bounds = [(0, 2*np.pi)]* self.dim() )

        return sol["x"], -1*perf(sol["x"])

    #todo: does it work???
    def optimum_phase_smart(self):
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

        if self.dim() == 0 :
            return 0 , self.performance()

        if self.diag:
            sample = phase_sample(step = 5)

            out = 0
            best = 0

            for phase in sample:

                diag_phase = np.repeat(phase, self.dim()) 
                self.rephase_graph(diag_phase)
                cur = self.locate_max()[1]

                if cur > best:
                    best = cur
                    out = phase

            return out, self.performance_diag( out)

        sample = phase_sample(step = 5)
        n_sample = [sample] * self.dim()
        grid = np.meshgrid( *n_sample)

        best = 0
        out = np.repeat(0, self.dim())

        it = np.nditer(grid[0], flags=['multi_index'])
        for val in it:

            i = it.multi_index
            
            phi_vec = np.empty(self.dim())
            for j in range(self.dim()):
                phi_vec[j] = grid[j][i]

            self.rephase_graph(phi_vec)
            cur = self.locate_max()[1]

            if cur > best:
                best = cur
                out = phi_vec

        return out, self.performance( out)
    
    def get_TC(self):
        return self.TIME_CONSTANT

    def set_TC(self, TC):
        self.TIME_CONSTANT = TC
        
    def max_search_time(self):
        """
        Return the time window size of the TC criterion
        """
        return self._graph.distance() * self.TIME_CONSTANT

    def get_gr(self):
        return self._graph

    def set_gr(self, gr):
        self._graph = gr


    def rephase_graph(self, phi_vec):
        self._graph.rephase(phi_vec)
        self._graph.update_eigen()

    def dim(self):
        """
        Number of free phases of the analyzed graph instance
        """
        return self.get_gr().get_phase_n()

    def get_label(self, mode = True):
        """
        Generate an automatic plot label that represent the internal state of the Analyzer
        """

        if not mode:
            return self.get_gr().code

        mode_label = self.mode

        if self.mode == "TC":
            mode_label = chr(957) + "=" + str( self.get_TC() )

        return self.get_gr().code + " " + mode_label

    def set_fix_phi(self, phi):
        self.fix_phi = phi

    def get_fix_phi(self):
        return self.fix_phi

    def get_solver(self):
        return self.solver

    def set_solver(self, solver):
        self.solver = solver

    def get_mode(self):
        return self.mode

    def set_mode(self, mode):
        self.mode = mode

    def get_opt_mode(self):
        return self.opt_mode

    def set_opt_mode(self, opt_mode):
        self.opt_mode = opt_mode

    def get_diag(self):
        return self.diag

    def set_diag(self, diag : bool):
        self.diag = diag




        
        

