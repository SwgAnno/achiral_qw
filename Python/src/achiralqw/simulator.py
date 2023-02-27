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


#todo: differentiate QutipSESolver and EigenSESolver


#general purpose class to handle the time evolution computation of the QW
class SESolver(object):

    def __init__(self, gr, qutip = False):
        self.gr = gr

        self.gr.update_eigen()
        self.set_qutip(qutip)

    
    # decompose and recompose localized state basis into matrix eigenvectors basis
    def decompose_localized(self, local_A) :
        return  self.gr.eig_vec.conj().T.dot(local_A)

    def recompose(self, eigen_A):
        return self.gr.eig_vec.dot(eigen_A)


    #return exponential map in eigenvector basis as column vector
    def exp_map(self,t):
        return np.exp(-1j * np.outer( self.gr.eig_val, t) )

    #evolve state in the localized basis/ get evolved probability amplitudes
    def evo_psi(self, psi, t):

        #carry on calculation in  H eigenvector basis
        A_t = self.decompose_localized(psi) *self.exp_map(t)

        return self.recompose( A_t)

    #evolve the desired state,
    #get results as sheer probability (p_psi)
    # or their respective derivatives (p_psi_prime)
    
    def evo_p_psi(self, psi, t):
        if self.qut:
            H = self.gr.get_h()

            E_list = []
            for i in range(self.gr.N):
                E_list.append(self.gr.get_projector(i))

            res = qt.sesolve(H, psi, t, E_list)

            return res.expect
        else:
            return np.abs(self.evo_psi(psi, t))**2

    def evo_p_psi_prime(self, psi, t):
   
        A_t = self.decompose_localized(psi)*self.exp_map(t)

        A_t_prime = self.gr.eig_val* A_t

        #get probability derivative as 2Re[(a)(a*)']
        out = self.recompose(A_t) * np.conjugate(self.recompose(A_t_prime))

        return -2* np.imag(out)


    #helper function to do optimization on
    #get localization probability on target site or its derivative
    # with the eigenvalue evolution method (p/p_prime_old)
    # or Qutip solvers (p/p_prime_qut)
    
    def target_p_old(self, t):
        return   self.evo_p_psi(self.gr.get_start_state(), t)[self.gr.target,:]
    
    def target_p_qut(self, t):
        #qutip has problems with non list input
        #and with evaulation times not starting with 0
        t, strip = format_qutip_time(t)
        
        psi_0 = self.gr.get_start_state(qut = True)
        H = self.gr.get_h()
        E = self.gr.get_projector()


        res = qt.sesolve(H, psi_0, t, [E])

        if strip:
            return res.expect[0][1:]
        else:
            return res.expect[0]

    def target_p_prime_old(self, t):
        return self.evo_p_psi_prime(self.gr.get_start_state(), t)[self.gr.target,:]

    def target_p_prime_qut(self, t):
        
        #qutip has problems with non list input
        #and with evaulation times not starting with 0
        t, strip = format_qutip_time(t)
            
        psi_0 = self.gr.get_start_state(qut = True)
        H = self.gr.get_h()
        E_prime = -1j * qt.commutator( self.gr.get_projector(), H)

        res = qt.sesolve(H, psi_0, t, [E_prime])

        if strip:
            return res.expect[0][1:]
        else:
            return res.expect[0]
    
    #not very useful getter/setter but it makes the implementation transparent
    def get_gr(self):
        return self.gr

    def set_gr(self, gr):
        self.gr = gr
        self.gr.update_eigen()

    def get_qutip(self):
        return self.qut

    def set_qutip(self, qutip):
        self.qut = qutip
        
        if not qutip:
            self.target_p = self.target_p_old
            self.target_p_prime = self.target_p_prime_old
        else:
            self.target_p = self.target_p_qut
            self.target_p_prime =  self.target_p_prime_qut

    def rephase_gr(self, phi_vec):
        self.gr.rephase(phi_vec)
        self.gr.update_eigen()


class EigenSESolver(SESolver) :

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


        #carry on calculation in  H eigenvector basis
        A_t = EigenSESolver._to_eiegn_basis(gr.eig_vec, psi) * EigenSESolver._exp_map(gr.eig_val, t)

        return EigenSESolver._to_state_basis(gr.eig_vec, A_t)

    #todo check correctness
    def evolve_state_deriv(self, gr : QWGraph, psi,t):
        
        return -1j * gr.mat.dot(psi)

    def evolve_state_p(self, gr, psi, t) :
        
        return np.power(np.abs(self.evolve_state(gr, psi, t)), 2)

    def evolve_state_p_deriv(self, gr, psi, t):
           
        A_t = EigenSESolver._to_eiegn_basis(gr.eig_vec, psi) * EigenSESolver._exp_map(gr.eig_val, t)

        A_t_prime = gr.eig_val* A_t

        #get probability derivative as 2Re[(a)(a*)']
        out = EigenSESolver._to_state_basis(gr.eig_vec, A_t) * np.conjugate(EigenSESolver._to_state_basis( gr.eig_vec, A_t_prime))

        return -2* np.imag(out)

    def evolve_default(self, gr : QWGraph, t):
        return   self.evolve_state(gr, gr.get_start_state(), t)[gr.target, :]
    
    def evolve_default_deriv(self, gr : QWGraph, t):
        return self.evolve_state_deriv(gr, gr.get_start_state(), t)[gr.target, :]

    def evolve_default_p(self, gr : QWGraph, t):
        return self.evolve_state_p(gr, gr.get_start_state(), t)[gr.target, :]

    def evolve_default_p_deriv(self, gr : QWGraph, t):
        return self.evolve_state_p_deriv(gr, gr.get_start_state(), t)[gr.target, :]

class QutipSESolver(SESolver):
    
    def __init__(delf):
        pass

    def evolve_state( self, gr : QWGraph, psi, t):
        H = gr.get_h()
        psi = qt.Qobj(psi)

        res = qt.sesolve(H, psi, t, [])     

        return res.states

    #todo check correctness
    def evolve_state_deriv(self, gr : QWGraph, psi,t):
        
        psi = qt.Qobj(psi)
        return -1j * gr.get_h() * psi 

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

    _modes = ["TC", "first"]
    _opt_modes = ["none", "smart", "fix"]
    _solver_modes = ["eigen", "qutip"]

    def __init__(self, gr = QWGraphBuilder.Line(2), event_s = 1, TC = 1, qutip = False, mode = "TC", opt_mode = "none",solver_mode = "eigen", diag = True):

        
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

    #get evolution on target site
    def evo_full(self, phi_vec = None, bounds =(0,10), step = .1):

        if phi_vec == None:
            phi_vec = np.repeat(0,self.dim())
        
        self.rephase_graph(phi_vec)

        sample =np.arange(bounds[0], bounds[1], step)

        return self.solver.evolve_default_p(self._graph, sample)



    #get stationary points in probability evolution on the target stite
    def deriv_roots(self):
        
        pass

    def locate_max(self):

        # minimize works fine for TC near 1 (aka few maxima actually there)
        # Not reliable for higher TC
        # therefore we split the search interval into
        # smaller pieces according to a reference "event size"
        
        if self.mode == "TC":      
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

        #digits are definitely note carefully researched, a check is needed
        if self.mode == "first" :
            start = .001
            end = self._graph.distance()/4
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

    #wrapper for locate_max() with desired phases
    def performance(self, phi_vec = None, t = False):

        if not np.any(phi_vec):
            phi_vec = np.repeat(0,self.dim())

        self.rephase_graph(phi_vec)

        if not t:
            return self.locate_max()[1]
        else :
            return self.locate_max()[0]

    #wrapper for locate_max() with desired equal phase
    def performance_diag(self, phi, t= False):

        #print(phi)

        self.rephase_graph( np.repeat( phi, self.dim() ))
        
        if not t:
            return self.locate_max()[1]
        else :
            return self.locate_max()[0]
                
    #get full sampled performance as a ndarray
    def performance_full(self, sample_step = 100, target = "p"):
        
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

    #equal phases setting transport performance
    def performance_full_diag(self, sample_step = 100, target = "p"):

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

    #brute force search for best phase
    def optimum_phase_yolo(self, step = 100):

        sample = phase_sample(step)

        if self.diag:
            perf = self.performance_full_diag(step)
        else :
            perf = self.performance_full(step)

        pos_max = np.argmax(perf)

        return pos_max/step*2*np.pi, perf[pos_max]
    
    #use scipy minimize for optimum phase
    def optimum_phase_minimize(self):

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

    #check for just +-1 +-i
    #todo: does it work???
    def optimum_phase_smart(self):

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

    def evaluate(self, target = "p"):
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
        target_t = (target != "p")

        #todo : t is boolean???????
        if self.diag:
            return self.performance_diag(best_phi, t = target_t)
        else:
            return self.performance(best_phi, t = target_t)
        

    def max_search_time(self):
        return self._graph.distance() * self.TIME_CONSTANT
    
    def get_TC(self):
        return self.TIME_CONSTANT

    def set_TC(self, TC):
        self.TIME_CONSTANT = TC
        
    def get_gr(self):
        return self._graph

    def set_gr(self, gr):
        self._graph = gr


    def rephase_graph(self, phi_vec):
        self._graph.rephase(phi_vec)
        self._graph.update_eigen()

    def set_fix_phi(self, phi):
        self.fix_phi = phi

    def get_fix_phi(self):
        return self.fix_phi

    #get qutip status
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

    def get_fix_phi(self):
        return self.fix_phi

    def set_fix_phi(self, phi):
        self.fix_phi = phi

    def get_diat(self):
        return self.diag

    def set_diag(self, diag : bool):
        self.diag = diag

    #number of free phases in the graph instance
    def dim(self):
        return self.get_gr().get_phase_n()

    #get automatic plot label according to target graph and analysis mode

    def get_label(self, mode = True):

        if not mode:
            return self.get_gr().code

        mode_label = self.mode

        if self.mode == "TC":
            mode_label = chr(957) + "=" + str( self.get_TC() )

        return self.get_gr().code + " " + mode_label



        
        

