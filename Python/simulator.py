from Graph import *
import numpy as np
from scipy import optimize as opt
import qutip as qt

#helper function to get a discrete sample of modulus 1 complex numbers
def phase_sample(step = 100):
    sample = np.linspace(0, np.pi*2, step)
    return np.exp(1j * sample)

class SESolver(object):

    def __init__(self, gr, qutip = False):
        self.gr = gr
        self.qut = qutip

        if not qutip:
            self.target_p = self.target_p_old
            self.target_p_prime = self.target_p_prime_old
        else:
            self.target_p = self.target_p_qut
            #todo implement target_p_prime_qut
            self.target_p_prime = None

    
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
    def target_p_old(self, t):
        return   self.evo_p_psi(self.gr.get_start_state(), t)[self.gr.target,:]
    
    def target_p_qut(self, t):
        psi_0 = self.gr.get_start_state(qut = True)
        H = self.gr.get_h()
        E = self.gr.get_projector()

        res = qt.sesolve(H, psi_0, t, [E])

        return res.expect[0]

    def target_p_prime_old(self, t):
        return self.evo_p_psi_prime(self.gr.get_start_state(), t)[self.gr.target,:]
    
    #not very useful getter/setter but it makes the implementation transparent
    def get_gr(self):
        return self.solver.gr

    def set_gr(self, gr):
        self.gr = gr

    def rephase_gr(self, phi_vec):
        self.gr.rephase(phi_vec)
        


class Analyzer(object):

    def __init__(self, gr, event_s = 2, TC = 1, qutip = False, mode = "TC"):
        self.solver = SESolver(gr, qutip)

        self.event_size = event_s
        self.TIME_CONSTANT = TC
        self.mode = "TC"

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
                return -1*self.solver.target_p(t)

            #todo update with qut p_prime implementation
            if self.solver.target_p_prime :
                def f_prime(t):
                    return -1*self.solver.target_p_prime(t)
            else:
                f_prime = None

            b_vec = np.linspace(start,end, (end-start)//self.event_size + 2)
            sol_vec = np.empty( len(b_vec)-1)
            
            for i in range( len(b_vec)-1):
                
                sol = opt.minimize(f, \
                                   x0 = (b_vec[i]+b_vec[i+1])/2, \
                                   bounds = [(b_vec[i],b_vec[i+1])], \
                                   jac = f_prime)
                
                sol_vec[i] = sol.x[0]
            
            probs = self.solver.target_p(sol_vec)
            return ( sol_vec[np.argmax(probs)], max(probs))
        
        if self.mode == "first" :
            start = .01
            end = self.solver.gr.N/4
            exp_scale = 1.5
            evt_sample_scale = .1
            sign_change_safe = -1e-8

            found = False
            res = 0
            while not found :
##                print("looking for first max in "+ str(start) +" - "+ str(end) )

                sample = np.arange(start, end, self.event_size* evt_sample_scale)
                deriv_evo = self.solver.target_p_prime(sample)

                for i in range(len(sample)-1):
                    if deriv_evo[i]*deriv_evo[i+1] < sign_change_safe :
##                        print("Found deriv sign inversion at  " \
##                              + str(sample[i]) +" - " \
##                              + str(sample[i+1]) )
##                        print("   -Derivative values at extrama:  " \
##                              + str(deriv_evo[i]) +" - " \
##                              + str(deriv_evo[i+1]) )

                        res = opt.root_scalar( self.solver.target_p_prime, \
                                               bracket = [sample[i], \
                                                          sample[i+1]])
                        found = True
                        break
                        #find solution
                
                start, end = end, (end + (end-start)*exp_scale)

            return (res.root, self.solver.target_p(res.root)[0])

    #wrapper for locate_max() with desired phases
    def performance(self, phi_vec):

        p = np.exp(1j * phi_vec)
        self.solver.rephase_gr(p)

        return self.locate_max()[1]

    #wrapper for locate_max() with desired 
    def performance_diag(self, phi):

        p = np.exp(1j * phi)
        self.solver.rephase_gr( np.repeat( p, self.dim() ))

        return self.locate_max()[1]
                
    #get full sampled performance as a ndarray
    def performance_full(self, sample_step = 100):
        
        sample = phase_sample(sample_step)
         
        n_sample = [sample] * self.dim()
        grid = np.meshgrid( *n_sample)
            
        out_shape = [sample_step] * self.dim()
        out = np.empty(sample_step**self.dim())
        out = np.reshape(out, out_shape)

        it = np.nditer(out, flags=['multi_index'])
        for val in it:

            i = it.multi_index
            
            phi_vec = []
            for j in range(self.dim()):
                phi_vec.append(grid[j][i])

            self.solver.rephase_gr(phi_vec)
            out[i] = self.locate_max()[1]

        return out

    #equal phases setting transport performance
    def performance_full_diag(self, sample_step = 100):

        sample = phase_sample(sample_step)

        out = np.empty(sample_step)
        for i in range(len(sample)):
            self.solver.rephase_gr( np.repeat( sample[i], \
                                               self.dim() ))
            out[i] = self.locate_max()[1]

        return out

    #brute force search for best phase
    def optimum_phase_yolo(self, step = 100):

        sample = phase_sample(step)
        perf = self.performance_full(step)

        pos_max = np.argmax(perf)

        return pos_max/step*2*np.pi
    
    #use scipy minimize for optimum phase
    def optimum_phase_minimize(self, diag = False):

        #minimize the inverse of perf function
        if diag:
            def perf(x):
                return -1*self.performance_diag(x)
        else:
            def perf(x):
                return -1* self.performance(x)

        sol = opt.minimize(perf, \
                           x0 = np.repeat( 0, self.dim() )    , \
                           bounds = [(0, 2*np.pi)]* self.dim() )

        return sol.x

    #check for just +-1 +-i
    def optimum_phase_smart(self):

        sample = phase_sample(step = 5)
        
        n_sample = [sample] * self.dim()
        grid = np.meshgrid( *n_sample)

        best = 0
        out = np.repeat(0, self.dim())

        it = np.nditer(grid[0], flags=['multi_index'])
        for val in it:

            i = it.multi_index
            
            phi_vec = []
            for j in range(self.dim()):
                phi_vec.append(grid[j][i])

            self.solver.rephase_gr(phi_vec)
            cur = self.locate_max()[1]

            if cur > best:
                best = cur
                out = phi_vec

        return np.angle(out)

    def max_search_time(self):
        return self.solver.gr.N * self.TIME_CONSTANT
    
    def get_TC(self):
        return self.TIME_CONSTANT

    def set_TC(self, TC):
        self.TIME_CONSTANT = TC
        
    def get_gr(self):
        return self.solver.gr

    def set_gr(self, gr):
        self.solver.set_gr(gr)

    #number of free phases in the graph instance
    def dim(self):
        return self.get_gr().get_phase_n()
        

        
#########

if __name__ == "__main__" :
    a = QWGraph.Ring(6)

    c = QWGraph.chain(a,5)

    test = Analyzer(c, qutip = False)

##    print(test.locate_max(mode = "first"))
##    print(test.locate_max())

##    print(test.performance_full(sample_step = 5))

    print(test.optimum_phase_smart())



        
        

