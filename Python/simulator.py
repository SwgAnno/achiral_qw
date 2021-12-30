from Graph import *
import numpy as np
from scipy import optimize as opt

#helper function to get a discrete sample of modulus 1 complex numbers
def phase_sample(step = 100):
    sample = np.linspace(0, np.pi*2, step)
    return np.exp(1j * sample)

class SESolver(object):

    def __init__(self, gr):
        self.gr = gr

    
    # decompose and recompose localized state basis into matrix eigenvectors basis
    def decompose_localized(self, local_A) :
        return  self.gr.eig_vec.conj().T.dot(local_A)

    def recompose(self, eigen_A):
        return self.gr.eig_vec.dot(eigen_A)


    #return exponential map in eigenvector basis as column vector
    def exp_map(self,t):
        return np.exp(-1j * np.outer( self.gr.eig_val, t) )

    #evolve state in the localized basis/ get evolved probability amplitudes
    def evo_l(self, l, t):

        #carry on calculation in  H eigenvector basis
        A_t = self.decompose_localized(l) *self.exp_map(t)

        return self.recompose( A_t)
    
    def evo_p_l(self, l, t):
        return np.abs(self.evo_l(l, t))**2

    def deriv_p_l(self, l, t):
   
        A_t = self.decompose_localized(l)*self.exp_map(t)

        A_t_prime = self.gr.eig_val* A_t

        #get probability derivative as 2Re[(a)(a*)']
        out = self.recompose(A_t) * np.conjugate(self.recompose(A_t_prime))

        return -2* np.imag(out)


    #helper function to do optimization on
    def target_p(self, t):
        return   self.evo_p_l(self.gr.get_start_state(), t)[self.gr.target,:]

    def target_p_prime(self, t):
        return self.deriv_p_l(self.gr.get_start_state(), t)[self.gr.target,:]
    
    #not very useful getter/setter but it makes the implementation transparent
    def get_gr(self):
        return self.solver.gr

    def set_gr(self, gr):
        self.gr = gr

    def rephase_gr(self, phi_vec):
        self.gr.rephase(phi_vec)
        


class Analyzer(object):

    def __init__(self, gr, event_s = 2, TC = 1):
        self.solver = SESolver(gr)

        self.event_size = event_s
        self.TIME_CONSTANT = TC

    #get stationary points in probability evolution on the target stite
    def deriv_roots(self):
        
        pass

    def locate_max(self, mode = "TC"):

        # minimize works fine for TC near 1 (aka few maxima actually there)
        # Not reliable for higher TC
        # therefore we split the search interval into
        # smaller pieces according to a reference "event size"
        
        if mode == "TC":      
            start = 0
            end = self.max_search_time()
            def f(t):
                return -1*self.solver.target_p(t)
            def f_prime(t):
                return -1*self.solver.target_p_prime(t)

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

    def performance(self, sample_step = 100):
        
        sample = phase_sample(sample_step)

        #todo: add support for n dimensionss

        #np.meshgrid

        out = np.empty(sample_step)
        for i in range(len(sample)):
            self.solver.rephase_gr( [sample[i]])
            out[i] = self.locate_max()[1]

        return out

    #brute force search for best phase
    def optimum_phase_yolo(self, step = 100)

        sample = phase_sample(step)
        perf = self.performance(step)

        pos_max = np.argmax(perf)

        return pos_max/step*2*np.pi

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

        
#########

if __name__ == "__main__" :
    a = QWGraph.Ring(6)

    test = Analyzer(a)

    print(test.locate_max())


        
        

