from Graph import *
import numpy as np

class SESolver(object):

    def __init__(self, gr):
        self.gr = gr

    
    # decompose and recompose localized state basis into matrix eigenvectors basis
    def decompose_localized(self, local_A) :
        return  self.gr.eig_vec.conj().T.dot(local_A)

    def recompose(self, eigen_A):
        return self.gr.eig_vec.dot(eigen_A)


    #return exponential map in eigenvector basis
    def exp_map(self,t):
        return np.exp(-1j * np.outer( self.gr.eig_val, t) )

    #evolve state in the localized basis/ get evolved probability amplitudes
    def evo_l_t(self, l, t):
        A_t = np.diag(self.decompose_localized(l)).dot(self.exp_map(t))

        return self.recompose(A_t)
    
    def evo_p_l_t(self, l, t):
        return np.abs(self.evo_l_t(l, t))**2
        

