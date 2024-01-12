from abc import abstractmethod
from achiralqw.graph import QWGraph
import numpy as np
import qutip as qt
    

#qutip has problems with non list input
#and with evaulation times not starting with 0
def format_qutip_time( t):
    if type(t) == float:
        t = [t]

    if t[0] != 0:
        if isinstance(t, (np.ndarray, np.generic)):
            out = np.insert(t,0,0.)
        else:
            out = t.copy()
            out.insert(0,0.)
        
        return out,True
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
    
    def __init__(self):
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
        
        #qutip has problems with non list input
        #and with evaulation times not starting with 0
        t, strip = format_qutip_time(t)

        H = gr.get_h()
        psi = qt.Qobj(psi)

        E_list = []
        for i in range(gr.N):
            E_list.append(gr.get_projector(i))

        res = qt.sesolve(H, psi, t, E_list)

        #print(t, res.expect)
        if strip:
            return [res.expect[i][1:] for i in range(gr.N)]
        else:
            return res.expect

    def evolve_state_p_deriv(self, gr, psi, t):
        
        #qutip has problems with non list input
        #and with evaulation times not starting with 0
        t, strip = format_qutip_time(t)

        H = gr.get_h()
        psi = qt.Qobj(psi)

        E_prime_list = []
        for i in range(gr.N):
            E_prime_list.append(-1j * qt.commutator( gr.get_projector(i), H))

        res = qt.sesolve(H, psi, t, E_prime_list)

        if strip:
            return [res.expect[i][1:] for i in range(gr.N)]
        else:
            return res.expect

    def evolve_default(self, gr : QWGraph, t):
        return   self.evolve_state(gr, gr.get_start_state(), t)[gr.target, :]
    
    def evolve_default_deriv(self, gr : QWGraph, t):
        return self.evolve_state_deriv(gr, gr.get_start_state(), t)[gr.target, :]

    def evolve_default_p(self, gr : QWGraph, t):
        return self.evolve_state_p(gr, gr.get_start_state(), t)[gr.target]

    def evolve_default_p_deriv(self, gr : QWGraph, t):
        return self.evolve_state_p_deriv(gr, gr.get_start_state(), t)[gr.target]



#############################################################
    
def evolution_grid(gr, solver = EigenSESolver(), phase_vec = None, bounds =(0,10), step = .1):
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
        phase_vec = np.repeat(0, gr.get_phase_n())

    gr.rephase(phase_vec)

    sample =np.arange(bounds[0], bounds[1], step)

    return solver.evolve_default_p(gr, sample)




        
        

