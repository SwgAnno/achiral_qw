import numpy as np
from numpy.linalg import eigh
import igraph as ig
from igraph import Graph

#get numpy basis vector

def basis(N,l):
    out = np.zeros(N, dtype = complex)
    out[l]= 1

    return out

#QW simulation oriented graph class

class QWGraph(Graph) :

    def __init__(self, N = 4,edges=None) :

        super().__init__(N, edges)

        self.N = N
        self.code = "e"
        
        self.init_mat()

        #it'd better be not mandatory
        self.update_eigen()
        
        self.re_coord = []
        self.start = 0
        self.target = 0

    #initialize adjacency matrix
    def init_mat(self) :
        self.mat = np.array ( self.get_adjacency().data, dtype = complex)

    #return localized start state for evolution
    def get_start_state(self):
        return basis(self.N, self.start)

    def update_eigen(self):
        self.eig_val, self.eig_vec = eigh(self.mat)
        

    def retrace_E(self, E):
        for i in range(self.N):
            self.mat[i][i] = E

        self.update_eigen()

    """def retrace_conn(self):
        for i in range(N):
            count = 0
            
            self.mat[i][i] = E"""


    #ring graph constructor helper
    def Ring(N, E = 2):

        #it's horrible but it seems to be the only way to inherit from
        # igraph constructor helpers
        reference = Graph.Ring(N)
        out = QWGraph(N, reference.get_edgelist())
        
        out.code = "C"+ str(N)

        out.init_mat()
        out.retrace_E(E)
        
        #trace settings

        out.update_eigen()

        out.start = 0
        out.target = int(N/2)
        
        if N!= 1 :
            out.re_coord.append( (0,N-1))

        return out


    #Line graph constructor

    def Line(N, E = 2):

        # See Ring constructor helper
        reference = Graph.Lattice([N], circular = False)
        out = QWGraph(N, reference.get_edgelist())
        
        code = "L"+ str(N)

        out.init_mat()
        out.retrace_E(E)

        #trace settings

        out.update_eigen()

        out.start = 0
        out.target = N-1
        

        return out

