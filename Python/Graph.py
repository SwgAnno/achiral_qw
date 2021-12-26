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

    def __init__(self, *args, **kwargs) :

        super().__init__(*args, **kwargs)

        #print(args)
        #print(kwargs)

        #set N as a globally available variable,
        #due to the __init__ arguments we must let
        # the super contructor do the work
        self.N = self.vcount()
        self.code = "e"
        
        self.init_mat()

        #it'd better be not mandatory
        self.update_eigen()
        
        self.re_coord = []
        self.start = 0
        self.target = 0

    #transformed union and disjoint union
    #to operation which return a connected structure
        
    def disjoint_union(self, other) :
        joint = (self.target, other.start + self.N)
        out = super().disjoint_union(other)

        out.init_mat()
        diag_buff = np.append( np.diag(self.mat) , np.diag(other.mat))
        out.retrace_buff(diag_buff)
        out.update_eigen()

        out.start = self.start
        out.target = other.target + self.N

        return out

    def union(self, other, byname='auto') :

        self.name_vs()
        other.name_vs( self.N-1)
        other.name_v( other.start, str(self.target))

        out = ig.Graph.union(self,other, byname)

        out.init_mat()
        diag_buff = np.append( np.diag(self.mat) , np.diag(other.mat))
        out.retrace_buff(diag_buff)
        out.update_eigen() 

        out.start = self.start
        out.target = other.target + self.N -1

        return out

    def name_vs(self, bias = 0) :

        name_l = []

        for i in range(self.N):
            name_l.append(str(i + bias))

        self.vs["label"] = name_l
        self.vs["name"] = name_l

    def name_v(self, i, lab = "") :
        if lab == "":
            lab = str(i)

        ###todo: fix the mess
        temp =  self.vs["label"]
        temp[i] = lab

        self.vs["label"] = temp
        self.vs["name"] = temp

        print(self.vs["name"][i])

    #initialize adjacency matrix
    def init_mat(self) :
        self.mat = np.array ( self.get_adjacency().data, dtype = complex)*-1

    #return localized start state for evolution
    def get_start_state(self):
        return basis(self.N, self.start)

    def update_eigen(self):
        self.eig_val, self.eig_vec = eigh(self.mat)
        

    def retrace_E(self, E):
        for i in range(self.N):
            self.mat[i][i] = E

        self.update_eigen()

    #reload values on the main diagonal
    def retrace_buff(self, buff):
        for i in range(self.N):
            self.mat[i][i] = buff[i]

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
        out = QWGraph( n= N, edges = reference.get_edgelist())
        
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
        out = QWGraph(n = N, edges = reference.get_edgelist())
        
        code = "L"+ str(N)

        out.init_mat()
        out.retrace_E(E)

        #trace settings

        out.update_eigen()

        out.start = 0
        out.target = N-1
        

        return out

