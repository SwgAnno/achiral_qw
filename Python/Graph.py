import numpy as np
from numpy.linalg import eigh
import igraph as ig

#get numpy basis vector

def basis(N,l):
    out = np.zeros((N,1), dtype = complex)
    out[l]= 1

    return out

#QW simulation oriented graph class

class QWGraph(object) :

    def __init__(self, N = 4) :

        self.N = N
        self.code = "e"
        
        self.mat = np.identity(N, dtype = complex)

        #it'd better be not mandatory
        self.update_eigen()
        
        self.re_coord = []
        self.start = 0
        self.target = 0

    #return localized start state for evolution
    def get_start_state(self):
        return basis(self.N, self.start)

    def update_eigen(self):
        self.eig_val, self.eig_vec = eigh(self.mat)

        self.eig_val = np.reshape(self.eig_val, (self.N,1))
        

    def retrace_E(self, E):
        for i in range(self.N):
            self.mat[i][i] = E

        self.update_eigen()

    def retrace_conn(self):
        ref = self.to_igraph()
        
        for i in range(self.N):
            self.mat[i][i] = ref.degree(i)
            print(ref.degree(i))

    def rephase(self, phi = [1j]) :
        if( len( self.re_coord) != len(phi)):
            print("rephase() error: wrong number of phases given")

        for i in range(len(phi)) :
            p = self.re_coord[i]
            self.mat[p[0]][p[1]] = -1*phi[i]
            self.mat[p[1]][p[0]] = -1*np.conjugate(phi[i])

        self.update_eigen()

    #ring graph constructor

    def join_link(self, other):

        out_N = self.N + other.N
        out = QWGraph(out_N)

        out.code = self.code + "+" + other.code

        #copy old adj mats
        out.mat[0:self.N, 0:self.N] = self.mat
        out.mat[self.N:out.N, self.N:out.N] = other.mat

        #add new link

        out.mat[self.target, self.N + other.start] = -1
        out.mat[self.N + other.start, self.target] = -1

        out.update_eigen()
        out.compute_re_coord()
        out.start = self.start
        out.target = self.N + other.target

        #todo add r_coord shift

        return out

    def __add__(self, other) :
        return QWGraph.join_link(self, other)

    def join_nolink(self, other):

        out_N = self.N + other.N -1
        out = QWGraph(out_N)

        out.code = self.code + "|" + other.code

        #copy first mat
        out.mat[0:self.N, 0:self.N] = self.mat

        #copy 2nd mat skipping start site
        mat2 = np.delete( other.mat, other.start,0 )
        mat2 = np.delete( mat2, other.start,1 )
        
        out.mat[self.N:out.N, self.N:out.N] = mat2

        #join start site information
        j_row = np.delete( other.mat, other.start,1)[other.start, :]
        j_col = np.delete( other.mat, other.start,0)[:, other.start]

        out.mat[self.target, self.N:out.N] = j_row
        out.mat[self.N:out.N, self.target] = j_col
        

        out.update_eigen()
        out.compute_re_coord()
        out.start = self.start
        out.target = self.N + other.target
        if(other.target > other.start):
            out.target -= 1

        #todo add r_coord shift

        return out

    def __or__(self, other) :
        return QWGraph.join_nolink(self, other)

    def chain(self, rep, space = 0, HANDLES = True):

        if HANDLES:
            out = QWGraph.Line(1) + self
        else:
            out = QWGraph.Line(0) + self

        for i in range(rep-1):
            out = out | (QWGraph.Line(space)+self)

        if HANDLES:
            out = out + QWGraph.Line(1)

        return out

    def __mul__(self, rep) :
        return self.chain(rep, HANDLES = False)

    def add_handles(self, size, mode = "both", fix =0):

        if mode == "both":
            dx = size
            sx = size
        elif mode == "l":
            dx = 0
            sx = size
        elif mode == "r":
            dx = size
            sx = 0
        elif mode == "fixl":
            dx = size
            sx = fix
        elif mode == "fixr":
            dx = fix
            sx = size

        return QWGraph.Line(sx) + self + QWGraph.Line(dx)

    def reverse(self ):
        self.start, self.target = self.target, self.start
        
    def compute_re_coord(self) :

        ref = self.to_igraph()
        tree = ref.spanning_tree()

        phase = ref-tree

        n_re_coord = []
        for e in phase.es :
            n_re_coord.append( e.tuple)

##        names = []
##        for i in range(self.N):
##            names.append(str(i))
##
##        tree.vs["label"] = names
##        ref.vs["label"] = names
##        ig.plot(ref -tree)

        self.re_coord = n_re_coord
        

    #return QWGraph instance from igraph reference
    #the adj mat is obviously real
    def from_igraph( ig, E = 2, ends = None) :

        out = QWGraph( ig.vcount())

        #todo assign name
        #out.code = ig["name"]
        
        out.mat = np.array ( ig.get_adjacency().data, dtype = complex)
        out.retrace_E(E)

        out.update_eigen()
        out.compute_re_coord()

        if not ends :
            out.start = 0
            out.target = out.N-1
        else :
            out.start , out.target = ends

        return out

    #return igraph instance representing the QWGraph
    #the process does not tranfer phase info
    def to_igraph(self) :

        ref = np.zeros( (self.N, self.N))

        #format adjacency matix for igraph imput
        
        for i in range(self.N) :
            for m in range(self.N) :
                if self.mat[i][m] != 0 :
                    ref[i][m] = 1
            ref[i][i] = 0
        
        out = ig.Graph.Adjacency(ref, mode = "undirected")

        #print(ref)
        return out

    def plot(self) :

        ref = self.to_igraph()

        names = []
        for i in range(self.N):
            names.append(str(i))

        cols = []
        for e in ref.es:
            if e.tuple in self.re_coord :
                cols.append("green")
            else :
                cols.append("black")

        ref.vs["label"] = names
        ref.es["color"] = cols
        ig.plot( ref )

    def Ring(N, HANDLES = False, E = 2):
        out = QWGraph(N)
        out.code = "C"+ str(N)

        if N==0 :
            return(out)

        out.retrace_E(E)

        for i in range(N):
            out.mat[i][(i+1)%N] = complex(-1)
            out.mat[(i+1)%N][i] = complex(-1)

        #trace settings

        out.update_eigen()
        out.compute_re_coord()

        out.start = 0
        out.target = int(N/2)

        if HANDLES :
            return out.add_handles(1)

        return out

    #Line graph constructor

    def Line(N, E = 2):
        out = QWGraph(N)
        out.code = "L"+ str(N)

        if N==0 :
            return(out)

        out.retrace_E(E)

        for i in range(N-1):
            out.mat[i][(i+1)%N] = complex(-1)
            out.mat[(i+1)%N][i] = complex(-1)

        #trace settings

        out.update_eigen()

        out.start = 0
        out.target = N-1
        

        return out

