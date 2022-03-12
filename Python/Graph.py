import numpy as np
from numpy.linalg import eigh
import igraph as ig
import qutip as qt

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
    def get_start_state(self, qut = False):
        return self.basis(self.start, qut)

    def update_eigen(self):
        self.eig_val, self.eig_vec = eigh(self.mat)

        self.eig_val = np.reshape(self.eig_val, (self.N,1))
        

    #Trace elements manipulation utilities
    def retrace_E(self, E):
        for i in range(self.N):
            self.mat[i][i] = E

        self.update_eigen()

    def retrace_conn(self):
        ref = self.to_igraph()
        
        for i in range(self.N):
            self.mat[i][i] = ref.degree(i)
            print(ref.degree(i))

    #utiliry method to save trace elements in a buffer
    #since converting from and to igraph instances
    #does not presereve that information
    def buffer_trace(self):
        out = [ self.mat[i][i] for i in range(self.N)]

        return out

    def retrace(self, T):
        for i in range(self.N):
            self.mat[i][i] = T[i]

        self.update_eigen()


    #assign a new value to the target rephasing links
    def rephase(self, phi = [1j]) :
        if( len( self.re_coord) != len(phi)):
            print("rephase() error: wrong number of phases given")
            return()

        for i in range(len(phi)) :
            p = self.re_coord[i]
            self.mat[p[0]][p[1]] = -1*phi[i]
            self.mat[p[1]][p[0]] = -1*np.conjugate(phi[i])

        self.update_eigen()

        
    #concatenate two graph creating the link
    #between first end site and second start site
    #(+ operator)
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

    
    #concatenate two graph merging the first end site and the second start site
    #( | [or] operator)
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

    #concatenate multiple graph units into a chain ( * operator)
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

    #add ("cut") the edges specified in the vector as tuples
    def cut(self, cut_vec):
        ref = self.to_igraph()
        ref.add_edges(cut_vec)

        temp = self.buffer_trace()

        out = QWGraph.from_igraph(ref, ends = (self.start, self.target))
        out.retrace(temp)

        return out

    #get automatically a new set of phased links
    #according to the spanning tree of the graph       
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

    #get a visual rapresentation of the graph (relies on igraph)
    def plot(self) :

        ref = self.to_igraph()

        names = []
        for i in range(self.N):
            names.append(str(i))

        cols = []
        for e in ref.es:
            if e.tuple in self.re_coord :
                cols.append("red")
            else :
                cols.append("black")

        v_cols = ["yellow"]* self.N
        v_cols[self.start] = "green"
        v_cols[self.target] = "red"

        ref.vs["label"] = names
        ref.vs["color"] = v_cols
        ref.es["color"] = cols
        ig.plot( ref , layout = ig.Graph.layout_fruchterman_reingold(ref))


    #Ring graph constructor
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

    #Multi path element graph constructor
    def Parallel(paths , p_len, E = 2):
        ref = ig.Graph()

        N = 2 + paths* (p_len-1)

        ref.add_vertices(N)

        for i in range(paths):
            ref.add_edge(0, i+1)
            ref.add_edge(N-2-i, N-1)

        for i in range(p_len-2):
            for m in range(paths):
                ref.add_edge(1 + paths*i + m,1 + paths*(i+1) + m)

        out = QWGraph.from_igraph(ref)
        out.retrace_E(E)
        out.update_eigen
        
        return out
        
    #get numpy basis vector
    def basis(self, i, qut = False):
        if qut:
            return qt.basis(self.N,i)
        
        out = np.zeros((self.N,1), dtype = complex)
        out[i]= 1

        return out

    #get QuTip projector operator on the Nth site
    def get_projector(self, i = None):
        if not i:
            i = self.target
            
        out_mat = np.zeros((self.N,self.N))
        out_mat[i,i] = 1

        out = qt.Qobj(out_mat)

        return out

    #get the hamiltonian of the system as a QuTip object
    def get_h(self):
        return qt.Qobj(self.mat)

    #get the number of registered free phases
    def get_phase_n(self):
        return len(self.re_coord)


############################

if __name__ == "__main__":
    a = QWGraph.Parallel(4,3)

    a.plot()

