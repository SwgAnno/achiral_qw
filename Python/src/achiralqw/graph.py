import numpy as np
from numpy.linalg import eigh
import igraph as ig
import qutip as qt
import matplotlib.pyplot as plt
import networkx as nx

#QW simulation oriented graph class

class QWGraph(object) :

    def __init__(self, N , mat = None, endpoints = None) :

        if mat is None :
            self.N = N
            self.code = "e"
            
            self.mat = np.zeros((N,N), dtype = complex)
        else :
            assert N == mat.shape[0]
            #mat is assumed to be a square numpy ndarray
            self.N = mat.shape[0]
            self.code = "m"

            self.mat = mat

        #it'd better be not mandatory
        
        self.compute_re_coord()

        if not endpoints == None :
            self.start = endpoints[0]
            self.target = endpoints[1]

            # ensure target site is the last diagonal entry
            if self.target != self.N -1 :
                #exchange row and columns
                self.mat[:,[self.target, self.N-1]] = self.mat[:,[self.N-1, self.target]]
                self.mat[[self.target, self.N-1],:] = self.mat[[self.N-1, self.target],:]

            if self.start  != 0:
                #revert change
                self.mat[:,[self.start, 0]] = self.mat[:,[0, self.start]]
                self.mat[[self.start, 0],:] = self.mat[[0, self.start],:]

        self.start = 0
        self.target = self.N-1

    def update_eigen(self):
        self.eig_val, self.eig_vec = eigh(self.mat)

        self.eig_val = np.reshape(self.eig_val, (self.N,1))
        

    #Trace elements manipulation utilities
    def retrace_E(self, E):
        for i in range(self.N):
            self.mat[i][i] = E

        

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

        


    #assign a new value to the target rephasing links
    #todo: check way to keep modulus constant
    def rephase(self, phi = [1j]) :
        
        if isinstance(phi, (list, np.ndarray)):
            if( len( self.re_coord) != len(phi)):
                print("rephase() error: wrong number of phases given")
                return()
        else :
            phi = [phi]

        for i in range(len(phi)) :
            p = self.re_coord[i]
            self.mat[p[0]][p[1]] = -1*phi[i]* np.abs(self.mat[p[0]][p[1]])
            self.mat[p[1]][p[0]] = np.conjugate(self.mat[p[0]][p[1]])

        

        
    #concatenate two graph creating the link
    #between first end site and second start site
    #(+ operator)
    def join_link(self, other):

        if self.N == 0 :
            return other
        if other.N == 0 :
            return self

        out_N = self.N + other.N
        out = QWGraph(out_N)

        out.code = self.code + "+" + other.code

        #copy old adj mats
        out.mat[0:self.N, 0:self.N] = self.mat
        out.mat[self.N:out.N, self.N:out.N] = other.mat

        #add new link

        out.mat[self.target, self.N + other.start] = -1
        out.mat[self.N + other.start, self.target] = -1

        
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

        outmat = np.zeros( (out_N, out_N), dtype= complex)

        # copy matrices so that target node of sel is the last diagonal entry
        if self.target != self.N -1 :
            #exchange row and columns
            self.mat[:,[self.target, self.N-1]] = self.mat[:,[self.N-1, self.target]]
            self.mat[[self.target, self.N-1],:] = self.mat[[self.N-1, self.target],:]

            outmat[0:self.N, 0:self.N] = self.mat

            #revert change
            self.mat[:,[self.target, self.N-1]] = self.mat[:,[self.N-1, self.target]]
            self.mat[[self.target, self.N-1],:] = self.mat[[self.N-1, self.target],:]
        else:
            outmat[0:self.N, 0:self.N] = self.mat

        #do the same thing with other but ensure start target is in first position
        if other.start != 0 :
            #exchange row and columns
            other.mat[:,[other.start, 0]] = other.mat[:,[0, other.start]]
            other.mat[[other.start, 0],:] = other.mat[[0, other.start],:]

            outmat[self.N-1: out_N, self.N-1: out_N] = other.mat

            #revert change
            other.mat[:,[other.start, 0]] = other.mat[:,[0, other.start]]
            other.mat[[other.start, 0],:] = other.mat[[0, other.start],:]
        else:
            outmat[self.N-1: out_N, self.N-1: out_N] = other.mat

    
        out = QWGraph(out_N, mat = outmat)
        out.compute_re_coord()
        out.code = self.code + "|" + other.code

        #account for the possible swap whte copying matrices
        out.start  = self.start                 if self.start != self.N-1   else self.target
        out.target = self.N + other.target-1    if other.target != 0        else self.N + other.start -1

        return out

    def __or__(self, other) :
        return QWGraph.join_nolink(self, other)

    #concatenate multiple graph units into a chain ( * operator)
    # speedup affects the couplings before adding the handles

    #todo: fix the fact that intra unit space does not get speedup
    def chain(self, rep, space = 0, speedup = 1, HANDLES = True):

        #create self copy with speedup

        su_mat = self.mat * speedup
        unit = QWGraph.Line(space, speedup = speedup) + QWGraph(self.N,mat = su_mat, endpoints = (self.start,self.target))
        
        #todo: df am i doing with P graph trace?
        temp = self.buffer_trace()
        unit.retrace(temp)

        new_code = self.code + "^{}".format(rep)
        out_N = (unit.N-1)*rep +1

        offset = 0
        if HANDLES:
            new_code = "h({})".format(new_code)

            outmat = np.zeros( (out_N+2,out_N+2), dtype= complex)
            outmat[0,1] = complex(-1)
            outmat[1,0] = complex(-1)
            outmat[-1,-2] = complex(-1)
            outmat[-2,-1] = complex(-1)

            offset = 1

        else :
            outmat = np.zeros((out_N,out_N) , dtype= complex)

        for i in range(rep):
            pos_0 = offset+ i*(unit.N-1)
            outmat[pos_0: (pos_0+ unit.N), pos_0: (pos_0+ unit.N)] = unit.mat

        out_N = out_N +2 if HANDLES else out_N
        out = QWGraph(out_N, outmat)
        out.compute_re_coord()
        out.code = new_code
        return out

    def __mul__(self, rep) :
        return self.chain(rep, HANDLES = False)

    def add_handles(self, size, mode = "both", fix =0):
        
        temp_code = self.code

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

        out = QWGraph.Line(sx) + self + QWGraph.Line(dx)
        out.code = "h(" + temp_code + ")"

        return out

    def reverse(self ):
        self.start, self.target = self.target, self.start

    #add ("cut") the edges specified in the vector as tuples
    def cut(self, cut_vec):

        #print(self.code, cut_vec)

        if type(cut_vec) is  tuple:
            cut_vec = [cut_vec]

        ref = self.to_igraph()
        ref.add_edges(cut_vec)

        temp = self.buffer_trace()

        out = QWGraph.from_igraph(ref, ends = (self.start, self.target))
        out.retrace(temp)

        return out

    def speedup(self, su):
        temp = self.buffer_trace()

        self.mat = self.mat * su
        self.retrace(temp)

    #get the equivalent graph in the krylov subspace relative to a given start state
    def krylov_transform( self, start_state = None):

        if start_state == None:
            start_state = self.get_start_state()

        #use Laczos algorithm to compute the new basis
        k_basis = []
        k_E = []
        k_A = []

        k_basis.append( start_state)
        k_A.append(0)
        l = 1
        for i in range(self.N):
            v = np.matmul( self.mat, k_basis[i])

            E = np.vdot(k_basis[i], v )

            v_orto = v - E* k_basis[i] - k_A[i]* k_basis[i-1]

            #print( "******************************")
            #print( v, E* k_basis[i], k_A[i]* k_basis[i-1])
            A =  np.linalg.norm( v_orto)
            v_orto = v_orto/ A

            #print(l, E, A)
            #print(v_orto)

            k_E.append(E)

            if A < 1e-14:
                break
            else :
                k_A.append(A)
                k_basis.append(v_orto)
                l += 1

        mat = np.zeros( (l,l), dtype = complex)

        for i in range(l-1):
            mat[i,i] = k_E[i]
            mat[i+1, i] = k_A[i+1]
            mat[i, i+1] = k_A[i+1]

        mat[l-1, l-1] = k_E[-1]

        return QWGraph(N = l, mat = mat)

    def eigen_basis( self, mode = "basis_plot"):

        if mode == "val":
            for j in range(len(self.eig_val)) :
                print(j, ") \t", self.eig_val[j])

        elif mode == "vec":
            for j in range(len(self.eig_vec)) :
                print(j, ") \t", self.eig_vec[j])

        elif mode == "plot_vec":

            x_range = np.arange(0,self.N)
            y_range = np.arange(0,len(self.eig_vec))

            modulus = np.zeros( ( self.N, len(self.eig_vec)) )
            phase = np.zeros( ( self.N, len(self.eig_vec)) )

            for j in range(len(self.eig_val)) :
                modulus[:,j] = np.abs( self.eig_vec[j])
                phase[:,j] = np.angle( self.eig_vec[j])

            fig, ax = plt.subplots( nrows = 1, ncols = 2, sharex = True, sharey = True)

            c1 = ax[0].pcolormesh(x_range, y_range, modulus, label = "mod")
            c2 = ax[1].pcolormesh(x_range, y_range, phase, label = "phase")

            fig.colorbar(c1, ax = ax[0])
            fig.colorbar(c2, ax = ax[1])

            ax[0].set_title("Projection modulus")
            ax[1].set_title("Relative phase")
            
            for stuff in ax :    
                stuff.set_xlabel('site')
                stuff.set_ylabel('eig_n')

            plt.show()

    def krylov_basis( self, start_state = None, mode = "x"):

        if start_state == None:
            start_state = self.get_start_state()

        #use Laczos algorithm to compute the new basis
        k_basis = []
        k_E = []
        k_A = []

        k_basis.append( start_state)
        k_A.append(0)
        l = 1
        for i in range(self.N):
            v = np.matmul( self.mat, k_basis[i])

            E = np.vdot(k_basis[i], v )

            v_orto = v - E* k_basis[i] - k_A[i]* k_basis[i-1]

            #print( "******************************")
            #print( v, E* k_basis[i], k_A[i]* k_basis[i-1])
            A =  np.linalg.norm( v_orto)
            v_orto = v_orto/ A

            #print(l, E, A)
            #print(v_orto)

            k_E.append(E)

            if A < 1e-8:
                break
            else :
                k_A.append(A)
                k_basis.append(v_orto)
                l += 1

        if mode == "e" :
            for e in k_E:
                print( e, "\t")
            return k_E

        if mode == "link" :
            for a in k_A:
                print( a, "\t")
            return k_A

        if mode == "basis_plot":

            modulus = np.zeros( (len(k_basis),self.N))
            phase = np.zeros( (len(k_basis),self.N))

            for j in range(len(k_basis)) :
                modulus[j,:] = np.abs( k_basis[j][:,0])
                phase[j,:] = np.angle( k_basis[j][:,0])

            x_range = np.arange(0,self.N)
            y_range = np.arange(0,len(k_basis))

            fig, ax = plt.subplots( nrows = 1, ncols = 2, sharex = True, sharey = True, figsize=(8, 4))

            c1 = ax[0].pcolormesh(x_range, y_range, modulus, label = "mod")
            c2 = ax[1].pcolormesh(x_range, y_range, phase, label = "phase")

            fig.colorbar(c1, ax = ax[0])
            fig.colorbar(c2, ax = ax[1])

            ax[0].set_title("Projection modulus")
            ax[1].set_title("Relative phase")
            
            for stuff in ax :    
                stuff.set_xlabel('site')
                stuff.set_ylabel('k_n')

            plt.show()

            return
        
        if mode == "link_plot":

            k_A = k_A[1:]
            
            fig, ax = plt.subplots()
            x = np.arange(1,len(k_A)+1)
            
            ax.scatter(x, k_A, label = "couplings")
            
            ax.set_xlabel('i')
            ax.set_ylabel('a')

            ax.legend()
            plt.show()

            return

        for i in range(len(k_basis)):
            vec = " "
            for elem in k_basis[i] :
                if mode == "x":
                    vec = vec + ("[]" if np.abs(elem)< 1e-10 else "[x]") + "\t"
                elif mode == "short_re":
                    vec = vec + "{: 0.2f}".format(np.abs(elem[0])) + "\t"
                else :
                    vec = vec + str(elem) + "\t"
            print( i, "|\t", vec)

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
    def from_igraph( ig, E = 0, ends = None) :

        out = QWGraph( ig.vcount())

        #todo assign name
        #out.code = ig["name"]
        
        out.mat = np.array ( ig.get_adjacency().data, dtype = complex)
        out.retrace_E(E)

        
        out.compute_re_coord()

        if not ends :
            out.start = 0
            out.target = out.N-1
        else :
            out.start , out.target = ends

        return out

    #return igraph object representing the QWGraph
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

    #return networkx object representing the QWGraph
    #the process does not tranfer phase info
    def to_networkx(self) :

        ref = np.zeros( (self.N, self.N))

        #format adjacency matix for igraph imput
        
        for i in range(self.N) :
            for m in range(self.N) :
                if self.mat[i][m] != 0 :
                    ref[i][m] = 1
            ref[i][i] = 0
        
        out = nx.from_numpy_matrix(ref)

        #print(ref)
        return out

    #construct a specific graph with the respective code number with variadic arguments

    def Graph_from_code(code , *args):

        if code == 0:
            return QWGraph.Line(args)
        if code == 1 :
            return QWGraph.Ring(args)
        if code == 2 :
            return QWGraph.Ring(args)

    def gfc(code, *args):
        return QWGraph.Graph_from_code(code, args)


    #Ring graph constructor
    def Ring(N, HANDLES = False, E = 0):

        code = "C"+ str(N)

        if N==0 :
            out = QWGraph(N)
            out.code = code
            return(out)

        outmat = np.identity(N, dtype=complex)* E

        for i in range(N):
            outmat[i][(i+1)%N] = complex(-1)
            outmat[(i+1)%N][i] = complex(-1)

        out = QWGraph(N, outmat, endpoints=(0, N//2))
        out.code = code


        #trace settings
        out.retrace_E(E)
        out.compute_re_coord()

        if HANDLES :
            return out.add_handles(1)

        return out

    #Line graph constructor
    def Line(N, E = 0, speedup = None):
        out = QWGraph(N)
        out.code = "P"+ str(N)

        if N==0 :
            return(out)

        out.retrace_E(E)

        for i in range(N-1):
            out.mat[i][(i+1)%N] = complex(-1)
            out.mat[(i+1)%N][i] = complex(-1)

            if speedup and i != 0 and i!= N-2 :
                out.mat[i][(i+1)%N] *= speedup
                out.mat[(i+1)%N][i] *= speedup


        #trace settings

        

        out.start = 0
        out.target = N-1
        

        return out

    #Multi path element graph constructor
    def Parallel(paths , p_len, E = 0):
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
        
        
        return out

    def SquareCut(E = 0):
        out = QWGraph.Ring(4, E = E)

        out = out.cut( (0,3))
        out.re_coord = [(1,3),(2,3)]
        out.code = "DiC4"

        return out


    #link distance between two given nodes
    def distance(self, start = None, to = None):
        if not start :
            start = self.start
        if not to :
            to = self.target
        
        graph = self.to_igraph()

        path = graph.get_shortest_paths(start,to)

        return len(path[0]) -1
        
    #get numpy basis vector
    def basis(self, i, qut = False):
        if qut:
            return qt.basis(self.N,i)
        
        out = np.zeros((self.N,1), dtype = complex)
        out[i]= 1

        return out

    #return localized start state for evolution
    def get_start_state(self, qut = False):
        return self.basis(self.start, qut)

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

    #get a visual rapresentation of the graph (relies on igraph)
    def plot(self) :

        ref = self.to_networkx()
        nx.draw_spring(ref, with_labels = True)

        # names = []
        # for i in range(self.N):
        #     names.append(str(i))

        # cols = []
        # for e in ref.es:
        #     if e.tuple in self.re_coord or e.tuple[::-1] in self.re_coord:
        #         cols.append("red")
        #     else :
        #         cols.append("black")

        # v_cols = ["yellow"]* self.N
        # v_cols[self.start] = "green"
        # v_cols[self.target] = "red"

        # ref.vs["label"] = names
        # ref.vs["color"] = v_cols
        # ref.es["color"] = cols
        # ig.plot( ref , layout = ig.Graph.layout_fruchterman_reingold(ref))


def get_list_x(gr_list, x_mode = "size"):
    out = []

    for i in range(len(gr_list)):
        if x_mode == "size":
            out.append(gr_list[i].N)
        elif x_mode == "dist":
            out.append( gr_list[i].distance())
        else :
            out.append(i)

    return out


############################

if __name__ == "__main__":
    a = QWGraph.Parallel(4,3)

    a.plot()
