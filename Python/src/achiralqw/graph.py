from typing import List
import numpy as np
from numpy.linalg import eigh
import igraph as ig
import qutip as qt
import networkx as nx

#QW simulation oriented graph class

class QWGraph(object) :
    """
    Achiral QW simulation oriented graph class
    """

    def __init__(self, N , mat = None, endpoints = None) :

        if mat is None :
            self.code = "e"
        else :
            self.code = "m"

        self.N = N
        self._init_mat(mat)
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

    def _init_mat(self, mat):
        """
        Initialize laplacian matrix, possibly creating a new one
        This method must be overwritten if the internal representation of the matrix changes from numpy ndarray
        """

        if mat is None :
            self.mat = np.zeros((self.N, self.N), dtype = complex)
        else :
            #mat is assumed to be a square numpy ndarray
            self.N = mat.shape[0]
            self.mat = mat

    def update_eigen(self):
        """
        Recompute and store eigenvalues and eigenvectors with numpy.linalg routine
        WARNING: this may be time consuming for large graphs!
        """
        self.eig_val, self.eig_vec = eigh(self.mat)

        self.eig_val = np.reshape(self.eig_val, (self.N,1))
        
    def retrace(self, T : List[float]):
        """
        Set new values to trace elements of Laplacian matrix
        """

        for i in range(self.N):
            self.mat[i][i] = T[i]

    def retrace_E(self, E : float):
        """
        Set new value to trace elements of Laplacian matrix
        """

        for i in range(self.N):
            self.mat[i][i] = E

        

    def retrace_conn(self):
        """
        Set trace elements with the outdegree of the corrisponding node
        """
        ref = self.to_igraph()
        
        for i in range(self.N):
            self.mat[i][i] = ref.degree(i)
            print(ref.degree(i))

    def _buffer_trace(self):
        """
        utiliry method to save trace elements in a buffer vector
        (useful for igraph objcet)
        """
        out = [ self.mat[i][i] for i in range(self.N)]

        return out

    def rephase(self, phi_vec = [0]) :
        """
        Assign new phase value to the phased cycle links
        Note: new phases values are inerpreted as being in radians!!!
        The operation does not change the coupling magnitude

        Args:
        phi_vec -> vector of phases, has to be of the same length as re_coord

        """
        if isinstance(phi_vec, (list, np.ndarray)):
            if( len( self.re_coord) != len(phi_vec)):
                raise ValueError("rephase() error: wrong number of phases given")
        else :
            # phi_vec is a single float: format to a numpy vector
            temp = phi_vec
            phi_vec = np.empty(1)
            phi_vec[0] = temp

        exp_vec = np.exp(1j * phi_vec)

        for i, edge in enumerate(self.re_coord) :
            self.mat[edge[0]][edge[1]] = -1* np.abs(self.mat[edge[0]][edge[1]])* exp_vec[i]
            self.mat[edge[1]][edge[0]] = np.conjugate(self.mat[edge[0]][edge[1]])

    def get_phases(self):
        """
        Return the vector of phases in radians of the phased links
        """

        out = np.empty(len(self.re_coord))

        for i,edge in enumerate(self.re_coord) :
            out[i] = np.angle(self.mat[edge[0]][edge[1]])

        return out

    def recouple(self, pos, value):
        """
        assign new value to an existing coupling
        """

        if self.mat[pos[0]][pos[1]] == 0:
            raise ValueError("Trying to recouple a non existing edge")

        self.mat[pos[0]][pos[1]] = value
        self.mat[pos[0]][pos[1]] = np.conjugate(value)

    def cut(self, cut_vec):
        """
        Create a new edge, rephasing links are going to be recomputed

        ##Arguments
        cut_vec : list[(int,int)] -> list of the new edges to insert in the graph
        """

        #print(self.code, cut_vec)

        if type(cut_vec) is  tuple:
            cut_vec = [cut_vec]

        ref = self.to_igraph()
        ref.add_edges(cut_vec)

        temp = self._buffer_trace()

        out = QWGraphBuilder.fromIgraph(ref, ends = (self.start, self.target))
        out.retrace(temp)

        return out

    def reverse(self ):
        """
        Exchange start and target site
        """
        self.start, self.target = self.target, self.start


    def speedup(self, su):
        """
        Multiply all non trace elements by a constant (speedup)
        """
        temp = self._buffer_trace()

        self.mat = self.mat * su
        self.retrace(temp)
    
    def join_link(self, other):
        """
        concatenate two graph creating the link between first end site and second start site
        (+ operator)
        """

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

    def join_nolink(self, other):
        """
        concatenate two graph merging the first end site and the second start site
        ( | [or] operator)
        """

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

        #account for the possible swap while copying matrices
        out.start  = self.start                 if self.start != self.N-1   else self.target
        out.target = self.N + other.target-1    if other.target != 0        else self.N + other.start -1

        return out

    def __or__(self, other) :
        return QWGraph.join_nolink(self, other)

    def chain(self, rep, space = 0, speedup = 1, HANDLES = True):
        """        
        concatenate multiple graph units into a chain
        ( * operator)
        speedup affects the couplings before adding the handles

        todo: fix the fact that intra unit space does not get speedup
        """

        #create self copy with speedup

        su_mat = self.mat * speedup
        unit = QWGraphBuilder.Line(space, speedup = speedup) + QWGraph(self.N,mat = su_mat, endpoints = (self.start,self.target))
        
        #todo: df am i doing with P graph trace?
        temp = self._buffer_trace()
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
        """
        Add external links and nodes on transport ends("handles")

        mode list:

        both -> add "size" lenght handles at both ends
        l -> only start site
        r -> only end site
        fixl -> set lenght on start site, "size" lenght on end site
        fixe -> set lenght on end site, "size" lenght on start site
        """
        
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

        out = QWGraphBuilder.Line(sx) + self + QWGraphBuilder.Line(dx)
        out.code = "h(" + temp_code + ")"

        return out

    def eigen_basis( self, mode = "vec"):
        """
        Pretty printing method with information on eigenvector basis

        List of modes

        val -> eigenvalues
        vec -> eigenvectors
        """

        if mode == "val":
            for j in range(len(self.eig_val)) :
                print(j, ") \t", self.eig_val[j])

        elif mode == "vec":
            for j in range(len(self.eig_vec)) :
                print(j, ") \t", self.eig_vec[j])

    def krylov_basis( self, start_state = None, mode = ""):
        """
        Compute krylov basis relative to the start state of the graph and eventually print some details

        List of modes:

        e           -> print krilov basis energies
        link        -> print krilov basis couplings
        x           -> represent krylov basis in the site basis as occipied sites []/[x]
        short_re    -> represent krylov basis in site basis with short formatted real part of projection
        default     -> simple pring of krylov basis in the site basis

        returns tuple(3):
        k_basis -> krylov basis vectors in site basis
        k_E     -> Krylov energies
        k_A     -> Krylov couplings

        """

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

        elif mode == "link" :
            for a in k_A:
                print( a, "\t")

        elif mode != "" :
            for i in range(len(k_basis)):
                vec = " "
                for elem in k_basis[i] :
                    if mode == "x":
                        vec = vec + ("[]" if np.abs(elem)< 1e-10 else "[x]") + "\t"
                    elif mode == "short_re":
                        vec = vec + "{: 0.2f}".format(np.abs(elem[0])) + "\t"
                    elif mode == "default":
                        vec = vec + str(elem) + "\t"
                print( i, "|\t", vec)

        return k_basis, k_E, k_A

    def compute_re_coord(self) :
        """
        Update rephasing link vector as missing link from a spanning tree of the graph
        """

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

    def krylov_transform( self, start_state = None):
        """
        Compute the krilov base and construct a new graph in the krylov subspace relative to a given start state
        """

        k_basis, k_E, k_A = self.krylov_basis(start_state)

        l = len(k_basis)

        mat = np.zeros( (l,l), dtype = complex)

        for i in range(l-1):
            mat[i,i] = k_E[i]
            mat[i+1, i] = k_A[i+1]
            mat[i, i+1] = k_A[i+1]

        mat[l-1, l-1] = k_E[-1]

        return QWGraph(N = l, mat = mat)   

    def to_igraph(self) :
        """
        Build Igraph object representing the QWGraph
        the process does not transfer phase info
        """

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

    def to_networkx(self):
        """
        Build networkx object representin the QWGtaph.
        The process does not transfer phase info
        """

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

    def distance(self, start = None, to = None):
        """
        Return distance in links betweeen two given nodes
        (Actually a wrapper of igraph get_shortest_paths)
        """
        if not start :
            start = self.start
        if not to :
            to = self.target
        
        graph = self.to_igraph()

        path = graph.get_shortest_paths(start,to)

        return len(path[0]) -1
        
    def basis(self, i, qut = False):
        """
        Return chosen basis vector as nupmy array
        """
        if qut:
            return qt.basis(self.N,i)
        
        out = np.zeros((self.N,1), dtype = complex)
        out[i]= 1

        return out

    def get_start_state(self, qut = False):
        """
        Get numpy vector representation of the localized start state
        """
        return self.basis(self.start, qut)

    def get_projector(self, i = None):
        """
        get QuTip projector operator on the Nth site
        """
        if not i:
            i = self.target
            
        out_mat = np.zeros((self.N,self.N))
        out_mat[i,i] = 1

        out = qt.Qobj(out_mat)

        return out

    def get_h(self):
        """
        get the Hamiltonian relative to the graph as QuTip object (for QuTip simulator)
        """
        return qt.Qobj(self.mat)

    def get_phase_n(self):
        """
        Get the number of free phases for the graph
        """
        return len(self.re_coord)

class SparseQWGraph( QWGraph):
    """
    Implementation of QWGraph with and underlying matrix representation
    """

    def _init_mat(self, mat):
        """
        Initialize laplacian matrix, possibly creating a new one.
        The matrix representation is a complex scipy sparse dictionary of key
        """

        if mat is None :
            self.mat = sparse.dok_matrix((self.N, self.N), dtype = complex)
        else :
            #mat is assumed to be a square scipy sparse dok matrix
            # we need efficient slicing and access to the elements
            assert isinstance(mat, sparse.dok_matrix)
            self.N = mat.shape[0]
            self.mat = mat

    def to_adjacency(mat):
        
        out = sparse.dok_matrix(mat.shape, dtype = "int")
        
        for k in mat.keys():
            if k[0] != k[1]:
                out[k] = 1
                
        return out

    def update_eigen(self):
        """
        Recompute and store eigenvalues and eigenvectors with numpy.linalg routine
        WARNING: this may be time consuming for large graphs!

        Apparently there is no sparse.linalg routine which retrieve all the eigenvectors
        Therefore we must converte to dense matrix an use te usual linalg.eigh
        """
        self.eig_val, self.eig_vec = (self.mat.todense())

        self.eig_val = np.reshape(self.eig_val, (self.N,1))

    def compute_re_coord(self) :
        """
        Update rephasing link vector as missing link from a spanning tree of the graph
        """

        mst = sparse_mst(self.mat)

        self.re_coord = []

        for edge in self.mat.keys():
            if edge[0] > edge[1] and my_mst[edge] == 0 :
                self.re_coord.append(edge)

##        names = []
##        for i in range(self.N):
##            names.append(str(i))
##
##        tree.vs["label"] = names
##        ref.vs["label"] = names
##        ig.plot(ref -tree)

        self.re_coord = n_re_coord

    def distance(self, start = None, to = None):
        """
        Return distance in links betweeen two given nodes
        (Actually a wrapper of igraph get_shortest_paths)
        """
        if not start :
            start = self.start
        if not to :
            to = self.target
        
        dist = sparse.csgraph.shortest_path(self.mat, method = "D", directed = False, unweighted= True, indices = start)

        return dist[target]

class QWGraphBuilder(object):

    @staticmethod
    def graphFromCode(code : str, *args) -> QWGraph:
        """
        Construct a specific graph with the respective code number with variadic arguments
        """

        if code == 0:
            return QWGraphBuilder.Line(args)
        if code == 1 :
            return QWGraphBuilder.Ring(args)
        if code == 2 :
            return QWGraphBuilder.Ring(args)

    @staticmethod
    def gfc(code : str, *args) -> QWGraph:
        """
        wrapper for GraphFromCode
        """
        return QWGraphBuilder.Graph_from_code(code, args)

    @staticmethod
    def fromIgraph( ig, E = 0, ends = None) :
        """
        return QWGraph instance from a igraph object
        note: the resulting Laplacian matrix is of course going to be real valued
        """

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

    #todo : finish this thing
    @staticmethod
    def fromNetworkx( nx_graph, E = 0, ends = None):
        raise NotImplementedError("to be implemented")

    @staticmethod
    def Ring(N : int, HANDLES : bool = False, E : float  = 0, COMPUTE_EIGEN = False, sparse = False) -> QWGraph:
        """
        Ring graph constructor
        """

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

        if COMPUTE_EIGEN:
            out.update_eigen()

        return out

    @staticmethod
    def Line(N : int , E : float = 0, speedup : float = None, COMPUTE_EIGEN = False) -> QWGraph:
        """
        Line graph constructor
        """
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

        if COMPUTE_EIGEN:
            out.update_eigen()

        return out

    @staticmethod
    def Parallel(paths : int , p_len : int , E : float = 0, COMPUTE_EIGEN = False) -> QWGraph:
        """
        Multi path element graph constructor
        """
        ref = ig.Graph()

        N = 2 + paths* (p_len-1)

        ref.add_vertices(N)

        for i in range(paths):
            ref.add_edge(0, i+1)
            ref.add_edge(N-2-i, N-1)

        for i in range(p_len-2):
            for m in range(paths):
                ref.add_edge(1 + paths*i + m,1 + paths*(i+1) + m)

        out = QWGraphBuilder.fromIgraph(ref)
        out.retrace_E(E)
        
        if COMPUTE_EIGEN:
            out.update_eigen()

        return out

    @staticmethod
    def SquareCut(E  : float = 0, COMPUTE_EIGEN = False) -> QWGraph:
        """
        C4 with extra start to and cut and straight rephasal links
        """
        out = QWGraphBuilder.Ring(4, E = E)

        out = out.cut( (0,3))
        out.re_coord = [(1,3),(2,3)]
        out.code = "DiC4"

        if COMPUTE_EIGEN:
            out.update_eigen()

        return out

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


def sparse_mst( smat ):
    #do it yourself implementation of kruskal algorithm for a 01 adjacency matrix
    #suited for dok matrices
    
    #union find parent list
    par = np.arange(0,smat.shape[0])
    
    out = sparse.dok_matrix(smat.shape, dtype = "int")
    
    #primitives for union-find
    def same(a,b) ->bool :
        
        p1 = par[a]
        while p1 != par[p1] :
            p1 = par[p1]
            par[a] = p1
        
        p2 = par[b]
        while p2 != par[p2] :
            p2 = par[p2]
            par[b] = p2
            
        #print( "same {} {}\t".format(a,b), par)
        return p1 == p2
        
    def join(a,b) -> None :
        
        p1 = par[a]
        while p1 != par[p1] :
            p1 = par[p1]
            par[a] = p1
        
        p2 = par[b]
        while p2 != par[p2] :
            p2 = par[p2]
            par[b] = p2
            
        par[p2] = p1
        
    for edge in smat.keys():
        if not same(edge[0], edge[1]):
            join(edge[0], edge[1])
            out[edge] = 1
            out[edge[1],edge[0]] = 1
            
    return out
