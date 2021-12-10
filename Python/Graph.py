import numpy as np
from numpy.linalg import eigh

#get numpy basis vector

def basis(N,l):
    out = np.zeros(N, dtype = complex)
    out[l]= 1

    return out

#QW simulation oriented graph class

class Graph(object) :

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
        

    def retrace_E(self, E):
        for i in range(self.N):
            self.mat[i][i] = E

        self.update_eigen()

    """def retrace_conn(self):
        for i in range(N):
            count = 0
            
            self.mat[i][i] = E"""





#ring graph constructor

def Ring(N, E = 2):
    out = Graph(N)
    out.code = "C"+ str(N)

    if N==0 :
        return(out)

    out.retrace_E(E)

    for i in range(N):
        out.mat[i][(i+1)%N] = complex(-1)
        out.mat[(i+1)%N][i] = complex(-1)

    #trace settings

    out.update_eigen()

    out.start = 0
    out.target = int(N/2)
    
    if N!= 1 :
        out.re_coord.append( (0,N-1))

    return out

#Line graph constructor

def Line(N, E = 2):
    out = Graph(N)
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

