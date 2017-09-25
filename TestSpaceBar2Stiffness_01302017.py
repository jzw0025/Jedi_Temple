# -*- coding: utf-8 -*-
import numpy as np
from scipy.sparse import *
from scipy.sparse.linalg import *

NodeCoordinates = np.array([[0.0,0.0,0.0],[10.0,0.0,0.0],[10.0,10.0,0.0]])
ElementNodes = [[0,1],[1,2],[0,2]]

def SpaceBar2Stiffness(coordinator, Em, A):
    """
    this calculate local space bar element.
    input coordinator: [[x1,y1,z1],[x2,y2,z2]]
    """
    x1, y1, z1 = coordinator[0][0], coordinator[0][1], coordinator[0][2]
    x2, y2, z2 = coordinator[1][0], coordinator[1][1], coordinator[1][2]
    x21, y21, z21 = x2-x1, y2-y1, z2-z1

    EA = Em*A
    LL = x21**2 + y21**2 + z21**2
    L = np.sqrt(LL)
    LLL = LL*L
    Ke = Em*A/LLL*np.array([[[x21*x21],[x21*y21],[x21*z21],[-x21*x21],[-x21*y21],[-x21*z21]],
                                             [[y21*x21],[y21*y21],[y21*z21],[-y21*x21],[-y21*y21],[-y21*z21]],
                                             [[z21*x21],[z21*y21],[z21*z21],[-z21*x21],[-z21*y21],[-z21*z21]],
                                             [[-x21*x21],[-x21*y21],[-x21*z21],[x21*x21],[x21*y21],[x21*z21]],
                                             [[-y21*x21],[-y21*y21],[-y21*z21],[y21*x21],[y21*y21],[y21*z21]],
                                             [[-z21*x21],[-z21*y21],[-z21*z21],[z21*x21],[z21*y21],[z21*z21]]])

    # return local stiffness matrix
    return Ke[:,:,0]

def SpaceTrussMasterStiffness(Node, Elem, EleMat, EleFab):
    NumElem = len(Elem)
    NumNode = len(Node)
    #K = np.zeros((3*NumNode,3*NumNode))
    K = csc_matrix((3*NumNode,3*NumNode))

    # element assembly
    for k in range(NumElem):
        ni = Elem[k][0]
        nj = Elem[k][1]
        print ni,nj
        # mapping the element local node into a global node
        # the python system maps the coordinator from 0 to (n-1)
        # the node should start from concurrent 3*i+(0,1,2)
        eftab = [3*ni, 3*ni+1, 3*ni+2, 3*nj, 3*nj+1, 3*nj+2]
        print eftab
        ncoor = np.array([Node[ni], Node[nj]])
        Em = EleMat[k]
        A = EleFab[k]
        #print ncoor
        Ke = SpaceBar2Stiffness(ncoor, Em, A)
        print Ke
        #print Ke
        neldof = len(Ke)
        for i in range(neldof):
            for j in range(neldof):
                ii = eftab[i]
                jj = eftab[j]
                #print ii,jj
                # K[ii][jj] = K[ii][jj] + Ke[i][j]
                K[ii,jj] = K[ii,jj] + Ke[i][j] # adapted for sparse matrix
    return K

def PrescDisplacementDOFTags(nodtag):
    numnod = len(nodtag)
    pdof = []
    k = 0
    for n in range(numnod):
        m = len(nodtag[n])
        for j in range(m):
            if nodtag[n][j] > 0:
                pdof.append(k+j)
        k+=m
    return pdof

def PrescDisplacementDOFValues(nodtag, nodval):
    numnod = len(nodtag)
    pval = []
    k = 0
    for n in range(numnod):
        m = len(nodtag[n])
        for j in range(m):
            if nodtag[n][j] > 0:
                pval.append(nodval[n][j])
        k+=m
    return pval

def ModifiedMasterStiffness(nodtag, K):
    #n = len(K)
    print K.shape
    n = K.shape[0]
    Kmod = K
    pdof = PrescDisplacementDOFTags(nodtag)
    np = len(pdof)
    for k in range(np):
        i = pdof[k]
        for j in range(n):
            Kmod[i,j] = Kmod[j,i] = 0
            Kmod[i,i] = 1
    return Kmod

def ModifiedNodeForces(nodtag, nodval, K, f):
    # the Stiffness matrix must be the original master stiffness without any modifications
    #n = len(K)
    print K.shape
    n = K.shape[0]
    pdof = PrescDisplacementDOFTags(nodtag)
    np = len(pdof)
    pval = PrescDisplacementDOFValues(nodtag,nodval)
    c = [1]*n
    fmod = f
    for k in range(np):
        i = pdof[k]
        c[i] = 0

    for k in range(np):
        i = pdof[k]
        d = pval[k]
        fmod[i] = d
        if d == 0:
            continue
        for j in range(n):
            print fmod[j]
            fmod[j] = fmod[j] - K[i,j]*c[j]*d
            print K[i,j], c[j], d, fmod[j]

    return fmod

if __name__ == "__main__":
    print "This is a test."
    import scipy.io as sio
    # this is a test case
    # 40. 60. 120. −40. −60. −120.
    # 60. 90. 180. −60. −90. −180.
    # 120. 180. 360. −120. −180. −360.
    # −40. −60. −120. 40. 60. 120.
    # −60. −90. −180. 60. 90. 180.
    # −120. −180. −360. 120. 180. 360.
    #K= SpaceBar2Stiffness(np.array([[0,0,0],[2,3,6]]), 343, 10)
    K= SpaceBar2Stiffness(np.array([[10,0.0,0.0],[10.0,10.0,0.0]]), 100, 0.5)

    EleMat = [1.0,0.5,2*np.sqrt(2)]
    EleFab = [100,100,100]
    Km =  SpaceTrussMasterStiffness(NodeCoordinates, ElementNodes, EleMat, EleFab)

    nodtag = [[1,1,1], [0,1,1], [0,0,1]]
    nodval = [[0,0,0], [0,0,0], [2,1,0]]
    #nodtag1 = [[1,1],[0,1],[0,0]]
    #nodval1 = [[1,2],[0,3],[0,0]]
    #print PrescDisplacementDOFTags(nodtag1)
    #print PrescDisplacementDOFValues(nodtag1,nodval1)
    #Km = np.ones((6,6))
    Km2 = Km
    print Km2.shape
    Kmod=ModifiedMasterStiffness(nodtag,Km2)
    print Kmod
    #print Km==Kmod
    fm=np.array([0,0,0,0,0,0,2,1,0.0])
    fm2 = np.copy(fm)

    fmod=ModifiedNodeForces(nodtag,nodval,Km,fm2);
    #print fmod

    #print np.linalg.solve(Kmod,fmod)
    print spsolve(Kmod,fmod)

    sio.savemat("/Users/junchaowei/Desktop/nodes.mat",{"par1":NodeCoordinates.T})
    sio.savemat("/Users/junchaowei/Desktop/truss.mat",{"par1":ElementNodes})
    sio.savemat("/Users/junchaowei/Desktop/values.mat",{"par1":fmod})
