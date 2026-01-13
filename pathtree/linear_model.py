import sys,os
TOOLSPATH='~/soft/src/dev/craft/gunfolds/tools/'
sys.path.append(os.path.expanduser(TOOLSPATH))

import ecj, bfutils
import bfutils as bfu
import graphkit as gk
import warnings
from statsmodels.tsa.api import VAR
from sympy.matrices import SparseMatrix
from progressbar import ProgressBar, Percentage, \
    Bar, RotatingMarker, ETA, FileTransferSpeed
from scipy import linalg, optimize
import numpy as np
import scipy
import ipdb

def symchol(M): # symbolic Cholesky
    B = SparseMatrix(M)
    t = B.row_structure_symbolic_cholesky()
    B = np.asarray(B)*0
    for i in range(B.shape[0]): B[i,t[i]] = 1
    return B

def G2SVAR(G):
    n = len(G)
    A,B = npG2SVAR(G)
    P,L,U = linalg.lu(B)
    A = linalg.inv(L).tolist()
    B = B.tolist()
    A = listplace(A, 0.0, 0.0)
    for i in range(0,n): A[i][i] = 1
    B = listplace(B, 0.0, 'e')
    for i in range(0,n): B[i][i] = 'e'
    return A,B,P

def G2AH(G):
    n = len(G)
    A,B = npG2SVAR(G)
    P,L,U = linalg.lu(B)
    A = linalg.inv(L).tolist()
    B = B.tolist()
    A = listplace(A, 0.0, 0.0)
    for i in range(0,n): A[i][i] = 1
    B = listplace(B, 0.0, 'e')
    for i in range(0,n): B[i][i] = 'e'
    return A,B,P

def bnf2CG(fname):
    d = eval(open(fname).read())
    G = {}
    for v in d:
        G[v] = {u: set((0,1)) for u in d[v]['pars']}
    G = ecj.tr(G)
    for v in G:
        ld = {u: set([(0,1)]) for u in G[v]}
        G[v] = ld
    return G

def npG2SVAR(G):
    n = len(G)
    A = bfu.graph2adj(G)
    B = np.tril(bfu.graph2badj(G))
    np.fill_diagonal(B,1)
    B = symchol(B)
    return A,B

def x2M(x, A, B, aidx, bidx):
    A[aidx] = x[:len(aidx[0])]
    B[bidx] = x[len(aidx[0]):]
    #B[(bidx[1],bidx[0])] = x[len(aidx[0]):]
    return A, B

def nllf(x, A, B, Y, aidx, bidx): # negative log likelihood
    A,B = x2M(x, A, B, aidx, bidx)
    T = Y.shape[1]
    X = Y[:,1:] - np.dot(A, Y[:,:-1])
    ldB = T*np.log(abs(1./linalg.det(B)))
    return ldB + 0.5*np.trace( np.dot(np.dot(B.T, B), np.dot(X,X.T)))

def nllf2(x, A, B, YY, XX, YX, T, aidx, bidx): # negative log likelihood
    A,B = x2M(x, A, B, aidx, bidx)
    AYX = np.dot(A, YX.T)
    S = YY - AYX - AYX.T + np.dot(np.dot(A,XX), A.T)
    ldB = T*np.log(abs(1./linalg.det(B)))
    return 0.5*np.dot(np.dot(B.T, B).T.flat,S.flat) + ldB
    #return ldB + 0.5*np.trace( np.dot(np.dot(B.T, B), S))

def VARbic(nllf, K, T):
    return 2*nllf + K*np.log(T)

def listplace(l, a, b):
    return [[b if x != a else x for x in row] if isinstance(row, list) else b if x != a else x for row in l]

# -------------------------------------------------------------------
# data generation
# -------------------------------------------------------------------

def randweights(n, c=0.1, factor=9):
    rw = np.random.randn(n)
    idx = np.where(abs(rw) < factor * c)
    if idx[0].size > 0:
        rw[idx] += np.sign(rw[idx]) * c * factor
    return rw

def transitionMatrix(cg, minstrength=0.1):
    A = gk.CG2adj(cg)
    edges = np.where(A==1)
    A[edges] = randweights(edges[0].shape[0], c=minstrength)
    l = linalg.eig(A)[0]
    c = 0
    pbar = ProgressBar(widgets=['Searching for weights: ', Percentage(), ' '], maxval=10000).start()
    while max(l*np.conj(l)) > 1:
        A[edges] = randweights(edges[0].shape[0], c=c)
        c += 1
        l = linalg.eig(A)[0]
        pbar.update(c)
    pbar.finish()
    return A

def sampleWeights(n, minstrength=0.1):
    r = scipy.randn(n)
    s = minstrength/np.min(np.abs(r))
    r = s*r
    return r

def transitionMatrix2(cg, minstrength=0.1):
    A = gk.CG2adj(cg)
    edges = np.where(A==1)
    A[edges] = sampleWeights(edges[0].shape[0], minstrength=minstrength)
    l = linalg.eig(A)[0]
    c = 0
    pbar = ProgressBar(widgets=['Searching for weights: ', Percentage(), ' '], maxval=10000).start()
    while max(l*scipy.conj(l)) > 1:
        A[edges] = sampleWeights(edges[0].shape[0], minstrength=minstrength)
        c += 1
        l = linalg.eig(A)[0]
        if c>pbar.maxval:
            raise ValueError
        pbar.update(c)
    pbar.finish()
    return A

def transitionMatrix3(cg, x0=None, minstrength=0.1):
    A = gk.CG2adj(cg)
    edges = np.where(A==1)

    try:
        s = x0.shape
        x = x0
    except AttributeError:
        A = initRandomMatrix(A, edges)
        x = A[edges]

    def objective(x):
        A[edges] = np.real(x)
        l = linalg.eig(A)[0]
        m = np.max(np.real(l*scipy.conj(l)))-0.99
        n = np.min(np.min(np.abs(x)),minstrength)-minstrength
        return m*m + 0.1*n*n

    o = np.zeros(len(edges))
    while np.min(np.abs(o[0])) < 0.8*minstrength:
        rpt = True
        while rpt:
            try:
                try:
                    o = optimize.fmin_bfgs(objective, x,
                                           gtol=1e-10, maxiter=100,
                                           disp=False, full_output=True)
                    A[edges]=np.real(o[0])
                    l = linalg.eig(A)[0]
                    if np.max(np.real(l*scipy.conj(l))) < 1:
                        rpt = False

                except:
                    rpt = True
            except Warning:
                x = scipy.randn(len(edges[0]))
                rpt = True
    A[edges]=np.real(o[0])
    return A

def initRandomMatrix(A, edges, maxtries=100, distribution='beta', stable=True):
    '''
    possible distributions:
    flat
    flatsigned
    beta
    normal
    uniform
    '''
    s = 2.0

    def init():
        if distribution=='flat':
            x = np.ones(len(edges[0]))
        elif distribution=='flatsigned':
            x = np.sign(scipy.randn(len(edges[0])))*scipy.ones(len(edges[0]))
        elif distribution=='beta':
            x = np.random.beta(0.5,0.5,len(edges[0]))*3-1.5
        elif distribution=='normal':
            x = scipy.randn(len(edges[0]))
        elif distribution=='uniform':
            x = np.sign(scipy.randn(len(edges[0])))*scipy.rand(len(edges[0]))
        else:
             raise ValueError('Wrong option!')
        return x

    def eigenvalue(A):
        l = linalg.eig(A)[0]
        s = np.max(np.real(l*scipy.conj(l)))
        return s

    x = init()
    A[edges] = x
    s = eigenvalue(A)
    alpha = np.random.rand()*(0.99-0.8)+0.8
    A = A/(alpha*s)
    s = eigenvalue(A)

    return A

def transitionMatrix4(g, minstrength=0.1, distribution='normal', maxtries=1000):
    A = gk.CG2adj(g)
    edges = np.where(A==1)
    s = 2.0
    c = 0
    pbar = ProgressBar(widgets=['Searching for weights: ',
                                Percentage(), ' '],
                       maxval=maxtries).start()
    while s > 1.0:
        minstrength -= 0.001
        A = initRandomMatrix(A, edges, distribution=distribution)
        x = A[edges]
        delta = minstrength/np.min(np.abs(x))
        A[edges] = delta*x
        l = linalg.eig(A)[0]
        s = np.max(np.real(l*scipy.conj(l)))
        c += 1
        if c > maxtries:
            return None
        pbar.update(c)
    pbar.finish()

    return A

def drawsamplesLG(A, nstd=0.1, samples=100):
    n = A.shape[0]
    data = scipy.zeros([n, samples])
    data[:,0] = nstd*scipy.random.randn(A.shape[0])
    for i in range(1,samples):
        data[:,i] = scipy.dot(A,data[:,i-1]) \
                    + nstd*scipy.random.randn(A.shape[0])
    return data


def getAgraph(n, mp=2, st=0.5, verbose=True):
    keeptrying = True
    while keeptrying:
        G = gk.rnd_cg(n, maxindegree=mp, force_connected=True)
        try:
            A = transitionMarix2(G, minstrength=st)
            keeptrying = False
        except ValueError as e:
            if verbose:
                print("!!! Unable to find strong links for a stable matrix !!!")
                print("*** trying a different graph")
    return {'graph':      G,
            'transition': A,
            'converges':  len(bfutils.call_undersamples(G))}

def getAring(n, density=0.1, st=0.5, verbose=True, dist='flatsigned', permute=False):
    keeptrying = True
    plusedges = bfutils.dens2edgenum(density,n)
    while keeptrying:
        G = gk.ringmore(n, plusedges, permute=permute)
        try:
            A = transitionMatrix4(G, minstrength=st, distribution=dist)
            try:
                s = A.shape
                keeptrying = False
            except AttributeError:
                keeptrying = True
        except ValueError:
            if verbose:
                print("!!! Unable to find strong links for a stable matrix !!!")
                print("*** trying a different graph")
    return {'graph':      G,
            'transition': A,
            'converges':  len(bfutils.call_undersamples(G))}



# -------------------------------------------------------------------
# estimation
# -------------------------------------------------------------------

def scoreAGraph(G, data, x0 = None):
    A,B = npG2SVAR(G)
    n = len(G)
    a_idx = np.where(A != 0)
    b_idx = np.where(B != 0)
    K = np.sum(len(a_idx[0])+len(b_idx[0])/2)

    if x0:
        o = optimize.fmin_bfgs(nllf, x0, args=(A, B, data, a_idx, b_idx),
                               disp=False, full_output=True)
    else:
        o = optimize.fmin_bfgs(nllf, scipy.randn(len(a_idx[0])+len(b_idx[0])),
                               args=(np.double(A), np.double(B),
                                     data, a_idx, b_idx),
                               disp=False, full_output=True)
    ipdb.set_trace()
    return 2*o[1] + K*np.log(data.shape[1]) #VARbic(o[1],K,data.shape[1])

def scoreAGraph2(G, data, x0 = None):
    A,B = npG2SVAR(G)
    n = len(G)
    a_idx = np.where(A != 0)
    b_idx = np.where(B != 0)

    T = data.shape[1]
    YY = np.dot(data[:,1:],data[:,1:].T)
    XX = np.dot(data[:,:-1],data[:,:-1].T)
    YX = np.dot(data[:,1:],data[:,:-1].T)

    K = scnpipy.sum(len(a_idx[0])+len(b_idx[0])/2)

    if x0:
        o = optimize.fmin_bfgs(nllf2, x0,
                               args=(np.double(A), np.double(B),
                                     YY,XX,YX,T,a_idx, b_idx),
                                     gtol=1e-12, maxiter=500,
                                     disp=False, full_output=True)
    else:
        o = optimize.fmin_bfgs(nllf2, scipy.randn(len(a_idx[0])+len(b_idx[0])),
                               args=(np.double(A), np.double(B),
                                     YY,XX,YX,T,a_idx, b_idx),
                                gtol=1e-12, maxiter=500,
                                disp=False, full_output=True)

    return 2*o[1] + K*np.log(data.shape[1]) #VARbic(o[1],K,data.shape[1])

def estimateG(G,YY,XX,YX,T,x0=None):
    A,B = npG2SVAR(G)
    K = np.sum(abs(A)+abs(B))
    a_idx = np.where(A != 0)
    b_idx = np.where(B != 0)
    try:
        s = x0.shape
        x = x0
    except AttributeError:
        x = scipy.randn(K)
    o = optimize.fmin_bfgs(nllf2, x,
                           args=(np.double(A), np.double(B),
                                 YY,XX,YX,T,a_idx, b_idx),
                           disp=False, full_output=True)
    A,B = x2M(o[0], np.double(A), np.double(B), a_idx, b_idx)
    return  A,B


def data2AB(data,x0=None):
    n = data.shape[0]
    T = data.shape[1]
    YY = np.dot(data[:,1:],data[:,1:].T)
    XX = np.dot(data[:,:-1],data[:,:-1].T)
    YX = np.dot(data[:,1:],data[:,:-1].T)

    model = VAR(data.T)
    r = model.fit(1)
    A = r.coefs[0,:,:]

    #A = np.ones((n,n))
    B = np.ones((n,n))
    np.fill_diagonal(B,0)
    B[np.triu_indices(n)] = 0
    K = np.int(scipy.sum(abs(B)))#abs(A)+abs(B)))

    a_idx = np.where(A != 0)
    b_idx = np.where(B != 0)
    np.fill_diagonal(B,1)

    try:
        s = x0.shape
        x = x0
    except AttributeError:
        x = np.r_[A.flatten(),0.1*scipy.randn(K)]
    o = optimize.fmin_bfgs(nllf2, x,
                           args=(np.double(A), np.double(B),
                                 YY,XX,YX,T,a_idx, b_idx),
                           gtol=1e-12, maxiter=500,
                           disp=False, full_output=True)
    ipdb.set_trace()
    A,B = x2M(o[0], np.double(A), np.double(B), a_idx, b_idx)
    B = B+B.T
    return  A,B

def amap(f, a):
     v = np.vectorize(f)
     return v(a)

def AB2intAB(A,B, th=0.09):
    A[amap(lambda x: abs(x) > th, A)] = 1
    A[amap(lambda x: abs(x) < 1, A)] = 0
    B[amap(lambda x: abs(x) > th, B)] = 1
    B[amap(lambda x: np.abs(x) < 1, B)] = 0
    return A,B

def intAB2graph(A,B):
    n = A.shape[0]
    g = {str(i):{} for i in range(1,n+1)}

    for i in range(n):
        for j in range(n):
            if A[j,i]: g[str(i+1)][str(j+1)] = set([(0,1)])

    for i in range(n):
        for j in range(n):
            if B[j,i] and j!=i:
                if str(j+1) in g[str(i+1)]:
                    g[str(i+1)][str(j+1)].add((2,0))
                else:
                    g[str(i+1)][str(j+1)] = set([(2,0)])
    return g

def data2graph(data,x0=None, th=0.09):
    A,B = data2AB(data,x0=x0)
    Ab,Bb = AB2intAB(A,B,th=th)
    return intAB2graph(Ab,Bb)

def data2VARgraph(data, pval=0.05):
    model = VAR(data.T)
    r = model.fit(1)
    A = r.coefs[0,:,:]
    n = A.shape[0]
    g = {str(i):{} for i in range(1,n+1)}

    for i in range(n):
        for j in range(n):
            if np.abs(A[j,i]) > pval: g[str(i+1)][str(j+1)] = set([(0,1)])
    return g

def data2VARgraph_model(data, pval=0.05):
    model = VAR(data.T)
    r = model.fit(1)
    A = r.coefs[0,:,:]
    n = A.shape[0]
    g = {str(i):{} for i in range(1,n+1)}

    for i in range(n):
        for j in range(n):
            if np.abs(A[j,i]) > pval: g[str(i+1)][str(j+1)] = set([(0,1)])
    return g, r
