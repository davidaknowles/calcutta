import time
import scipy.sparse as sp
import scipy.sparse.linalg as sp_linalg
import numpy as np
from sksparse.cholmod import cholesky, cholesky_AAt
import scipy.linalg

## note these only work for sparse invL_b
def get_diag_direct(invL_b):
    return (invL_b.T @ invL_b).diagonal()

def get_diag_fast(invL_b): 
    return (invL_b.multiply(invL_b)).sum(0) 

def get_diag_faster(invL_b):
    temp = invL_b.copy()
    temp.data *= temp.data
    return temp.sum(0)

def get_diag_all(invL_b): 
    v1 = get_diag_direct(invL_b)
    v2 = get_diag_fast(invL_b)
    v3 = get_diag_faster(invL_b)
    assert( np.abs(v1-v2).mean() < 1e-6 )
    assert( np.abs(v1-v3).mean() < 1e-6 )
    return np.array(v3)

def get_diag_invA(invA_b,b): 
    return (b.T @ invA_b).diagonal()

def get_diag_invA_fast(invA_b,b): 
    return b.multiply(invA_b).sum(0)

N = 5000
P = 1000
M = 8000

X = sp.rand( N, P, density = 0.005 )
#A = sp.csc_matrix(np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]]))
A = X.T @ X # P x P
A.nnz / np.prod(A.shape) # density of A
 
# b = sp.rand( 5, 1, density = 0.1 )
b = sp.rand( P, M, density = 0.01 ).tocsc()

# want b'inv(A)b = g'g: 6x6
# or really just the diagional 

eps = 1e-4
sp_eps = sp.diags(np.full(P,eps))


### Luke 
# Z = [ sp_linalg.lsqr(X.T, np.array(b.todense())[:,i])[0] for i in range(M) ]
# should be able to do 
# Z = solve(X',b)
# Zp = X @ solve(X, Z)
# (Z * Zp).sum(1)
# but don't really have the solvers I would need

### things that work

start_time = time.time()
ch = cholesky(A, beta=eps)
invL_b = ch.solve_L(b, use_LDLt_decomposition = False)
get_diag_faster(invL_b) # 0.36s
time.time() - start_time # 2.4s

start_time = time.time()
ch = cholesky(A, beta=eps)
invA_b = ch(b) # equivalent to solve(A,b)
get_diag_invA_fast(invA_b,b) # 0.04s
time.time() - start_time # 4.05s 

start_time = time.time()
ch_AAt = cholesky_AAt(X.T, beta = eps)
invL_b_AAt = ch_AAt.solve_L(b, use_LDLt_decomposition = False)
get_diag_faster(invL_b_AAt)
time.time() - start_time # 2.25s

start_time = time.time()
ch_AAt = cholesky_AAt(X.T, beta = eps)
invA_b_AAt = ch_AAt(b) # equivalent to solve(A,b)
get_diag_invA_fast(invA_b_AAt,b)
time.time() - start_time # 4.27s

start_time = time.time()
invA_b_spsolve = sp_linalg.spsolve(A + sp_eps, b) # docs say can't handle multiple b, but actually does. need small eps
get_diag_invA_fast(invA_b_spsolve,b)
time.time() - start_time # 16s


### Things that don't work: either answer is wrong or can't handle multiple b

invA_LU = sp_linalg.splu(A + sp_eps) # only works for square matrices, so can't do solve(X,b)
invA_b_LU = invA_LU.solve(b) # can't handle sparse b
get_diag_invA(invA_b_LU,b)


invX_b = sp_linalg.spsolve(X.T.tocsc() + sp_eps,b) # X must be square. Handles multiple b. 
#invX_b.dot(invX_b) # way off
(invX_b*invX_b).sum(0) # way off

invX_b = sp_linalg.lsqr(X, b)[0] # can't handle multiple b
#invX_b.dot(invX_b) # way off

invX_b = sp_linalg.cg(X + sp_eps, b)[0] # can't handle multiple b
invX_b.dot(invX_b) # way off

invA_b = sp_linalg.cg(A + sp_eps, b)[0] # can't handle multiple b
b.dot(invA_b) 

def cg_lo_solve(X, y, eps = 0.): 
    # Solve the equation (X'X) * b = y using cg

    def XtX_operator(x):
        return X.T.dot(X.dot(x)) + eps * x

    XtX_linear_operator = sp_linalg.LinearOperator((X.shape[1], X.shape[1]), matvec=XtX_operator)

    #return sp_linalg.cg(XtX_linear_operator, y)[0] # can't handle multiple b
    # sp_linalg.gmres(XtX_linear_operator, y)[0] # can't handle multiple b
    return sp_linalg.gmres(XtX_linear_operator, y)[0] # scipy.linalg.solve doesn't work with LinearOperator

invA_b = cg_lo_solve(X,b,eps=eps)
b.dot(invA_b)


### Bekas

# want diag(b'inv(A)b), which is MxM
S = 1000
# v = np.random.randn(M,S)
v = np.random.rand(M,S)*2. - 1.

# invA_b_v = sp_linalg.spsolve(A + sp_eps, b @ v) 
invA_b_v = sp_linalg.spsolve(A + sp_eps, b @ v) 

bT_invA_b_v = b.T @ invA_b_v

(v * bT_invA_b_v).sum(1) / (v * v).sum(1)  # seems pretty high variance
