################################################################################
'''Module "helpers.py", July 2010, is part of
BETR-Research by Harald von Waldow <hvwaldow@chem.ethz.ch>, which is
based on BETR-Global by Matt MacLeod <matthew.macleod@chem.ethz.ch>

This module contains various helper functions'''
################################################################################
from numpy import *
from globalz import *
import scipy.sparse as sp
import sys

def sparse_slice(matrix,columns,rows):
    ''' slice a generic sparse matrix''' 
    spformat=A.getformat()
    if max(columns) > shape(martix)[0] or max(rows) > shape(martix)[1]:
        raise IndexError
    row_mask= sp.coo_matrix((ones(len(rows)),(rows, arange(0,len(rows)))), shape = (dim, len(rows))).T
    col_mask= sp.coo_matrix((ones(len(columns)),(columns, arange(0,len(columns)))), shape = (dim, len(columns)))
    
    return (row_mask * matrix * col_mask).asformat(spformat)

def findemptycells(A):
    """ returns a list of indices with zero rows resp. columns"""
    occupied=A.tocsr().indices
    idx=[i for i in range(shape(A)[0]) if i not in occupied]
    return(idx)

#    B=A.todense()
#    for i in range(0,B.shape[0]):
#        if (all(B[i,:] == 0)):
#            idx.append(i)
#            if (not all(B[:,i] == 0)):
#                sys.exit("findemptycells: empty cell receiving some ! Abort !")
#    return(idx)

def findnonemptycells(A):      # function added by HW, 17.01.2011
    """ returns a list of indices with zero rows resp. columns"""
    idx=[]
    for i in range(0,A.shape[0]):
        if (any(A[i,:] != 0)):
            idx.append(i)
    return(idx)

def compressmat(A, killidx):
    ''' deletes rows and columns in killidx from matrix'''
    spformat=A.getformat()
    dim= shape(A)[0]
    keepidx =[i for i in range(shape(A)[0]) if not i in killidx]
    X = sp.coo_matrix((ones(len(keepidx)),(keepidx, arange(0,len(keepidx)))), shape = (dim, len(keepidx)))
    return (X.T * A * X).asformat(spformat)

#def compressmat(A, killidx):
#    ''' deletes rows and columns in killidx from matrix'''
#    spformat=A.getformat()
#    B=sp.coo_matrix(delete(delete(A.todense(),killidx,0), killidx,1))
#    return B.asformat(spformat)
    
    
def compressvec(v, killidx):
    """ Takes vector v and returns weeded vector"""
    vw = delete(v, killidx)
    return vw

def expandmatrix(A,killidx):
    '''re-expands previosly compressed matrix'''
    spformat=A.getformat()
    dim= shape(A)[0]+len(killidx)
    keepidx =[i for i in range(shape(A)[0]) if not i in killidx]
    X = sp.coo_matrix((ones(len(keepidx)),(keepidx, arange(0,len(keepidx)))), shape = (dim, len(keepidx)))
    return (X * A * X.T).asformat(spformat) 
    

#def expandmatrix(A, killidx):
#    spformat=A.getformat()
#    '''re-expands previosly compressed matrix'''
#    killidx=sort(killidx)
#    B=A.toarray()
#    for i in killidx:
#        B = insert(B,i,[0],axis=0)
#        B = insert(B,i,[0],axis=1)
#    B = sp.coo_matrix(matrix(B))
#    return(B.asformat(spformat))
    
    
    

def expandvec(v, killidx):
    """ Re - inserts zeros at the killidx - positions"""
    ve=copy(v)
    killidx=sort(killidx)
    for i in killidx:
        ve=insert(ve,i,0)
    return ve

def cell2regcomp(cell,m):
    '''  translates cell index into region-number and compartment ID '''
    reg=floor(cell/m.nocomp) + 1
    comp=mod(cell,m.nocomp) + 1
    return (reg.astype(int),comp.astype(int))

def tocell(reg,comp,m):
    ''' translates region number and compartment ID into cell-index'''
    cell=m.nocomp*(reg-1)+comp-1 
    return cell

###############################################################################
