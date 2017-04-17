################################################################################
'''Module "mkbigD.py", July 2010, is part of
BETR-Research by Harald von Waldow <hvwaldow@chem.ethz.ch>, which is
based on BETR-Global by Matt MacLeod <matthew.macleod@chem.ethz.ch>

The function mkflowmat(m) takes a model object (m) as argument and returns
a list of sparse matrices which contain D-values from inter-cell flows
for each timestep.

The function mkprocmat(m) takes a model object (m) as argument and returns
a list of sparse matrices which contain D-values from intra-cell flows
for each timestep.

The function mkbigD(m) is called from BETRS.py, calls mkflowmat and mkprocmat
and returns the final list of system-matrices'''
################################################################################
import numpy
from numpy import *
import inspect
from globalz import *
from helpers import *
import sys
import copy
import scipy.sparse as sp

def mkflowmat(m):
    '''takes a model object (m) as argument and returns a list of sparse
    matrices which contain D-values from inter-cell flows for each timestep.'''
    fdict=copy.deepcopy(m.Dflow)
    for tp in fdict.keys(): # check whether we have all compartments
        for c in tp:
            if c not in m.compdict.keys():
                print('''flows.py: Compartment %i not in compartment-list.
                Aborting !''') % (c)
                sys.exit(1)
    # remove intra-cell flows
    # this is for example oceanic sinking flux, which is dealt with
    # by "betr_ocean_sinkflux" in "processes.py".
    for [k,v] in fdict.items():
        intercell=list(where(v[:,0]!=v[:,1])[0])
        intracell=list(where(v[:,0]==v[:,1])[0])
        fdict[k]=v[intercell,:]
     # construct list of large sparse matrices
    matlist=[]
    for t in arange(0,m.nots):
        matlist.append(sp.coo_matrix(matrix(zeros((m.matdim,m.matdim),
                                                  dtype=float32))))  # original: float64 
    for [k,v] in fdict.items():
        if len(v)==0: continue
        toidx=array([(x-1)*m.nocomp+k[1]-1 for x in v[:,1]])
        fromidx=array([(x-1)*m.nocomp+k[0]-1 for x in v[:,0]])
        ij=(toidx.astype(int),fromidx.astype(int))
        count=0
        for ts in arange(2,v.shape[1]):
            matlist[count]=matlist[count]\
                         +sp.coo_matrix((v[:,ts],ij),shape=(m.matdim,m.matdim))
            count+=1
    return(matlist)

def mkprocmat(m):
    ''' returns a list of sparse matrices, one for each timestep, that
    contain D-values for intra-cell processes'''
    for tp in m.Dproc.keys(): # check whether we have all compartments
        tp=(tp[0],tp[1])
        for c in tp:
            if c not in m.compdict.keys():
                print('''mkprocmat.py: Compartment %i not in compartment-list.
                Aborting !''') % (c)
                sys.exit(1)
    # construct list of large sparse matrices
    matlist=[]
    for t in arange(0,m.nots):
        matlist.append(sp.coo_matrix(matrix(zeros((m.matdim,m.matdim),
                                                  dtype=float64))))
    for [k,v] in m.Dproc.items():
        toidx=array([x*m.nocomp+k[1]-1 for x in range(0,m.nocells)])
        fromidx=array([x*m.nocomp+k[0]-1 for x in range(0,m.nocells)])
        ij=(toidx.astype(int),fromidx.astype(int))
        count=0
        for ts in arange(0,v.shape[1]):
            matlist[count]=matlist[count]\
                         +sp.coo_matrix((v[:,ts],ij),shape=(m.matdim,m.matdim))
            count+=1
    return(matlist)

def mkbigD(m, track_flows=False):
    ''' combines matrices for inter-cell flows and intra-cell processes.
    constructs diagonals. Returns list of sparse csc matrices.
    option track_flows:  # SSchenker
        
        * Steinlin type flux integration 
        * creates a copy of the diagonals of BigD to follow mass fluxes
        * Final BigD matrices look like:
            
            -d 0 0 + 0 0 0 0    +    : Influx
            0 -d 0 0 0 0 0 0    -d   : Outflux
            + + -d 0 0 0 0 0    +d   : mirrored Outfluxes 
            0 0 + -d 0 0 0 0
            +d 0 0 0 0 0 0 0 
            0 +d 0 0 0 0 0 0
            0 0 +d 0 0 0 0 0
            0 0 0 +d 0 0 0 0
            
        * killidx should still work as usual but the concentration
          vector needs to be expanded
            
    '''
    matlist_flow=mkflowmat(m)
    matlist_proc=mkprocmat(m)
    # check for same shape of all matrices
    matshp=matlist_flow[0].shape
    for m in matlist_flow+matlist_proc:
        if m.shape != matshp:
            print("mkbigD: inconsistent matrix dimensions. Aborting!\n")
            sys.exit(1)
	
    # add matrices for intercell transport  and intracell processes 
    matlist=[]
    for ts in arange(0,len(matlist_flow)):
        matlist.append(matlist_flow[ts]+matlist_proc[ts])
	
    # construct proper diagonals
    for ts in arange(0,len(matlist)):
        # seperate matrices from their diagonals
        diagonal=diag(matlist[ts].todense())
        matlist[ts]=matlist[ts]-sp.spdiags(diagonal,0,matshp[0],matshp[1],
                                        format="csc")
        # put losses by transport to other cells on diagonal and
        # put back intra-cell loss with negative sign
        loss=sum(matlist[ts].todense(), axis=0)
        matlist[ts]=matlist[ts]\
                     -sp.spdiags(loss,0,matshp[0],matshp[1],format="csc")\
                     -sp.spdiags(diagonal,0,matshp[0],matshp[1],format="csc")
                     
        if track_flows==True:
            diagonals = -sp.spdiags(loss,0,matshp[0],matshp[1],format="csc")\
                    -sp.spdiags(diagonal,0,matshp[0],matshp[1],format="csc")
            matlist[ts] = sp.bmat([[ matlist[ts], None],[diagonals,\
            sp.csc_matrix(numpy.zeros(matshp))]])
                 
    return(matlist)

def mkZVinv(m):
    ''' uses zdict and vdict to construct a list of
    sparse (csr) diagonal matrices A = 1/ZV'''
    # construct list of large sparse matrices
    mlist=[]
    diags=zeros((m.nots,m.matdim))
    for c in m.compdict.keys():
        idx=tocell(arange(1, m.nocells+1),c,m)
        d=m.zdict[c]['bulk']*m.vdict[c]['bulk']
        for t in arange(0,m.nots):
            diags[t,idx]=d[:,t]
    for d in diags:
        v=1/compressvec(d,m.killidx)
        mlist.append(sp.spdiags(v,0,len(v),len(v), format='csr'))
    return(mlist)
