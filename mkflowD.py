################################################################################
'''Module "mkflowsD", July 2010, is part of
BETR-Research by Harald von Waldow <hvwaldow@chem.ethz.ch>, which is
based on BETR-Global by Matt MacLeod <matthew.macleod@chem.ethz.ch>

This module calculates D-values for inter-cell transport processes of a
particular model parametrization.'''
################################################################################
from numpy import *
import inspect
from globalz import *
import sys
import copy

def mkflowD(model):
    fdict=copy.deepcopy(model.flowdict)
    for f in fdict.keys():
        fromcells=fdict[f][:,0].astype(int)
        zvals=model.zdict[f[0]]['bulk'][fromcells-1,:]
        fdict[f][:,2:]=fdict[f][:,2:]*zvals
    return(fdict)
