################################################################################
'''Module "output.py", July 2010, is part of
BETR-Research by Harald von Waldow <hvwaldow@chem.ethz.ch>, which is
based on BETR-Global by Matt MacLeod <matthew.macleod@chem.ethz.ch>

This module contains output-routines for BETR-Research'''
################################################################################
from numpy import *
import sys
import os
import csv
import cPickle
#import write_ncfile
#reload(write_ncfile)
#try:
from write_ncfile import *
#except ImportError:
#	from write_ncfile_old import *

from helpers import *
from helpers3 import *
from secondary_emission_reconstruction import reconstruct_dyn_res

def write_cpk_file(fn, out):
    
    fncpk=fn+'_out.cpk'
    if not os.path.exists(os.path.dirname(fncpk)):
        os.mkdir(os.path.dirname(fncpk))
    else:
        if os.path.exists(fncpk):
            print('Attention: overwriting %s\n') % (fncpk)
    f=open(fncpk,'w')
    cPickle.dump(out,f)
    f.close()

def write_output_summary(m, fn, nchemical, nrun, nchemdb, nseasonalparfile, nconstantparfile,\
                       ncompartmentfile, nflowdir, nprocessfile, ncontrolfile):
    if not os.path.exists(os.path.dirname(fn)):
        os.mkdir(os.path.dirname(fn))    
    writer = csv.writer(open(fn, 'w'), delimiter = ' ')
    writer.writerows([["BETR-Global 3.0"],
                      ['runID', nrun], 
                      ['seasfile', nseasonalparfile],
                      ['chemdata', nchemdb], 
                      ['chemnr', nchemical],
                      ['compfile', ncompartmentfile], 
                      ['constparfile', nconstantparfile],
                      ['flowdirectory', nflowdir], 
                      ['procfile', nprocessfile],
                      ['contfile', ncontrolfile]])      
      
def write_output_summary2(m, fn, nemfile):
    writer = csv.writer(open(fn, 'a'), delimiter = ' ')
    writer.writerow(['emisfile', nemfile])      
    
def write_output_summary3(m, fn, nsolvfile):
    writer = csv.writer(open(fn, 'a'), delimiter = ' ')
    writer.writerow(['solvfile', nsolvfile])    
    
def write_output_ss(m, fn, units, netcdf, cpk):     # parameter cpk added by HW
    ''' Output for steady-state results'''
    nlat=int(sqrt(m.matdim/m.nocomp/2))                  # added by HW
    nlon=int(sqrt(m.matdim/m.nocomp/2)*2)                # added by HW
    out={}
    for c in m.compdict.keys():
        out[c]={}
        idx=arange(c-1,m.matdim,len(m.compdict.keys()))
        ## deal with zero-volumes / Z-values
        varr=mean(m.vdict[c]['bulk'],axis=1)
        zv_arr=mean(m.zdict[c]['bulk'],axis=1)*varr
        zerovolidx=where(varr==0)
        varr[zerovolidx]=-9999.99
        zerozvidx=where(zv_arr==0)
        zv_arr[zerozvidx]=-9999.99
        if 'mol' in units:
            out[c]['mol']=m.ss_res[idx]
        if 'kg' in units:
            out[c]['kg']=m.ss_res[idx]*m.chemdict['molmass']/1000.0
        if 'mol_per_m3' in units:
            out[c]['mol_per_m3']=m.ss_res[idx]/varr
            out[c]['mol_per_m3'][zerovolidx]=0
        if 'kg_per_m3' in units:
            out[c]['kg_per_m3']=m.ss_res[idx]*m.chemdict['molmass']/1000.0/varr
            out[c]['kg_per_m3'][zerovolidx]=0
        if 'Pa' in units:
            out[c]['Pa']=m.ss_res[idx]/zv_arr
            out[c]['Pa'][zerozvidx]=0
            
    if not out:
        print("WARNING: empty output requested !\n")
        
    if not os.path.exists(os.path.dirname(fn)):         # moved upwards
        os.mkdir(os.path.dirname(fn))
            
    if cpk:                                             # condition added by HW
        fncpk=fn+'_out.cpk'
        #if not os.path.exists(os.path.dirname(fncpk)):
            #os.mkdir(os.path.dirname(fncpk))
        #else:
        if os.path.exists(fncpk):
            print('Attention: overwriting %s\n') % (fncpk)
        f=open(fncpk,'w')
        cPickle.dump(out,f)
        f.close()
    if netcdf:
        for u in units:
            fnnc=fn+'_'+u+'.nc'
            if os.path.exists(fnnc):
                print('Attention: overwriting %s\n') % (fnnc)
            #outarray=zeros((m.nocomp,nlon,nlat)) #12,24
            outarray=zeros((m.nocomp,nlat,nlon)) #12,24 # Modified by RKG, 28.07.2014
            for c in m.compdict.keys():
                #outarray[c-1,:,:]=out[c][u].reshape(nlon,nlat) #12,24
                outarray[c-1,:,:]=out[c][u].reshape(nlat,nlon) #12,24 # Modified by RKG, 28.07.2014
            writenc(outarray,fnnc,False,True,1,varname='V',unit=u)

def write_output_dyn(m, fn, units, netcdf, cpk):     # parameter cpk added by HW
    ''' Output for dynamic results'''
    nlat=int(sqrt(m.matdim/m.nocomp/2))                  # added by HW
    nlon=int(sqrt(m.matdim/m.nocomp/2)*2)                # added by HW
    out={}
    #print 'write_output_dyn: self.dyn_res:', type(m.dyn_res),m.dyn_res.shape # RKG, 03.07.2014
    timesteps=m.dyn_res.shape[1]
    #print 'write_output_dyn: timesteps = ', timesteps # RKG, 03.07.2014
    periods=int((timesteps-1)/float(m.nots))
    if periods != (timesteps-1)/float(m.nots):
        sys.exit('timesteps of output {0:d} not multiple of '
                 +'seasons ({1:d}): Aborting !'.format(timesteps,m.nots))
    for c in m.compdict.keys():
        out[c]={}
        idx=arange(c-1,m.matdim,len(m.compdict.keys()))
        ## deal with zero volumes / Z-values
        varr=hstack((ones((idx.shape[0],1)),
                     tile(m.vdict[c]['bulk'], (1,periods))))
        zv_arr=hstack((ones((idx.shape[0],1)), 
                       tile(m.zdict[c]['bulk'], (1,periods))))*varr
        zerovolidx=where(varr==0)
        varr[zerovolidx]=-9999.99
        zerozvidx=where(zv_arr==0)
        zv_arr[zerozvidx]=-9999.99

        if 'mol' in units:
            out[c]['mol']=m.dyn_res[idx]
        if 'kg' in units:
            out[c]['kg']=m.dyn_res[idx]*m.chemdict['molmass']/1000.0
        if 'mol_per_m3' in units:
            out[c]['mol_per_m3']=m.dyn_res[idx]/varr
            out[c]['mol_per_m3'][zerovolidx]=0
        if 'kg_per_m3' in units:
            out[c]['kg_per_m3']=m.dyn_res[idx]*m.chemdict['molmass']/1000.0/varr
            out[c]['kg_per_m3'][zerovolidx]=0
        if 'Pa' in units:
            out[c]['Pa']=m.dyn_res[idx]/zv_arr
            out[c]['Pa'][zerozvidx]=0
    
    if not out:
        print("WARNING: empty output requested !\n")
        
    if not os.path.exists(os.path.dirname(fn)):   # moved upwards
        os.mkdir(os.path.dirname(fn))
        
    if cpk:                          # condition added by HW
        fncpk=fn+'_out.cpk'
        #if not os.path.exists(os.path.dirname(fncpk)):
            #os.mkdir(os.path.dirname(fncpk))
        #else:
        if os.path.exists(fncpk):
            print('Attention: overwriting %s\n') % (fncpk)
        f=open(fncpk,'w')
        cPickle.dump(out,f)
        f.close()
    if netcdf:
        for u in units:
            fnnc=fn+'_'+u+'.nc'
            if os.path.exists(fnnc):
                print('Attention: overwriting %s\n') % (fnnc)
            outarray=zeros((timesteps,m.nocomp,nlat,nlon))# 12,24
            for c in m.compdict.keys():
                for t in arange(0,timesteps):
                    outarray[t,c-1,:,:]=out[c][u][:,t].reshape(nlat,nlon) #12,24
            writenc(outarray, fnnc, True, True, 1, varname='V',unit=u)


def write_output_se(m, fn, units, netcdf, cpk, scenario="air"):     # parameter cpk added by HW
    ''' Output for secondary emissions results'''
#    from secondary_emission_reconstruction import reconstruct_dyn_res
    units = ["pe", "se", 'pe_frac', 'se_frac']
    dyn_res_pe, dyn_res_se =  reconstruct_dyn_res(m, scenario=scenario,nseasons=12)
    
    frac_se = dyn_res_se / (dyn_res_pe + dyn_res_se)
    frac_pe = dyn_res_pe / (dyn_res_pe + dyn_res_se)
    
    index_NaN = isnan(frac_se)
    frac_se[index_NaN] = 0
    frac_pe[index_NaN] = 0

    nlat=int(sqrt(m.matdim/m.nocomp/2))                  # added by HW
    nlon=int(sqrt(m.matdim/m.nocomp/2)*2)                # added by HW
    out={}
    timesteps=m.dyn_res.shape[1]
    #print 'write_output_se: timesteps = ', timesteps # RKG, 03.07.2014
    periods=int((timesteps-1)/float(m.nots))
    if periods != (timesteps-1)/float(m.nots):
        sys.exit('timesteps of output {0:d} not multiple of '
                 +'seasons ({1:d}): Aborting !'.format(timesteps,m.nots))
    for c in m.compdict.keys():
        out[c]={}
        idx=arange(c-1,m.matdim,len(m.compdict.keys()))
        ## deal with zero volumes / Z-values
        varr=hstack((ones((idx.shape[0],1)),
                     tile(m.vdict[c]['bulk'], (1,periods))))
        zv_arr=hstack((ones((idx.shape[0],1)), 
                       tile(m.zdict[c]['bulk'], (1,periods))))*varr
        zerovolidx=where(varr==0)
        varr[zerovolidx]=-9999.99
        zerozvidx=where(zv_arr==0)
        zv_arr[zerozvidx]=-9999.99
        
        out[c]['pe']=dyn_res_pe[idx]
        out[c]['se']=dyn_res_se[idx]
        out[c]['pe_frac']=frac_pe[idx]
        out[c]['se_frac']=frac_se[idx]
        

#        if 'mol' in units:
#            out[c]['mol']=m.dyn_res[idx]
#        if 'kg' in units:
#            out[c]['kg']=m.dyn_res[idx]*m.chemdict['molmass']/1000.0
#        if 'mol_per_m3' in units:
#            out[c]['mol_per_m3']=m.dyn_res[idx]/varr
#            out[c]['mol_per_m3'][zerovolidx]=0
#        if 'kg_per_m3' in units:
#            out[c]['kg_per_m3']=m.dyn_res[idx]*m.chemdict['molmass']/1000.0/varr
#            out[c]['kg_per_m3'][zerovolidx]=0
#        if 'Pa' in units:
#            out[c]['Pa']=m.dyn_res[idx]/zv_arr
#            out[c]['Pa'][zerozvidx]=0
    
    if not out:
        print("WARNING: empty output requested !\n")
        
    if not os.path.exists(os.path.dirname(fn)):   # moved upwards
        os.mkdir(os.path.dirname(fn))
        
    if cpk:                          # condition added by HW
        fncpk=fn+'_out.cpk'
        #if not os.path.exists(os.path.dirname(fncpk)):
            #os.mkdir(os.path.dirname(fncpk))
        #else:
        if os.path.exists(fncpk):
            print('Attention: overwriting %s\n') % (fncpk)
        f=open(fncpk,'w')
        cPickle.dump(out,f)
        f.close()
    if netcdf:
        for u in units:
            fnnc=fn+'_'+u+'.nc'
            if os.path.exists(fnnc):
                print('Attention: overwriting %s\n') % (fnnc)
            outarray=zeros((timesteps,m.nocomp,nlat,nlon))# 12,24
            for c in m.compdict.keys():
                for t in arange(0,timesteps):
                    outarray[t,c-1,:,:]=out[c][u][:,t].reshape(nlat,nlon) #12,24
            writenc(outarray, fnnc, True, True, 1, varname='V',unit=u)
            
            
def write_output_dyn_tmp(m, fn, timesteps):                             # function added by HW
    nreg=m.matdim/m.nocomp
    out_tmp = zeros([nreg*7, 3])#288*7
    out_tmp[:, 0] = tile(range(1, nreg+1), 7)#289
    out_tmp[:, 1] = repeat(range(1, 8), nreg)#288
    
    #timesteps=m.dyn_res.shape[1]
    periods=int((timesteps-1)/float(m.nots))
    if periods != (timesteps-1)/float(m.nots):
        sys.exit('timesteps of output {0:d} not multiple of '
                 +'seasons ({1:d}): Aborting !'.format(timesteps,m.nots))
    
    for c in m.compdict.keys():
        idx=arange(c-1,m.matdim,len(m.compdict.keys()))
        ## deal with zero volumes / Z-values
        varr=hstack((ones((idx.shape[0],1)),
                     tile(m.vdict[c]['bulk'], (1,periods))))
        zv_arr=hstack((ones((idx.shape[0],1)),
                        tile(m.zdict[c]['bulk'], (1,periods))))*varr
        zerovolidx=where(varr==0)
        varr[zerovolidx]=-9999.99
        zerozvidx=where(zv_arr==0)
        zv_arr[zerozvidx]=-9999.99
        
        out_tmp[(c-1)*nreg : (c-1)*nreg+nreg, 2] = m.dyn_res[idx, timesteps-1] # unit = mol   #288 
          
    savetxt(fn, out_tmp, fmt = ['%i', '%i', '%.16e'], delimiter=' ')  # changed from .8e to .16e
    
    
def write_output_end_txt(m, fn, timesteps):                             # function added by HW
    nreg=m.matdim/m.nocomp
    out_tmp = zeros([nreg*7, 3])#288*7
    out_tmp[:, 0] = tile(range(1, nreg+1), 7)#289
    out_tmp[:, 1] = repeat(range(1, 8), nreg)#288
    
    #timesteps = m.dyn_res.shape[1]    
    periods=int((timesteps-1)/float(m.nots))
    if periods != (timesteps-1)/float(m.nots):
        sys.exit('timesteps of output {0:d} not multiple of '
                 +'seasons ({1:d}): Aborting !'.format(timesteps,m.nots))
    
    for c in m.compdict.keys():
        idx=arange(c-1,m.matdim,len(m.compdict.keys()))
        ## deal with zero volumes / Z-values
        varr=hstack((ones((idx.shape[0],1)),
                     tile(m.vdict[c]['bulk'], (1,periods))))
        zv_arr=hstack((ones((idx.shape[0],1)),
                        tile(m.zdict[c]['bulk'], (1,periods))))*varr
        zerovolidx=where(varr==0)
        varr[zerovolidx]=-9999.99
        zerozvidx=where(zv_arr==0)
        zv_arr[zerozvidx]=-9999.99
        
        out_tmp[(c-1)*nreg : (c-1)*nreg+nreg, 2] = m.dyn_res[idx, timesteps-1] # unit = mol   #288  #
          
    savetxt(fn, out_tmp, fmt = ['%i', '%i', '%.16e'], delimiter=' ')   # changed from .8e to .16e

def write_output_ss_txt(m, fn):                             # function added by RKG, 04.01.2014
    nreg=m.matdim/m.nocomp
    out_tmp = zeros([nreg*7, 3])#288*7
    out_tmp[:, 0] = tile(range(1, nreg+1), 7)#289
    out_tmp[:, 1] = repeat(range(1, 8), nreg)#288
    
    #timesteps = m.dyn_res.shape[1]    
    #periods=int((timesteps-1)/float(m.nots))
    #if periods != (timesteps-1)/float(m.nots):
    #    sys.exit('timesteps of output {0:d} not multiple of '
    #             +'seasons ({1:d}): Aborting !'.format(timesteps,m.nots))
    
    for c in m.compdict.keys():
        idx=arange(c-1,m.matdim,len(m.compdict.keys()))
        ## deal with zero volumes / Z-values
    #    varr=hstack((ones((idx.shape[0],1)),
    #                 tile(m.vdict[c]['bulk'], (1,periods))))
    #    zv_arr=hstack((ones((idx.shape[0],1)),
    #                    tile(m.zdict[c]['bulk'], (1,periods))))*varr
    #    zerovolidx=where(varr==0)
    #    varr[zerovolidx]=-9999.99
    #    zerozvidx=where(zv_arr==0)
    #    zv_arr[zerozvidx]=-9999.99
        
        out_tmp[(c-1)*nreg : (c-1)*nreg+nreg, 2] = m.ss_res[idx] # unit = mol   #288  #
          
    savetxt(fn, out_tmp, fmt = ['%i', '%i', '%.16e'], delimiter=' ')   

    
def write_bigDlist_txt(bigDlist, fn, m):             # function added by HW
    nreg=m.matdim/m.nocomp
    #dtypeDlist = [('regfrom', int), ('compfrom', int), ('regto', int), ('compto', int), 
                      #(str(1), float), (str(2), float), (str(3), float), (str(4), float), (str(5), float), (str(6), float),
                      #(str(7), float), (str(8), float), (str(9), float), (str(10), float), (str(11), float), (str(12), float)]
    Dflow = zeros((4*nreg*nreg, 16))
    fromreg = hstack([sort(range(1,nreg+1)*nreg)]*4)
    toreg = hstack([range(1,nreg+1)*nreg]*4)
    Dflow[:,range(0,4)] = transpose([fromreg,
                                     hstack([[1]*nreg*nreg, [2]*nreg*nreg, [4]*nreg*nreg, [5]*nreg*nreg]),
                                     toreg,
                                     hstack([[1]*nreg*nreg, [2]*nreg*nreg, [4]*nreg*nreg, [5]*nreg*nreg])])
    Dflow = Dflow[Dflow[:,0] != Dflow[:,2], :]
    for month in range(0,12):
        Dflow[:,month + 4] = bigDlist[month][tocell(Dflow[:,2], Dflow[:,3], m), tocell(Dflow[:,0], Dflow[:,1], m)]
    Dflow = Dflow[findnonemptycells(Dflow[:, range(4,16)]), :]
    
    Ddiag = zeros((7*nreg, 16))
    fromreg = sort(range(1,nreg+1)*7)
    fromcomp = range(1,8)*nreg
    Ddiag[:,range(0,4)] = transpose([fromreg, fromcomp, fromreg, fromcomp])
    for month in range(0,12):
        Ddiag[:,month + 4] = bigDlist[month][tocell(Ddiag[:,2], Ddiag[:,3], m), tocell(Ddiag[:,0], Ddiag[:,1], m)]
    Ddiag = Ddiag[findnonemptycells(Ddiag[:, range(4,16)]), :]
    
    Dint = zeros((nreg*7*7, 16))
    fromreg = sort(range(1,nreg+1)*7*7)
    fromcomp = hstack([sort(range(1,8)*7)]*nreg)
    tocomp = range(1,8)*7*nreg
    Dint[:,range(0,4)] = transpose([fromreg, fromcomp, fromreg, tocomp])
    Dint = Dint[Dint[:,1] != Dint[:,3],:]
    Dint = Dint[Dint[:,0] == Dint[:,2],:]   # unneccesary
    for month in range(0,12):
        Dint[:,month + 4] = bigDlist[month][tocell(Dint[:,2], Dint[:,3], m), tocell(Dint[:,0], Dint[:,1], m)]
    Dint = Dint[findnonemptycells(Dint[:, range(4,16)]), :]
    
    bigDlist2 = vstack([Dflow, Ddiag, Dint])
    #bigDlist2.dtype = dtype
    #bigDlist2 = sort(bigDlist2, order = 'regfrom')
    savetxt(fn, bigDlist2, 
            fmt = ['%i', '%i', '%i', '%i', '%.8e', '%.8e', '%.8e', '%.8e', '%.8e', '%.8e', '%.8e', '%.8e', '%.8e', '%.8e', '%.8e', '%.8e'],            
            delimiter=' ')


def write_output_dflux(m,
                       fn,
                       netcdf=False,
                       units=['mol'],
                       nlat= 12,
                       nlon=24,
                       ):
    '''Output all files related to the integrated fluxes'''                                
#    possiblefiles = ["flow", "airwater", "deg", "sed", "dep_air", "dep_water"]
    
    out = mk_allfluxes(m)
    #print 'write_output_dflux: out =', type(out), out # RKG, 03.07.2014
    #print 'write_output_dflux: out.keys() =', out.keys() # RKG, 03.07.2014
    k = out.keys()[0]
    #print 'write_output_dflux: k = ',k # RKG, 03.07.2014
    timesteps=len(out[k])
    #print 'write_output_dflux: timesteps = ', timesteps # RKG, 03.07.2014
    #print 'write_output_dflux: timesteps = ', type(out[k]) # RKG, 25.08.2014
    #print 'out["air_to_soil"][0]', shape(out["air_to_soil"][0]),out["air_to_soil"][0] # RKG, 25.08.2014
#    if "flow" in proplist:
#        mk_flow(m)
#    if "airwater" in proplist:
#        mk_airwater(m)
#    if "deg" in proplist:
#        mk_deg(m)
#    if "sed" in proplist:
#        mk_sed(m)
#    if "dep_air" in proplist:
#        mk_dep_air(m)
#    if "dep_water" in proplist:
#        mk_dep_water(m)
    write_cpk_file(fn, out)
    if netcdf:
        for u in out.keys():
            fnnc=fn+'_'+u+'.nc'
            if os.path.exists(fnnc):
                print('Attention: overwriting %s\n') % (fnnc)
            outarray=zeros((timesteps,1,nlat,nlon))
            for t in arange(0,timesteps):
                    outarray[t,0,:,:]=out[u][t].reshape(nlat,nlon)
            writenc(outarray, fnnc, True, True, 1, varname='V',unit=u)
    if netcdf: # The section below is added by RKG, 25.08.2014, to write the deposition fluxes in units of ng/m2/h
        areal_rate = ["air_to_veg", "veg_to_air", "air_to_freshwater", "freshwater_to_air",
                      "air_to_ocean", "ocean_to_air", "air_to_soil", "soil_to_air",        
                      "air_to_all", "all_to_air", "air_to_all_net"]
        for u in out.keys():
            if (u in areal_rate):
                for t in arange(0,timesteps):
                    for cell in arange(0,m.nocells):
                        #print t, cell, t%12
                        if out[u][t][cell] > 0.0:
                            if u in ["air_to_all", "all_to_air", "air_to_all_net", "upper_air_to_lower_air"]: # divide by the whole region area
                                out[u][t][cell] = out[u][t][cell] * m.chemdict['molmass'] * 1.0E+9 / m.par['A'][cell][t%12] / 730.0 # Change 'mol' to 'ng/m2-all/h'
                            if u in ["air_to_veg", "veg_to_air"]: # divide by the vegetated area
                                out[u][t][cell] = out[u][t][cell] * m.chemdict['molmass'] * 1.0E+9 / (m.par['A'][cell][t%12]*m.par['perc6'][cell][t%12]*m.par['perc3'][cell][t%12]) / 730.0 # Change 'mol' to 'ng/m2-veg/h'
                                #print u, m.par['perc6'][cell][t%12], m.par['perc3'][cell][t%12]
                            if u in ["air_to_freshwater", "freshwater_to_air"]: # divide by the freshwater area
                                out[u][t][cell] = out[u][t][cell] * m.chemdict['molmass'] * 1.0E+9 / (m.par['A'][cell][t%12]*m.par['perc4'][cell][t%12]) / 730.0 # Change 'mol' to 'ng/m2-freshwater/h'
                                #print u, m.par['perc4'][cell][t%12]
                            if u in ["air_to_ocean", "ocean_to_air"]: # divide by the ocean area
                                out[u][t][cell] = out[u][t][cell] * m.chemdict['molmass'] * 1.0E+9 / (m.par['A'][cell][t%12]*m.par['perc5'][cell][t%12]) / 730.0 # Change 'mol' to 'ng/m2-ocean/h'
                                #print u, m.par['perc5'][cell][t%12]
                            if u in ["air_to_soil", "soil_to_air"]: # divide by the soil area
                                out[u][t][cell] = out[u][t][cell] * m.chemdict['molmass'] * 1.0E+9 / (m.par['A'][cell][t%12]*m.par['perc6'][cell][t%12]) / 730.0 # Change 'mol' to 'ng/m2-soil/h
                                #print u, m.par['perc6'][cell][t%12]
                fnnc=fn+'_arate_'+u+'.nc'
                if os.path.exists(fnnc):
                    print('Attention: overwriting %s\n') % (fnnc)
                outarray=zeros((timesteps,1,nlat,nlon))
                for t in arange(0,timesteps):
                    outarray[t,0,:,:]=out[u][t].reshape(nlat,nlon)
                writenc(outarray, fnnc, True, True, 1, varname=u,unit='ng/m2/h')
        molar_rate = ["flow_out", "flow_in", "flow_in_lower_air", "flow_in_others","upper_air_to_lower_air", # added by RKG, 28.08.2014
                      "air_to_veg", "veg_to_air", "air_to_freshwater", "freshwater_to_air", # added by RKG, 09.04.2014
                      "air_to_ocean", "ocean_to_air", "air_to_soil", "soil_to_air",        
                      "air_to_all", "all_to_air", "air_to_all_net",
                      "flow_out_lower_air", "flow_in_upper_air", "flow_out_upper_air", # added by RKG, 27.11.2014
                      "flow_in_freshwater", "flow_out_freshwater",
                      "flow_in_ocean", "flow_out_ocean",
                      "flow_out_air", "flow_in_air",
                      "deg_upper_air", "deg_lower_air", "deg_air", "deg_vegetation",
                      "deg_freshwater", "deg_ocean", "deg_soil", "deg_sediment",
                      "lower_air_to_upper_air",
                      "air_to_sediment", "sediment_to_air", "air_to_sediment_net",
                      "flow_in_soil", "flow_out_soil", "flow_in_sediment", "flow_out_sediment", # added by RKG, 17.12.2014
                      "soil_to_freshwater", "freshwater_to_soil", "soil_to_freshwater_net",
                      "soil_to_ocean", "ocean_to_soil", "soil_to_ocean_net",
                      "sediment_to_freshwater", "freshwater_to_sediment", "sediment_to_freshwater_net",
                      "sediment_to_ocean", "ocean_to_sediment", "sediment_to_ocean_net", 
                      "sed_upper_air", "sed_lower_air", "sed_air", "sed_vegetation",
                      "sed_freshwater", "sed_ocean", "sed_soil", "sed_sediment", # last 2 lines added by RKG, 07.02.2015
                      "freshwater_to_ocean"] # last line added by RKG, 02.06.2015

      # commenting out the section for writing mrate files. It contains an error. RKG, 09.09.2015
      #  for u in out.keys():
       #     if (u in molar_rate):
        #        for t in arange(0,timesteps):
         #           for cell in arange(0,m.nocells):
          #              #print t, cell, t%12
           #             if out[u][t][cell] > 0.0:
            #                out[u][t][cell] = out[u][t][cell] / 730.0 # Change 'mol' to 'mol/h' ##THIS IS ERRONEOUS!!! RKG, 09.09.2015
                                                                                                 ## out array has values in units of ng/m2/h. See code above!!!
             #   fnnc=fn+'_mrate_'+u+'.nc'
              #  if os.path.exists(fnnc):
               #     print('Attention: overwriting %s\n') % (fnnc)
               # outarray=zeros((timesteps,1,nlat,nlon))
                #for t in arange(0,timesteps):
                 #   outarray[t,0,:,:]=out[u][t].reshape(nlat,nlon)
                #writenc(outarray, fnnc, True, True, 1, varname=u,unit='mol/h')

        
    # Go through all compartments and summarize the expected properties
##############################################################################3
#SSchenker    The following functions are needed for Linear Approx. Fluxes
#             and for the follow_flow mode (Correct fluxes)
    
def mk_allfluxes(m, nseasons=12):
    '''Calculates the net in/outflow from all compartments at every timestep'''
    air = [0,1]
    stationary = [2,5,6]
    water = [3,4]
    flux_key = normalizefluxkey(m.flux_key)
    
    fluxdict = {}
    fluxes=["flow", "air_to_water", "deg", "sed", "air_to_x", "water_to_x",
            "air_to_veg", "veg_to_air", "air_to_freshwater", "freshwater_to_air", # last 3 lines added by RKG, 25.08.2014
            "air_to_ocean", "ocean_to_air", "air_to_soil", "soil_to_air",         # These are all gross fluxes
            "air_to_all", "all_to_air", "air_to_all_net",                         # except: "air_to_all_net"
            "upper_air_to_lower_air",                                             # this one added by RKG, 28.08.2014
            "flow_in", "flow_out", "flow_in_lower_air",# "flow_in_others"]         # these are added by RKG, 28.08.2014
            "flow_out_lower_air", "flow_in_upper_air", "flow_out_upper_air",
            "flow_in_freshwater", "flow_out_freshwater",
            "flow_in_ocean", "flow_out_ocean",
            "flow_out_air", "flow_in_air",
            "deg_upper_air", "deg_lower_air", "deg_air", "deg_vegetation",
            "deg_freshwater", "deg_ocean", "deg_soil", "deg_sediment",
            "lower_air_to_upper_air",
            "air_to_sediment", "sediment_to_air", "air_to_sediment_net",        #last 8 lines added by RKG, 27.11.2014
            "flow_in_soil", "flow_out_soil", "flow_in_sediment", "flow_out_sediment",
            "soil_to_freshwater", "freshwater_to_soil", "soil_to_freshwater_net",
            "soil_to_ocean", "ocean_to_soil", "soil_to_ocean_net",
            "sediment_to_freshwater", "freshwater_to_sediment", "sediment_to_freshwater_net",
            "sediment_to_ocean", "ocean_to_sediment", "sediment_to_ocean_net",  # last 5 lines added by RKG, 17.12.2014
            "sed_upper_air", "sed_lower_air", "sed_air", "sed_vegetation",
            "sed_freshwater", "sed_ocean", "sed_soil", "sed_sediment", # last 2 lines added by RKG, 07.02.2015
            "freshwater_to_ocean"] # this line added by RKG, 02.06.2015

    for key in fluxes:
        fluxdict[key] = zeros((len(m.flux_res), m.nocells))    
        
    compmap = mk_compmap(m)
    ts0 = time.time()
    ts1 = time.time()
    for ts in range(len(m.flux_res)):
        season = ts % nseasons
        year = ts / nseasons
        if season ==0:
            print ("Processing year %i \t" %(year)),
        
        fluxmat = m.flux_res[ts]
        x,y = fluxmat.nonzero()
        for i in range(len(x)):
            value = fluxmat[x[i],y[i]]
            
            xcell, xcomp = x_to_cellcomp( x[i], compmap, m.nocomp)
            ycell, ycomp = x_to_cellcomp( y[i], compmap, m.nocomp)
            
            if x[i]==y[i]:
                try :
                    match_dict = m.flux_key[season][xcell][xcomp]
                except KeyError:
                    print [[season],[xcell],[xcomp]]
                for mnkey in match_dict.keys():
                    if type(mnkey) is tuple:
                        #print "mk_allfluxes:", i, "mnkey: ", mnkey, "xcell: ",xcell, "xcomp: ", xcomp
                        fluxdict["flow"][ts][xcell] += flux_key[season][xcell][xcomp][mnkey] * value
                        fluxdict["flow"][ts][mnkey[1]-1] -= flux_key[season][xcell][xcomp][mnkey] * value
                    if mnkey  == "deg":
                        fluxdict["deg"][ts][xcell] += flux_key[season][xcell][xcomp][mnkey] * value
                    if mnkey == "sed":
                        fluxdict["sed"][ts][xcell] += flux_key[season][xcell][xcomp][mnkey] * value
                    
                    if xcomp in air and mnkey in water:
                        fluxdict["air_to_water"][ts][xcell] -= flux_key[season][xcell][xcomp][mnkey] * value
                    if xcomp in water and mnkey in air:
                        fluxdict["air_to_water"][ts][xcell] += flux_key[season][xcell][xcomp][mnkey] * value
                    
                    if xcomp in air and mnkey in stationary:
                        fluxdict["air_to_x"][ts][xcell] -= flux_key[season][xcell][xcomp][mnkey] * value
                    if xcomp in stationary and mnkey in air:
                        fluxdict["air_to_x"][ts][xcell] += flux_key[season][xcell][xcomp][mnkey] * value
                    
                    if xcomp in water and mnkey in stationary:
                        fluxdict["water_to_x"][ts][xcell] -= flux_key[season][xcell][xcomp][mnkey] * value
                    if xcomp in stationary and mnkey in water:
                        fluxdict["water_to_x"][ts][xcell] += flux_key[season][xcell][xcomp][mnkey] * value
                # new calculations added by RKG, 25.08.2014
                    if xcomp in air and mnkey in [2]:
                        fluxdict["air_to_veg"][ts][xcell] -= flux_key[season][xcell][xcomp][mnkey] * value
                    if xcomp in [2] and mnkey in air:
                        fluxdict["veg_to_air"][ts][xcell] -= flux_key[season][xcell][xcomp][mnkey] * value
                    if xcomp in air and mnkey in [3]:
                        fluxdict["air_to_freshwater"][ts][xcell] -= flux_key[season][xcell][xcomp][mnkey] * value
                    if xcomp in [3] and mnkey in air:
                        fluxdict["freshwater_to_air"][ts][xcell] -= flux_key[season][xcell][xcomp][mnkey] * value
                    if xcomp in air and mnkey in [4]:
                        fluxdict["air_to_ocean"][ts][xcell] -= flux_key[season][xcell][xcomp][mnkey] * value
                    if xcomp in [4] and mnkey in air:
                        fluxdict["ocean_to_air"][ts][xcell] -= flux_key[season][xcell][xcomp][mnkey] * value
                    if xcomp in air and mnkey in [5]:
                        fluxdict["air_to_soil"][ts][xcell] -= flux_key[season][xcell][xcomp][mnkey] * value
                    if xcomp in [5] and mnkey in air:
                        fluxdict["soil_to_air"][ts][xcell] -= flux_key[season][xcell][xcomp][mnkey] * value
                    if xcomp in air and mnkey in [2,3,4,5,6]:
                        fluxdict["air_to_all"][ts][xcell] -= flux_key[season][xcell][xcomp][mnkey] * value
                    if xcomp in [2,3,4,5,6] and mnkey in air:
                        fluxdict["all_to_air"][ts][xcell] -= flux_key[season][xcell][xcomp][mnkey] * value
                        
                    if xcomp in air and mnkey in [2,3,4,5,6]:
                        fluxdict["air_to_all_net"][ts][xcell] -= flux_key[season][xcell][xcomp][mnkey] * value
                    if xcomp in [2,3,4,5,6] and mnkey in air:
                        fluxdict["air_to_all_net"][ts][xcell] += flux_key[season][xcell][xcomp][mnkey] * value
                # new calculations added by RKG, 28.08.2014
                    if xcomp in [0] and mnkey in [1]:
                        fluxdict["upper_air_to_lower_air"][ts][xcell] -= flux_key[season][xcell][xcomp][mnkey] * value
                    if type(mnkey) is tuple:
                        #print "mk_allfluxes:", i, "mnkey: ", mnkey, "xcell: ",xcell, "xcomp: ", xcomp
                        fluxdict["flow_out"][ts][xcell] -= flux_key[season][xcell][xcomp][mnkey] * value
                        fluxdict["flow_in"][ts][mnkey[1]-1] -= flux_key[season][xcell][xcomp][mnkey] * value
                        if xcomp in [1]:
                            #print "flow_in_and_out_lower_air"
                            #print "mk_allfluxes:", i, "mnkey: ", mnkey, "xcell: ",xcell, "xcomp: ", xcomp
                            fluxdict["flow_in_lower_air"][ts][mnkey[1]-1] -= flux_key[season][xcell][xcomp][mnkey] * value
                # new calculations added by RKG as per LC's request, 28.11.2014
                            fluxdict["flow_out_lower_air"][ts][xcell] -= flux_key[season][xcell][xcomp][mnkey] * value
                        if xcomp in [0]:
                            fluxdict["flow_in_upper_air"][ts][mnkey[1]-1] -= flux_key[season][xcell][xcomp][mnkey] * value
                            fluxdict["flow_out_upper_air"][ts][xcell] -= flux_key[season][xcell][xcomp][mnkey] * value
                        if xcomp in [3]:
                            #print "flow_in_and_out_freshwater"
                            if mnkey[1]-1 == xcell:
                                print "mk_allfluxes:", i, "mnkey: ", mnkey, "xcell: ",xcell, "xcomp: ", xcomp
                            fluxdict["flow_in_freshwater"][ts][mnkey[1]-1] -= flux_key[season][xcell][xcomp][mnkey] * value
                            fluxdict["flow_out_freshwater"][ts][xcell] -= flux_key[season][xcell][xcomp][mnkey] * value
                        if xcomp in [4]:
                            #print "flow_in_and_out_ocean"
                            if mnkey[1]-1 == xcell:
                                print "mk_allfluxes:", i, "mnkey: ", mnkey, "xcell: ",xcell, "xcomp: ", xcomp
                            fluxdict["flow_in_ocean"][ts][mnkey[1]-1] -= flux_key[season][xcell][xcomp][mnkey] * value
                            fluxdict["flow_out_ocean"][ts][xcell] -= flux_key[season][xcell][xcomp][mnkey] * value
                        if xcomp in [5]: # added by RKG, 17.12.2014
                            fluxdict["flow_in_soil"][ts][mnkey[1]-1] -= flux_key[season][xcell][xcomp][mnkey] * value
                            fluxdict["flow_out_soil"][ts][xcell] -= flux_key[season][xcell][xcomp][mnkey] * value
                        if xcomp in [6]: # added by RKG, 17.12.2014
                            fluxdict["flow_in_sediment"][ts][mnkey[1]-1] -= flux_key[season][xcell][xcomp][mnkey] * value
                            fluxdict["flow_out_sediment"][ts][xcell] -= flux_key[season][xcell][xcomp][mnkey] * value    
                        #else:
                        #    fluxdict["flow_in_others"][ts][mnkey[1]-1] -= flux_key[season][xcell][xcomp][mnkey] * value
                        if xcomp in air:
                            fluxdict["flow_out_air"][ts][xcell] -= flux_key[season][xcell][xcomp][mnkey] * value
                            fluxdict["flow_in_air"][ts][mnkey[1]-1] -= flux_key[season][xcell][xcomp][mnkey] * value

                    if mnkey  == "deg":
                        if xcomp == 0:
                            fluxdict["deg_upper_air"][ts][xcell] += flux_key[season][xcell][xcomp][mnkey] * value
                        if xcomp == 1:
                            fluxdict["deg_lower_air"][ts][xcell] += flux_key[season][xcell][xcomp][mnkey] * value
                        if xcomp in air:
                            fluxdict["deg_air"][ts][xcell] += flux_key[season][xcell][xcomp][mnkey] * value
                        if xcomp == 2:
                            fluxdict["deg_vegetation"][ts][xcell] += flux_key[season][xcell][xcomp][mnkey] * value
                        if xcomp == 3:
                            fluxdict["deg_freshwater"][ts][xcell] += flux_key[season][xcell][xcomp][mnkey] * value
                        if xcomp == 4:
                            fluxdict["deg_ocean"][ts][xcell] += flux_key[season][xcell][xcomp][mnkey] * value
                        if xcomp == 5:
                            fluxdict["deg_soil"][ts][xcell] += flux_key[season][xcell][xcomp][mnkey] * value
                        if xcomp == 6:
                            fluxdict["deg_sediment"][ts][xcell] += flux_key[season][xcell][xcomp][mnkey] * value

                    if xcomp in [1] and mnkey in [0]:
                        fluxdict["lower_air_to_upper_air"][ts][xcell] -= flux_key[season][xcell][xcomp][mnkey] * value

                    if xcomp in air and mnkey in [6]:
                        fluxdict["air_to_sediment"][ts][xcell] -= flux_key[season][xcell][xcomp][mnkey] * value
                    if xcomp in [6] and mnkey in air:
                        fluxdict["sediment_to_air"][ts][xcell] -= flux_key[season][xcell][xcomp][mnkey] * value 
                    if xcomp in air and mnkey in [6]:
                        fluxdict["air_to_sediment_net"][ts][xcell] -= flux_key[season][xcell][xcomp][mnkey] * value
                    if xcomp in [6] and mnkey in air:
                        fluxdict["air_to_sediment_net"][ts][xcell] += flux_key[season][xcell][xcomp][mnkey] * value  
                            
              # new calculations added by RKG as per LC's request, 17.12.2014
                    if xcomp in [5] and mnkey in [3]:
                        fluxdict["soil_to_freshwater"][ts][xcell] -= flux_key[season][xcell][xcomp][mnkey] * value
                    if xcomp in [3] and mnkey in [5]:
                        fluxdict["freshwater_to_soil"][ts][xcell] -= flux_key[season][xcell][xcomp][mnkey] * value 
                    if xcomp in [5] and mnkey in [3]:
                        fluxdict["soil_to_freshwater_net"][ts][xcell] -= flux_key[season][xcell][xcomp][mnkey] * value
                    if xcomp in [3] and mnkey in [5]:
                        fluxdict["soil_to_freshwater_net"][ts][xcell] += flux_key[season][xcell][xcomp][mnkey] * value

                    if xcomp in [5] and mnkey in [4]:
                        fluxdict["soil_to_ocean"][ts][xcell] -= flux_key[season][xcell][xcomp][mnkey] * value
                    if xcomp in [4] and mnkey in [5]:
                        fluxdict["ocean_to_soil"][ts][xcell] -= flux_key[season][xcell][xcomp][mnkey] * value 
                    if xcomp in [5] and mnkey in [4]:
                        fluxdict["soil_to_ocean_net"][ts][xcell] -= flux_key[season][xcell][xcomp][mnkey] * value
                    if xcomp in [4] and mnkey in [5]:
                        fluxdict["soil_to_ocean_net"][ts][xcell] += flux_key[season][xcell][xcomp][mnkey] * value

                    if xcomp in [6] and mnkey in [3]:
                        fluxdict["sediment_to_freshwater"][ts][xcell] -= flux_key[season][xcell][xcomp][mnkey] * value
                    if xcomp in [3] and mnkey in [6]:
                        fluxdict["freshwater_to_sediment"][ts][xcell] -= flux_key[season][xcell][xcomp][mnkey] * value 
                    if xcomp in [6] and mnkey in [3]:
                        fluxdict["sediment_to_freshwater_net"][ts][xcell] -= flux_key[season][xcell][xcomp][mnkey] * value
                    if xcomp in [3] and mnkey in [6]:
                        fluxdict["sediment_to_freshwater_net"][ts][xcell] += flux_key[season][xcell][xcomp][mnkey] * value

                    if xcomp in [6] and mnkey in [4]:
                        fluxdict["sediment_to_ocean"][ts][xcell] -= flux_key[season][xcell][xcomp][mnkey] * value
                    if xcomp in [4] and mnkey in [6]:
                        fluxdict["ocean_to_sediment"][ts][xcell] -= flux_key[season][xcell][xcomp][mnkey] * value 
                    if xcomp in [6] and mnkey in [4]:
                        fluxdict["sediment_to_ocean_net"][ts][xcell] -= flux_key[season][xcell][xcomp][mnkey] * value
                    if xcomp in [4] and mnkey in [6]:
                        fluxdict["sediment_to_ocean_net"][ts][xcell] += flux_key[season][xcell][xcomp][mnkey] * value
                    # new calculations added by RKG to check river to ocean flows, 02.06.2015
                    if xcomp in [3] and mnkey in [4]:
                        #print "xcell:",xcell,"freshwater to ocean:", flux_key[season][xcell][xcomp][mnkey] * value
                        fluxdict["freshwater_to_ocean"][ts][xcell] -= flux_key[season][xcell][xcomp][mnkey] * value 
                    #if xcomp in [4] and mnkey in [3]:
                        #print "xcell:",xcell,"ocean to freshwater:", flux_key[season][xcell][xcomp][mnkey] * value
                    # new calculations added by RKG as per LC's request, 07.02.2015
                    if mnkey  == "sed":
                        if xcomp == 0:
                            fluxdict["sed_upper_air"][ts][xcell] += flux_key[season][xcell][xcomp][mnkey] * value
                        if xcomp == 1:
                            fluxdict["sed_lower_air"][ts][xcell] += flux_key[season][xcell][xcomp][mnkey] * value
                        if xcomp in air:
                            fluxdict["sed_air"][ts][xcell] += flux_key[season][xcell][xcomp][mnkey] * value
                        if xcomp == 2:
                            fluxdict["sed_vegetation"][ts][xcell] += flux_key[season][xcell][xcomp][mnkey] * value
                        if xcomp == 3:
                            fluxdict["sed_freshwater"][ts][xcell] += flux_key[season][xcell][xcomp][mnkey] * value
                        if xcomp == 4:
                            fluxdict["sed_ocean"][ts][xcell] += flux_key[season][xcell][xcomp][mnkey] * value
                        if xcomp == 5:
                            fluxdict["sed_soil"][ts][xcell] += flux_key[season][xcell][xcomp][mnkey] * value
                        if xcomp == 6:
                            fluxdict["sed_sediment"][ts][xcell] += flux_key[season][xcell][xcomp][mnkey] * value

                        
                    
                    
                
                        
                
                        
# SSchenker: flows are now calculated from the compartment they origin
#            considering the off diagonal elements is no longer nescessary
#
#            elif x[i]!=y[i] and xcell!=ycell:
#                    fluxdict["flow"][ts][xcell] +=  value
        print ("."),
        if season == (nseasons-1):
            print('  [%.3f s]')  % (time.time()-ts1)
            ts1=time.time()
    print ("Time in mk_allfluxes(): %.3f min " % ((time.time()-ts0)/60.))
    return fluxdict        
                      
                
                
#def  x_to_cellcomp(x, cmap, ncomp ):
#    x = cmap[x]
#    cell = x/ncomp
#    comp = x % ncomp
#    return cell, comp    
#
#def normalizefluxkey(flux_key):
#    for i in range(len(flux_key)):
#        for ii in range(len(flux_key[i])):
#            for iii in range(len(flux_key[i][ii])):
#                fsum=0.0
#                for key in flux_key[i][ii][iii].keys():
#                    fsum += flux_key[i][ii][iii][key]
#                if fsum !=0.0:
#                    for key in flux_key[i][ii][iii].keys():
#                        flux_key[i][ii][iii][key] = flux_key[i][ii][iii][key]/fsum
#    return flux_key
#                    
#    
#    
#def mk_airwater(m):
#    raise NotImplementedError    
#    
#def mk_deg(m):
#    raise NotImplementedError    
#    
#def mk_sed(m):
#    raise NotImplementedError    
#    
#def mk_dep_air(m):
#    raise NotImplementedError    
#    
#def mk_dep_water(m):
#    raise NotImplementedError    
#    
#def mk_compmap(m):
#    nr=0
#    compmap = []
#    for i in range(m.nocells*m.nocomp):
#        if i in m.killidx:
#            nr+=1
#        else:
#            compmap.append(nr)
#            nr+=1
#            
#    return compmap 
#    
#            
#################################################################################

        
