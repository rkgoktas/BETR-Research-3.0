## This is a run for BETR-Reserach with annually changing env.-parameters ####
import os
import time
import BETRS
reload(BETRS)
from BETRS import *

t_s = time.time() # Start time

"""
SET OPTIONS FOR THE FAST SOLVER AND THE FLUX_INTEGRATION AND
PRIMARY/SECONDARY_EMISSION MODE
"""

use_odespy = True       #SWITCH ON/OFF FAST SOLVER
track_fluxes = True     #SWITCH ON/OFF FLUX INTEGRATION
track_se = True         #SWICH ON/OFF TRACKING OF SECONDARY EMISSIONS

## options
"""
Change to current Directory to ensure that relative paths are set correctly

Tends to cause problem under Windows otherwise

"""

abspath=os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)


runID = runID = ['PCB153_ivm2_corrdiffstrato_hr'] # output names 
years = [range(1930, 1935)]*len(runID)    # range of modeling run (years  A.D.)

emisdir = ['PCB153_HIGH_0375_mol_per_hour']*len(runID)  # emission inventory ('Emissions/annual/') 
seasparfile = ['seasonal_parameters_20C3M.a2.comp.ivm.moh_add2mix_hr.txt']*len(runID)#  seasonally varying parameters ('Environment/)
constparfile = ['const_parameters_20C3M.comp.ivm_corrdiffstrato_hr.txt']*len(runID)  # seasonally constant parameters ('Environment/')
flowdirectory = ['OcMix.At100_20C3M.a2.comp.ivm_hr']*len(runID)  # flows in the atmosphere, ocean and fresh water ('Flows/)

chemdata = ['chemicals_v4.txt']*len(runID)  # chemical properties ('Chemicals/')
chemnr = [13]*len(runID)  # selection of chemical from chemical properties files
compfile = ['compartments_4T.txt']*len(runID)   # compartments used in the model ('Environment/') 

procfile = ['processes_default.txt']*len(runID)   # processes used in the model  ('Processes/')
contfile = ['control_default.txt']*len(runID)      # some options ('Control/')
solvfile = ['solvparams_default.txt']*len(runID)    # options for ODE solver  ('Solver/')
mkendfile = False


for v in [years, emisdir, seasparfile, chemdata, chemnr, constparfile, compfile, flowdirectory, procfile, contfile, solvfile]:
    if len(v) != len(runID):
        sys.exit('Warning: one of your input lists is not of same length as number of runs specified')        
    

## now run the model

for i in range(0, len(runID)):
    print('\n\nStarting run ' + runID[i])
    ## model first year and write temporary result to text file
    print('\n\nBETR run ' + runID[i] + ' for year ' + str(years[i][0]))
    m=Model(chemical = chemnr[i], 
            run = runID[i],
            chemdb = chemdata[i],
            #seasonalparfile = os.path.join('annual', seasdir[i], str(years[i][0]) + '.txt'),
            seasonalparfile = seasparfile[i],    
            constantparfile = constparfile[i], 
            compartmentfile = compfile[i], 
            #flowdir = os.path.join('annual', flowdirectory[i], str(years[i][0])),
            flowdir = os.path.join(flowdirectory[i]),
            processfile = procfile[i], 
            controlfile=contfile[i],
            track_flows = (track_fluxes or track_se)
           )
    
    m.update_emissions(os.path.join('annual', emisdir[i], str(years[i][0]) + '.txt'))
    m.update_solver(solvparamfile = solvfile[i])
    m.solve_dyn(1,use_odespy=use_odespy) 
    if mkendfile == True:
            m.output_end_txt(filename = 'endstate.' + str(years[i][0]) + '.txt')
    result = m.dyn_res                                                          # save model result from year 1
    flux_rkg = m.flux_res # experimenting to get the whole fluxes written, RKG, 04.07.2014
    m.output_dyn_tmp(filename = runID[i] + '_tmp.txt')                             # save end state of year 1 as text file
    
    ## model year 2 to n
    for y in years[i][1:]:
        del(m)
        print('\n\nBETR run ' + runID[i] + ' for year ' + str(y))
        m=Model(chemical = chemnr[i], 
            run = runID[i],
            chemdb = chemdata[i],
            #seasonalparfile = os.path.join('annual', seasdir[i], str(y) + '.txt'),
            seasonalparfile = seasparfile[i],
            constantparfile = constparfile[i], 
            compartmentfile = compfile[i], 
            #flowdir = os.path.join('annual', flowdirectory[i], str(y)),
            flowdir = os.path.join(flowdirectory[i]),
            processfile = procfile[i], 
            controlfile=contfile[i],
            track_flows = (track_fluxes or track_se)
            )       
        
        m.update_emissions(os.path.join('annual', emisdir[i], str(y) + '.txt')) # update solver with end state of year - 1 from text file
        m.update_solver(solvparamfile = solvfile[i], initfile = runID[i] + '_tmp.txt')
        m.solve_dyn(1,use_odespy=use_odespy)
        if mkendfile == True:
            m.output_end_txt(filename = 'endstate.' + str(y) + '.txt')
        if y == years[i][len(years[i])-1]: # If this is the last year, save the end state.
            m.output_end_txt(filename = 'endstate.' + str(y) + '.txt')
        result = hstack((result, m.dyn_res[:, 1:13]))                           # save model result from year y
        flux_rkg = hstack((flux_rkg, m.flux_res)) # experimenting to get the whole fluxes written, RKG, 04.07.2014
        m.output_dyn_tmp(filename = runID[i] + '_tmp.txt')                         # save end state of year y as text file
    
    
    ## nc-file output
    m.dyn_res = result           
    m.output_dyn(filename = 'OUT', units = ['mol','kg','mol_per_m3','kg_per_m3','Pa'], netcdf = True, cpk = False)
    os.remove(os.path.join('Solver', runID[i] + '_tmp.txt'))                    # remove temporary file
    if track_fluxes:
        m.flux_res = flux_rkg # experimenting to get the whole fluxes written, RKG, 04.07.2014
        m.output_fluxes(netcdf=True)
    if track_se:
        m.output_se(cpk=False, netcdf=True,scenario="air")

t_e = time.time() # End time

print 'Simulation completed!'
print runID
print 'TOTAL SIMULATION TIME: %f minutes = %f hours.' % ( (t_e - t_s) / 60.0 , (t_e - t_s) / 3600.0 )

