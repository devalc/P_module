# -*- coding: utf-8 -*-
"""
Created on Tue Apr 03 12:28:10 2018

This file has all the formulations of P module in SWAT, as explained in 2009 
theoretical description document, put into a block of python code.

@author: Chinmay Deval
"""

import os
import numpy as np
import pandas as pd

#############################################################################
                        #Initial P states
#############################################################################


"""
#############################################################################
                        ###Organic P###
#############################################################################
                        
Users are allowed to specify soluble P and organic P contained in humic substances 
for all soil layers 
#The concentration of solution phosphorus in all layers is initially set to 5 
mg/kg soil ( representative of unmanaged land under native vegetation). 
A  concentration  of  25  mg/kg  soil  in  the  plow  layer  is  considered 
representative of cropland."""

"""
Organic phosphorus levels are assigned assuming that the N:P ratio for 
humic materials is 8:1
humic_orgN: Concentration of humic organic N in the layer (mg/kg) 
humic_orgP: Concentration of humic organic P in the layer (mg/kg) 
"""

humic_orgP = 0.125*humic_orgN


"""
fresh_orgP_surf: the phosphorus in the fresh organic pool in the top 10mm (kg P/ha)  
rsd_surf: material in the residue pool for the top 10mm of soil (kg/ha)
"""


fresh_orgP_surf= 0.0003*rsd_surf


"""
#############################################################################
                        ###Mineral P###
#############################################################################
active_minP: Amount of phosphorus in the active mineral pool (mg/kg)
solP: Amount of phosphorus in solution (mg/kg)
PAI: phosphorus availability index
stable_minP: Concentration of phosphorus in the stable mineral pool (mg/kg) 

"""

active_minP = solP*(1-PAI/PAI)
stable_minP = 4*active_minP # This is true when both the pools are in equilibrium


"""
While  SWAT  allows  nutrient  levels  to  be  input  as  concentrations,  it 
performs all calculations on a mass basis.
massP: 
concP: concentration of phosphorus in a layer (mg/kg or ppm), 
b: the bulk density of the layer (Mg/m3)
lyr_dpth: the depth of the layer (mm).  
"""
massP = concP * b * lyr_dpth / 100 


#############################################################################
            #MINERALIZATION & DECOMPOSITION/IMMOBILIZATION
#############################################################################

"""
Two  sources  are considered  for  mineralization:  the  fresh  organic  
P  pool  associated  with  crop residue and microbial biomass and the active
organic P pool associated with soil humus.  Mineralization  and  decomposition 
are  allowed  to  occur  only  if  the temperature of the soil layer
is above 0 degreeC. Mineralization and decomposition are dependent on water 
availability and temperature.  Two  factors  are  used  in  the  mineralization
and  decomposition equations to account for the impact of temperature and water
on these processes. 

NCTF: nutrient cycling temperature factor for each layer [not allowed to be smaller than 0.1]
NCWF: nutrient cycling water factor for each layer [not allowed to be smaller than 0.05]
soilT: temperature  of each layer in degreeC 
wc: soil water content for a given layer on a given day (mm)
fc: water content of a given layer at field capacity on a given day (mm)
"""

NCTF = 0.9 * soilT/ soilT + exp(9.93-0.312*soilT) + 0.1

NCWF = wc/fc


"""
HUMUS MINERALIZATION

active_orgP:  amount of phosphorus in the active organic pool (kg P/ha)
stable_orgP:  amount of phosphorus in the stable organic pool (kg P/ha)
minP_humicorgP: the phosphorus mineralized from the humus active organic P pool (kg  P/ha)
B_min: rate  coefficient  for  mineralization  of  the  humus  active organic nutrients
Phosphorus mineralized from the humus active organic pool is added to the
solution P pool in the layer
"""

active_orgP = humic_orgP * (active_orgN/active_orgN + stable_orgN)
stable_orgP = humic_orgP * (stable_orgN/active_orgN + stable_orgN)

minP_humicorgP = 1.4* B_min * np.sqrt(NCTF*NCWF) * active_orgP



"""
Decomposition and mineralization 

Allowed only in first soil layer and controlled by a decay rate constant that 
is updated daily. 
The decay rate constant is  calculated  as  a  function  of  the  C:N  ratio 
and  C:P  ratio  of  the  residue, temperature and soil water content. 

E_cn: C:N ratio of the residue n the soil layer 
rsd: residue in layer ly (kg/ha)
0.58: fraction of residue that is carbon
fresh_orgN: nitrogen in the fresh organic pool in layer (kg N/ha)
NO3: amount of nitrate in layer (kg N/ha). 
E_cp: C:P ratio of the residue n the soil layer
solP : amount of phosphorus in solution in layer (kg P/ha)
fresh_orgP: hosphorus in the fresh organic pool in layer (kg P/ha).
"""

E_cn = 0.58 * rsd /fresh_orgN + NO3

E_cp = 0.58 * rsd /fresh_orgP + solP


"""
DECAY RATE CONSTANT

The decay rate constant defines the fraction of residue that is decomposed. 

d_rate_const: decay rate constant 

"""

d_rate_const =  B_rsd * NCRC * np.sqrt(NCTF*NCWF)


a = np.exp(-0.693*(E_cn - 25)/25)
b = np.exp(-0.693*(E_cp - 200)/200)
c = [a,b,1]

NCRC = min(c)


"""
minP_freshorgP : Mineralization from the residue fresh organic P pool 
decP_freshorgP : decomposition from the residue fresh organic P pool 
"""

minP_freshorgP = 0.8 * d_rate_const * fresh_orgP

decP_freshorgP = 0.2 * d_rate_const * fresh_orgP


#############################################################################
            #SORPTION OF INORGANIC P
#############################################################################

"""
SWAT assumes a rapid equilibrium exists between solution  P  and  an  “active”
mineral  pool. slow  reaction  is simulated  by  the  slow  equilibrium  
assumed  to  exist  between  the  “active”  and “stable”  mineral  pools. 

Equilibration between the solution and active mineral pool is governed by 
the phosphorus availability index (PAI)

solnP_f: amount  of phosphorus in solution after fertilization and incubation
solnP_i: amount  of phosphorus in solution before fertilization 
fert_minP: the amount of soluble P fertilizer added to the sample

"""

PAI = solnP_f - solnP_i/fert_minP

"""
MOVEMENT BETWEEN ACTIVE MINERAL POOL AND SOLUTION

The  movement  of  phosphorus  between  the  solution  and  active  mineral 
pools is governed by the equilibration equations: 

P_trans_sol_active_P: Amount of phosphorus transferred between the soluble and 
active mineral pool (kg P/ha). Positive value indicates transfer from 
solution to active mineral pool and vice a versa.

solP:  Amount of phosphorus in solution (kg P/ha)

active_minP: amount of phosphorus in the active mineral pool (kg P/ha)
    
"""
#solP =0.2
#active_minP = 0.02
#PAI = 0.6
#minP=0.4

if solP > minP * (PAI/(1-PAI)):
    P_trans_sol_active_P = 0.1*(solP - active_minP * (PAI/(1-PAI)))
    print P_trans_sol_active_P
elif solP < minP * (PAI/(1-PAI)):
    P_trans_sol_active_P = 0.6*(solP - active_minP * (PAI/(1-PAI)))
    print P_trans_sol_active_P


"""
MOVEMENT BETWEEN STABLE MINERAL POOL AND SOLUTION

When not in equilibrium, the movement of phosphorus between the active 
and stable mineral pools is governed by the equations: 
"""

if stable_minP < 4 * active_minP:
    P_trans_sol_stable_P = B_eqp * (4 * active_minP - stable_minP)
elif stable_minP > 4 * active_minP:
    P_trans_sol_stable_P = 0.1 * B_eqp * (4 * active_minP - stable_minP)


#############################################################################
            #LEACHING
#############################################################################

"""
Due to the low mobility of phosphorus, SWAT allows soluble P to leach 
only from the top 10 mm of soil into the first soil layer.

P_perc : Amount of phosphorus moving from the top 10 mm into the first 
soil layer (kg P/ha)

solP_10mm : amount of phosphorus in solution in the top 10 mm (kg P/ha)
W_prec_surf:  amount of water percolating to the first soil layer 
from the top 10 mm on a given day (mm)
B_d: bulk density of the top 10mm  (Mg/m3) (assumed  to  be  equivalent  to  
bulk  density  of  first  soil  layer)
dpth: depth of the “surface” layer (10 mm)
Kd_perc: phosphorus percolation  coefficient  (m3/Mg)


The  phosphorus  percolation  coefficient  is  the ratio  of  the  phosphorus
concentration  in  the  surface  10  mm  of  soil  to  the concentration of 
phosphorus in percolate. 
"""


P_perc = (solP_10mm * W_prec_surf) / (10 * B_d * dpth * Kd_perc)









