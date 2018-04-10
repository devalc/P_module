# -*- coding: utf-8 -*-
"""
Created on Tue Apr 03 12:28:10 2018
Updated on Mon Apr 09 10:36:28 2018

@author: chinmay deval


This file contains all the formulations in the P module of SWAT, as explained 
in 2009 SWAT theoretical description document. Other reference materials used 
include 
1. SWAT 2005 theoretical document,
2. "Phosphorus Modeling in Soil and Water Assessment Tool (SWAT) Model" chapter in
the book Modelling Phosphorus in the Environment by Radcliffe etal 2007 
and 
3. "DEVELOPMENT OF A COUPLED WEPP-WQ MODEL" chapter from Phd thesis 
"CLIMATE CHANGE IMPACTS ON SOIL EROSION AND NUTRIENT LOSSES IN THE 
GREAT LAKES REGION" by Wang, 2015.
4. http://graham.umich.edu/media/files/How_SWAT_models_P.pdf

 

"""

#############################################################################
                        #import modules 
#############################################################################

import pandas as pd
import numpy as np
#import matplotlib
import matplotlib.pyplot as plt



#############################################################################
                        # Read data
#############################################################################

data= pd.read_excel('C:/Chinmay/Github/P_module/dummy_data.xlsx')


#############################################################################
                    # convert values to mm 
#############################################################################

data['Percolation(mm)'] = (data['Percolation (m^3)'] /data['Area (m^2)'])*1000
data['Runoff(mm)'] = data['Runoff (m^3)']/ data['Area (m^2)'] *1000
data['Latera;(mm)'] = data['Lateral (m^3)']/ data['Area (m^2)'] *1000


#############################################################################
                    # Initialize soil
#############################################################################

#Bulk density values should fall between 1.1 and 1.9 Mg/m3.

#layer bulk density (Mg/m3)
B_d = 1.3 

#The WEPP model makes use of up to ten soil layers, the top two with thicknesses
#of 100 mm each and subsequent lower layers with thicknesses of 200 mm each 
#(Wang, 2015: coupled WEPP-WQ module thesis)


#lyr_dpth: the depth of the layer (mm).
lyr_dpth = 100

#wc: soil water content for a given layer on a given day (mm)
#wc = 

#fc: water content of a given layer at field capacity on a given day (mm)
#fc =

#soilT: temperature  of each layer in degreeC
# soilT =  


#############################################################################
                    #Initialize P
#############################################################################

#Organic phosphorus levels are assigned assuming that the N:P ratio for 
#humic materials is 8:1
#Users may define the concentration of  organic phosphorus (dry  weight  basis)
#contained  in  humic  substances  for  all  soil  layers  at  the  beginning  
#of  the  simulation.  If  the  user does  not  specify  initial  organic  P  
#concentrations,  initialize  levels  of  organic  phosphorus  using  the 
#following  equations (see Chapter  3:2 of SWAT the  Theoretical Documentation).

#humic_orgN1: Concentration of humic organic N in the layer (mg/kg) 
#humic_orgP1: Concentration of humic organic P in the layer (mg/kg) 
#orgC: amount of organic carbon in the layer (%)

            ### humic phosphorus
orgC = 50 ###################
humic_orgN1 = 10**4 * (orgC/14)
humic_orgP1 = 0.125*humic_orgN1
humic_orgP = humic_orgP1* B_d * lyr_dpth/100 #converts conc (mg/kg, etc) to mass (Kg P /ha)
humic_orgN = humic_orgN1* B_d * lyr_dpth/100 #converts conc (mg/kg, etc) to mass (Kg P /ha)


            ### fresh organic phosphorus
#P,N in fresh organic pool is set to 0 in all layers except top 10mm if soil.

#fresh_orgP: P in the fresh organic pool in layer (kg P/ha)
#fresh_orgN: N in the fresh organic pool in layer (kg N/ha)


fresh_orgP = 2 ##################
fresh_orgN = 20 ##################



            ###Soluble Phosphorus
#solP: Amount of phosphorus in solution (mg/kg)

#Users  may  define  the  concentration  of  solution  P  (dry weight  basis)  
#for  all  soil  layers  at  the  beginning  of  the  simulation.  If  the  user 
#does  not  specify  initial  solution  P concentrations,  SWAT  will  initialize
#the  concentration  to  5 mg P/kg soil in all soil layers for unmanaged land under
#native vegetation and 25 mg P kg−1 soil for cropland conditions (Neitsch et al. 2001).

solP1 = 5
solP = solP1 * B_d * lyr_dpth/100 #converts conc (mg/kg, etc) to mass (Kg P /ha)



            ### Nitrate concentration in soil layer
#NO3: amount of nitrate in layer (kg N/ha)

#Users may define the concentration of nitrate (dry weight basis) for all soil
#layers at the beginning of the simulation. If  the  user  does  not specify 
#initial  nitrate  concentrations, SWAT  initializes levels of  nitrate
#using  the  following  equations (see Chapter  3:2 of SWAT the 
#Theoretical Documentation)


NO3 = 7 * np.exp(-lyr_dpth/1000) 

            ### Init active and stable Organic Nitrogen

#active_orgN1: concentration  of  nitrogen  in  the  active  organic  pool (mg/kg) #### Kg/ha used
#stable_orgN1: concentration  of  nitrogen  in  the  stable  organic  pool (mg/kg) #### Kg/ha used
#Fr: The fraction of humic nitrogen in the active pool


Fr = 0.02 #default
active_orgN = humic_orgN * Fr
stable_orgN = humic_orgN * (1-Fr)


            ### Init active and stable Organic Phosphorus

#active_orgP:  amount of phosphorus in the active organic pool (kg P/ha)
#stable_orgP:  amount of phosphorus in the stable organic pool (kg P/ha)


active_orgP = humic_orgP * (active_orgN/active_orgN + stable_orgN)
stable_orgP = humic_orgP * (stable_orgN/active_orgN + stable_orgN)

#############################################################################
    # Transformation Factors and Rate Constants
#############################################################################

            ###Mineralization

#Two  sources  are considered  for  mineralization:  the  fresh  organic  
#P  pool  associated  with  crop residue and microbial biomass, and the active
#organic P pool associated with soil humus.  
#Mineralization  and  decomposition are  allowed  to  occur  only  if  
#the temperature of the soil layer is above 0 degreeC. 
#Mineralization and decomposition are dependent on water availability and 
#temperature. 
#Two  factors (NCTF and NCWF) are  used  in  the  mineralization and  
#decomposition equations to account for the impact of temperature and water on 
#these processes. 


#NCTF: nutrient cycling temperature factor for each layer [not allowed to be smaller than 0.1]
#NCWF: nutrient cycling water factor for each layer [not allowed to be smaller than 0.05]
#soilT: temperature  of each layer in degreeC 
#wc: soil water content for a given layer on a given day (mm)
#fc: water content of a given layer at field capacity on a given day (mm)


#NCTF = 0.9 * soilT/ soilT + np.exp(9.93-0.312*soilT) + 0.1
NCTF = 0.1 #########################
#NCWF = wc/fc
NCWF=0.05  #########################


            ###Decomposition

#Decomposition is allowed only in first soil layer 
#and controlled by a decay rate constant that is updated daily. 
#The decay rate constant is  calculated  as  a  function  of  the  C:N  ratio 
#and  C:P  ratio  of  the  residue, temperature and soil water content. 
#E_cn: C:N ratio of the residue n the soil layer 
#rsd: residue in layer ly (kg/ha)
#0.58: fraction of residue that is carbon
#fresh_orgN: nitrogen in the fresh organic pool in layer (kg N/ha)
#NO3: amount of nitrate in layer (kg N/ha). 
#E_cp: C:P ratio of the residue n the soil layer
#solP : amount of phosphorus in solution in layer (kg P/ha)
#fresh_orgP: hosphorus in the fresh organic pool in layer (kg P/ha).


rsd=0.2 ############################### (same as RSD_COVCO Residue cover factor for 
                              ######### computing fraction of cover (0.1 –0.5)) ?????

E_cn = 0.58 * rsd /fresh_orgN + NO3

E_cp = 0.58 * rsd /fresh_orgP + solP


            ###DECAY RATE CONSTANT
            
#The decay rate constant defines the fraction of residue that is decomposed. 


#d_rate_const: decay rate constant 
#B_rsd: The  fraction  of  residue  which  will  decompose  in  a  day assuming 
#optimal  moisture,  temperature,  C:N  ratio  and  C:P ratio (default=0.05).
# NCRC: nutrient cycling residue composition factor

a = np.exp(-0.693*(E_cn - 25)/25)
b = np.exp(-0.693*(E_cp - 200)/200)
c = [a,b,1]

NCRC = min(c)

B_rsd = 0.05

d_rate_const =  B_rsd * NCRC * np.sqrt(NCTF*NCWF)


#############################################################################
                    # Transformations
#############################################################################

            ### HUMUS MINERALIZATION

#minP_humicorgP: the phosphorus mineralized from the humus active organic P pool (kg  P/ha)
#B_min: rate  coefficient  for  mineralization  of  the  humus  active organic nutrients
#Phosphorus mineralized from the humus active organic pool is added to the
#solution P pool in the layer
           
B_min = 0.0003 ##############################

minP_humicorgP = 1.4* B_min * np.sqrt(NCTF*NCWF) * active_orgP



            ### FRESH ORGANIC P (PLANT RESIDUE) MINERALIZATION and DECAY 

#minP_freshorgP : Mineralization from the residue fresh organic P pool 
#decP_freshorgP : decomposition from the residue fresh organic P pool


minP_freshorgP = 0.8 * d_rate_const * fresh_orgP
decP_freshorgP = 0.2 * d_rate_const * fresh_orgP


            ###SORPTION OF INORGANIC P


#SWAT assumes a rapid equilibrium exists between solution  P  and  an  “active”
#mineral  pool. slow  reaction  is simulated  by  the  slow  equilibrium  
#assumed  to  exist  between  the  “active”  and “stable”  mineral  pools. 
#Equilibration between the solution and active mineral pool is governed by 
#the phosphorus availability index (PAI)


#solnP_f: amount  of phosphorus in solution after fertilization and incubation
#solnP_i: amount  of phosphorus in solution before fertilization 
#fert_minP: the amount of soluble P fertilizer added to the sample
#PAI is calulated as:

    #{NOTE_TO_SELF: This value, given the formulation, seems to be a determinant
      #of Ag soils}
    
#PAI = solnP_f - solnP_i/fert_minP 

#if the value is not provided, default value is set to 0.4"""
PAI = 0.4


        ###MOVEMENT BETWEEN ACTIVE MINERAL POOL AND SOLUTION


#The  movement  of  phosphorus  between  the  solution  and  active  mineral 
#pools is governed by the equilibration equations: 
#P_trans_sol_active_P: Amount of phosphorus transferred between the soluble and 
#active mineral pool (kg P/ha). Positive value indicates transfer from 
#solution to active mineral pool and vice a versa.
#solP:  Amount of phosphorus in solution (kg P/ha)
#active_minP: amount of phosphorus in the active mineral pool (kg P/ha)
#

active_minP = solP * (1-PAI/PAI)
stable_minP = 4 * active_minP

if solP > active_minP * (PAI/(1-PAI)):
    P_trans_sol_active_P = 0.1*(solP - active_minP * (PAI/(1-PAI)))
    #print P_trans_sol_active_P
elif solP < active_minP * (PAI/(1-PAI)):
    P_trans_sol_active_P = 0.6*(solP - active_minP * (PAI/(1-PAI)))
    #print P_trans_sol_active_P
    
    
        ###MOVEMENT BETWEEN STABLE AND ACTIVE MINERAL POOL
#the movement of phosphorus between the active and stable mineral pools is 
#governed by the equations:
#B_eqp:: slow equilibrium constant set to 0.0006 per day

B_eqp= 0.0006 
if stable_minP < 4 * active_minP:
    P_trans_sol_stable_P = B_eqp * (4 * active_minP - stable_minP)
elif stable_minP > 4 * active_minP:
    P_trans_sol_stable_P = 0.1 * B_eqp * (4 * active_minP - stable_minP)    
    

#############################################################################
                    # LOSSES
############################################################################# 

                    ###P LEACHING 


#P_perc : Amount of phosphorus moving from the top 10 mm into the first 
#soil layer (kg P/ha)
#solP_10mm : amount of phosphorus in solution in the top 10 mm (kg P/ha)
#W_prec_surf:  amount of water percolating to the first soil layer 
#from the top 10 mm on a given day (mm)
#B_d: bulk density of the top 10mm  (Mg/m3) (assumed  to  be  equivalent  to  
#bulk  density  of  first  soil  layer) Bulk density values should fall 
#between 1.1 and 1.9 Mg/m3. ]
#dpth: depth of the “surface” layer (10 mm)
#Kd_perc: phosphorus percolation  coefficient  (m3/Mg)


W_perc_surf = data['Percolation(mm)']
#solP_10mm = 1
#B_d = 1.3 
#dpth = 200
Kd_perc = 10 #The value can range from 10.0 to 17.5. default is 10

P_perc = (solP * W_perc_surf) / (10 * B_d * lyr_dpth * Kd_perc)

data['P_percolation(Kg P /ha)'] = P_perc



            ### Phosphorus Uptake by Plants
    
#The  model assumes that plant uptake of P comes from the labile P pool (ie. solP) 
#P_uptake: plant P demand (kg /ha)[ Potential P uptake]
#Bio_optp:  expected  amount  of  P  content  in  plant  biomass  at a  given  
#plant  stage
#Bio_p: actual  amount  of  P  content  in  plant biomass

Bio_optp = 2.5 ####################RANDOM
Bio_p = 1.7 #######################RANDOM

P_uptake = 1.5 *(Bio_optp - Bio_p) 

         ### Phosphorus  uptake  from different soil lyr depths 


#P_uptake: potential P uptake by the plant to soil depth (Kg/ha)
#lyr_dpth: soil  depth  from  the  surface  (mm)
#lyr_dpth_r: rooting depth (mm)
#B_p: distribution parameter for P uptake default set to 20


B_p = 20.0
lyr_dpth_r = 20

P_uptake = (P_uptake / 1-np.exp(-B_p))*(1-np.exp(-B_p*lyr_dpth/lyr_dpth_r))
P_uptake


#The  P uptake  for  a  soil  layer  is  calculated  as  a  difference
#between  P  uptake  at  the  lower and upper boundary of that soil layer.

#P_actual: actual amount of P removed from soil #actUptake
#P_uptake: 
#P_demand: P  uptake  demand  not  met  by  overlying  soil  layers  (kg  P  /ha) # potUptake
#P_sol: amount  of  labile  P  present  in  the  soil  (kg  P / ha). 


P_demand = 1.1
P_actual = min((P_uptake + P_demand), solP)
P_actual



#P_stress: P stress for a given day.If a sufficient amount of P is not available
#in the soil for optimum plant growth, plants may experience P stress.
#
#P_stress varies non-linearly between 0, optimal P content and 1 when P content
#of the plant is <= 50% of optimal value 
#
#
#phi_P: scaling factor
#Bio_P: actual  P  content  of  plant  biomass (kg P /ha)
#Bio_optp: optimum P content of plant biomass (kg P /ha)


phi_P = 200 * (Bio_p /Bio_optp -0.5)
P_stress = 1 - phi_P /(phi_P + np.exp(3.535-0.02597*phi_P))


        ###PHOSPHORUS MOVEMENT IN SURFACE RUNOFF
        

#P_surf: the  amount  of  soluble  P  transported  by  surface  runoff  (kg  P  ha−1)
#solP: the amount of labile P (p in solution) in the top 10 mm (kg P /ha)
#Q_surf: the amount of surface runoff on a given day as depth (mm)
#B_d: bulk density of the top 10 mm of the soil (equivalent to B_d of first soil layer)
#dpth: depth of surface soil layer (mm)
#Kd_surf: phosphorus soil partitioning coefficient (m3/mg)


Kd_surf = 175.0 ### default value
Q_surf= data['Runoff(mm)']

P_surf = solP * Q_surf/ B_d * lyr_dpth * Kd_surf

data['SolubleP(Kg P /ha)'] = P_surf

        ###PHOSPHORUS TRANSPORTED BY SEDIMENT
        
#mass  of  P  transported  with  sediment  to  the  stream  
#uses a loading function developed by McElroy et al. (1976) and Williams and 
#Hann (1978).

#SedP_surf:amount of P transported with sediment to the main channel in surface
#runoff (kg P/ ha)
#SedP_conc: concentration of P attached to sediment in top 
#10 mm (g P metric ton soil−1)
#sed: sediment yield on a given day (metric tons)
#area: HRU area (ha) (OFE in this case? or outlet for that matter, if we want 
#values at the wshed outlet)
#E_ps: P enrichment ratio: The ratio of the concentration of P transported with 
#the sediment to the 
#concentration of P in the soil surface layer is  defined  as  the  P  
#enrichment  ratio 
#sed_conc_sq: concentration  of  sediment  in  surface  runoff  (mg  sediment /m3)


sed = data['Sediment (tons)']
area= data['Area (m^2)'] * 1e-4

# Enrichment ratio is calculated for each storm event using formula by (Menzel 1980)
sed_conc_sq = sed / 10 * area * Q_surf
E_ps = 0.78 * (sed_conc_sq)**(- 0.2468)

sedP_conc = 100*(active_minP + stable_minP + humic_orgP + fresh_orgP/B_d * lyr_dpth)
SedP_surf = 0.001* sedP_conc *(sed/area) * E_ps

data['Sediment P (Kg P /ha)'] = SedP_surf


###############################################################################################################
#
#
#
#
#
##############################################################################
#        #functions to display first and last n rows in dataframe
#"""just for evaluation purpose, remove this block once done"""
##############################################################################
#def front(self, n):
#    return self.iloc[:, :n]
#
#def back(self, n):
#    return self.iloc[:, -n:]
#
#pd.DataFrame.front = front
#pd.DataFrame.back = back
#data.front(50).back(50)
#
##############################################################################
#                                # Some plots
##############################################################################
#
#plt.plot( 'Date', 'P_percolation(Kg P /ha)', data=data, color='red', label = 'Perc_P')
#plt.xlabel('Time')
#plt.ylabel('Phosphorus [Kg P/ha]')
#plt.legend(prop={'size': 8})
#
#plt.plot( 'Date', 'SolubleP(Kg P /ha)', data=data, color='blue', label = 'SolP')
#plt.xlabel('Time')
#plt.ylabel('Phosphorus [Kg P/ha]')
#plt.legend(prop={'size': 8})
#
#plt.plot( 'Date', 'Sediment (tons)', data=data, color='green', label = 'Sed')
#plt.plot( 'Date', 'Sediment P (Kg P /ha)', data=data, color='blue', label = 'SedP')
#plt.xlabel('Time')
#plt.ylabel('Phosphorus [Kg P/ha]')
#plt.legend(prop={'size': 8}) 
#
#plt.plot( 'Date', 'Runoff (m^3)', data=data, color='blue', label = 'Q_surf')
#plt.show()
#
