# -*- coding: utf-8 -*-
"""
Created on Tue Apr 10 15:26:53 2018

@author: Chinmay
"""
import pandas as pd
import numpy as np


def read_wepp_dat(filepath):
    data = pd.read_excel(filepath)
    data['Percolation(mm)'] = (data['Percolation (m^3)'] /data['Area (m^2)'])*1000
    data['Runoff(mm)'] = data['Runoff (m^3)']/ data['Area (m^2)'] *1000
    data['Lateral;(mm)'] = data['Lateral (m^3)']/ data['Area (m^2)'] *1000
    return data

def conc_to_mass(conc, B_d, Toplyr_dpth):
    """converts conc (mg/kg, etc) to mass (Kg P /ha)
       Toplyr_dpth: depth of the layer (mm).
       B_d: layer bulk density (Mg/m3)"""
    mass = conc * B_d * Toplyr_dpth /100
    return mass


def initP(orgC, fresh_orgN, fresh_orgP, solP):
    """orgC: percent carbon
       fresh_orgN: fresh organic N (kg N/ha)
       fresh_orgP: fresh organic P (kg P/ha)
       solP: soluble P (mg P/ kg)"""
    #orgC = 50 ###################
    humic_orgN1 = 10**4 * (orgC/14)
    humic_orgP1 = 0.125*humic_orgN1
    humic_orgN = conc_to_mass(humic_orgN1, 1.3, 100)
    humic_orgP = conc_to_mass(humic_orgP1 , 1.3, 100)
    solP = conc_to_mass(solP,1.3,100)
    return humic_orgN, humic_orgP, fresh_orgN, fresh_orgP, solP


## Soil water characteristics (SWC)
##    fc: field capacity (mm)
##    wc: water content (mm)
##   B_d: Bulk Density (Mg/m3) 
#SWC = namedtuple('SWC', 'fc wc B_d')
#
## Soil water characteristics in the rooting zone.
##    wc: soil water content (in mm)
##    vwc: soil water content (volumetric)
##    critical: soil water content threshold, below which
##              plant water stress occurs
##    sat: saturation point
##    fc: field capacity
##    pwp: permanent wilting point
#RootZone = namedtuple('RootZone', 'wc vwc critical sat fc pwp')
#    
#class SoilLayers(object):
#    
#    LyrDpthAct = 0.0 # used to determine layer depth
#    
#    def __init__(self):
#        """Initialize the SoilLayer object."""
#        #thick - thickness of the soil layer (mm)
#        self.thick = 0.0 
#        #wc - water content (mm)
#        self.wc = 0.0
#        #accthick - cumulative thickness (mm)
#        self.accthick = 0.0
#        #depth - depth of layer from soil surface (mm)
#        self.depth = 0.0
#        #swc - soil water characteristics 
#        self.swc = SWC(0.0, 0.0, 0.0)
#        self.prev = None
#        self.next = None
#        