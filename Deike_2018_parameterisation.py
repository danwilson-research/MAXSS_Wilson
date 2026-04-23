# -*- coding: utf-8 -*-
"""
Created on Wed Jun  1 11:06:01 2022

This is an implementation of Deike et.al 2018.

Deike, L., & Melville, W. K. (2018).
Gas transfer by breaking waves,
Geophysical Research Letters, 45, 10,482–10,492.
https://doi.org/10.1029/2018GL078758

https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2018GL078758?casa_token=o10woxAFDCYAAAAA%3ADNcM2uCmE9CzDRnesz12Y7j_x7_R472sk9TzEAbcQLQewMnAmAMlA1kFwd78PyYMpUjYyrVtgEIwwww

Important note:this parameterisation is a parameterisation for higher solubility gases like CO2 and DMS.
Use with caution for less soluble gases.

Important note: This parameterisation raises friction velocity and significant
wave height to fractional powers e.g. 5/3 and 4/3.
This means that to capture the variability in these two parameters as we do by taking the
second moment of the wind e.g <u10>^2 for other parameterisations. It is necessary 
to calculate the moment to this power for both these variables beforehand 
and you must pass that as an input in the configuration file.
e.g. <u*>^1.666 <Hs>^1.333. Here these variables are crudely called 
friction_velocity_1_666 and sig_wv_ht_1_333

@author: rps207
"""

from fluxengine.core.rate_parameterisation import KCalculationBase
from fluxengine.core.datalayer import DataLayer

#Deike 2018
class k_Deike2018(KCalculationBase):
    def __init__(self):
        self.name = self.__class__.__name__;
    
    def input_names(self):
        return ["friction_velocity", "sig_wv_ht","solubility_skin", "scskin"];
        
    def output_names(self):
        return ["k"];
    
    def __call__(self, data):
        function = "(rate_parameterisation.py: k_Deike2018.__call__)";
        print("%s Using the Deike et al., 2018 (D18) k parameterisation" % (function))
        try:
            #for ease of access, simply assign attributes to each input/output.
            for name in self.input_names() + self.output_names():
                setattr(self, name, data[name].fdata);
            data["k"].standardName="gas_transfer_velocity_of_carbon_dioxide" ;
            data["k"].longName="Deike et al., 2018 (D18) gas transfer velocity";
        except KeyError as e:
            print("%s: Required data layer for selected k parameterisation was not found." % function);
            print(type(e), e.args);
            return False;
        
        #determine the Deike et al., 2018  k relationship
        for i in arange(len(self.k)):   
            self.k[i] = DataLayer.missing_value
            if ( (self.friction_velocity[i] != DataLayer.missing_value) and (self.sig_wv_ht[i] != DataLayer.missing_value) and (self.solubility_skin[i] != DataLayer.missing_value)and (self.scskin[i] != DataLayer.missing_value) and (self.scskin[i] > 0.0) ):
                
                # This is the implementation of Deike et.al 2018 (their equation 12)
                self.A_nb[i]=1.55e-4 # Deike note this constant is the non breaking constant from Fairall et al., 2011
                self.A_b[i]=1e-5
                g=9.81
                gamma=5/3
                zeta=4/3
                
                #this is the equation for kw 
                #- note to get the right side of their equation 12 in terms of kw divide through by  sqrt(self.scskin[i]/660.0)          
                #self.k[i] =(self.A_nb[i]*self.friction_velocity[i]* sqrt(660.0/self.scskin[i]))+((self.A_b[i]/self.solubility_skin[i])*((pow(self.friction_velocity[i],gamma))*(pow(pow(g*self.sig_wv_ht[i],0.5),zeta))))
                
                #this is the equation for k660- the right side of equation 12
                self.k[i] =(self.A_nb[i]*self.friction_velocity[i])+(((self.A_b[i]* sqrt(self.scskin[i]/660.0))/self.solubility_skin[i])*((pow(self.friction_velocity[i],gamma))*(pow(pow(g*self.sig_wv_ht[i],0.5),zeta))))

            else:
                self.k[i] = DataLayer.missing_value
        
        return True;


